################   Basic Multivariate Analysis   ###############
SA_uni <- function(physeq, cat, name){
	data <- data.frame(otuTable(physeq))
	samData <- data.frame(samData(physeq))
	tree <- tre(physeq)
	taxonomy <- data.frame(taxTab(physeq))

	ID <- sample.names(physeq)
	feature <- species.names(physeq)
	group <- samData[,cat]

	unidist <- UniFrac(physeq, weighted=TRUE, normalized=FALSE, parallel=TRUE, fast=TRUE)
	write.csv(as.matrix(unidist),paste(name,".dist.w.nn.uni.csv",sep=""))

	SA_basic(data, group, ID, name, uni=TRUE, unidist)
}

SA_basic <- function(data, group, ID, name, uni=FALSE, unidist=0){
	bcdist <- vegdist(data,"bray")
	write.csv(as.matrix(bcdist),paste(name,".dist.bc.csv",sep=""))

	##################   PCA, COA, NMDS   ###################
	library(ade4)
	library(rgl)
		pca <- dudi.pca(data, scannf=FALSE, nf=3)
		coa <- dudi.coa(data, scannf=FALSE, nf=3)
		nmds2.bc <- metaMDS(bcdist, k=2)
		nmds3.bc <- metaMDS(bcdist, k=3)

		hclusterplot(bcdist, group, paste(name,".hcluster.bc.png", sep=""))
		nmdsplot(nmds2.bc, group, paste(name,".nmds2.bc.png",sep=""), cex=1, label=ID)
		nmdsplot(nmds3.bc, group, paste(name,".nmds3.bc.png",sep=""), cex=1, label=ID)
#		Sys.sleep(10)

		if(uni==TRUE){
			nmds2.uni <- metaMDS(unidist, k=2)
			nmds3.uni <- metaMDS(unidist, k=3)
			hclusterplot(unidist,group,paste(name,".hcluster.uni.png",sep=""))
			nmdsplot(nmds2.uni,group, paste(name,".nmds2.uni.png",sep=""),cex=1, label=ID)
			nmdsplot(nmds3.uni,group, paste(name,".nmds3.uni.png",sep=""),cex=1, label=ID)
#			Sys.sleep(10)
		}

		pcplot(pca, group, paste(name,".pca.png",sep=""), cbox=0.5, label=ID)
#		Sys.sleep(10)
		pcplot(coa, group, paste(name,".coa.png",sep=""), cbox=0.5, label=ID)
#		Sys.sleep(10)

	detach("package:rgl")

	#########################################################
	########   Multi-test    #########
		#######  ANOSIM, MRPP, perMANOVA  #######
		anosim.bc <- anosim(bcdist, group, permutations=1000)
		mrpp.bc <- mrpp(bcdist, group, permutations=1000)
		per.bc <- adonis(bcdist~group, permutations=1000)
#		perMANOVA.bc <- adonis(data ~ group, method="bray", permutations=1000)
#		perMANOVA$aov.tab

		print(anosim.bc)
		print(mrpp.bc)
		print(per.bc)

		if(uni==TRUE){
			anosim.uni <- anosim(unidist, group, permutations=1000)
			mrpp.uni <- mrpp(unidist, group, permutations=1000)
			per.uni <- adonis(unidist~group, permutations=1000)

			print(anosim.uni)
			print(mrpp.uni)
			print(per.uni)
		}

	#########################################################
	##################  Diversity Index  ####################
	library(BiodiversityR)
		shannon <- diversity(data, index="shannon", MARGIN=1, base=exp(1))
		simpson <- diversity(data, index="simpson", MARGIN=1, base=exp(1))
#		k <- sample(nrow(data), 9)
#		rendiv <- renyi(data[k,])
#		plot(rendiv)

		richness<- estimateR(round(data))
#		Srar <- rarefy(Data_C, ceiling(min(rowSums(Data_C))))
#		sac <- specaccum(Data)
#		plot(sac, ci.type="polygon", ci.col="yellow")

		eco <- t(rbind(ID=ID, group=as.vector(group), shannon, simpson, richness))
		write.csv(eco, paste(name,".diver.csv",sep=""))
	detach("package:BiodiversityR")

	#########################################################
	##################         LDA       ####################
	library(lattice)
	library(rgl)
	library(MASS)
		lda <- lda(group~., data, CV=FALSE)
		lda.hat <- predict(lda,decision.values=TRUE)
		tabTrain <- confusion(group, lda.hat$class)
		denplot(lda.hat, group, paste(name,".lda.density.png",sep=""))
		ldaplot(lda.hat, group, paste(name,".lda.xyplot.png",sep=""))
#		Sys.sleep(10)
	detach("package:lattice")
	detach("package:rgl")

	#########################################################
	##################  Lin's CCC heatmap  ##################
#	library(ClassDiscovery)
#	library(ClassComparison)
#	library(epiR)
#	library(proxy)
#		lin <- function(x1, x2, ci="z-transform", conf.level=0.95){
#			tmp <- epi.ccc(x1, x2, ci, conf.level)
#			return(tmp$rho.c[,1])
#		}
#		pr_DB$set_entry(FUN=lin, names="Lin")
#		method="Lin"
#		heatmaplot(data, group, ID, method=method, paste(name,".",method,".heatmap.png",sep=""))
#	detach("package:ClassDiscovery")
#	detach("package:ClassComparison")
#	detach("package:epiR")
#	detach("package:proxy")
}

################   Basic Feature Selection   ###############
FS_uni <- function(physeq, cat, name, kernel="linear"){
	data <- data.frame(otuTable(physeq))
	samData <- data.frame(samData(physeq))
	tree <- tre(physeq)
	taxonomy <- data.frame(taxTab(physeq))

	ID <- sample.names(physeq)
	feature <- species.names(physeq)
	group <- samData[,cat]
	FS_basic(data, group, ID, name, kernel=kernel)
}

FS_basic <- function(data, group, ID, name, kernel="linear"){
	########################################################
	################  basic test selection  ################
	library(coin)
		parattest <- para_t.test(data,group,paired=FALSE,alternative="two.sided")
		parachitest <- para_chi.test(data,group,paired=FALSE,alternative="two.sided")
		paranpttest <- para_npt.test(data,group,paired=FALSE,alternative="two.sided")
		write.csv(cbind(parattest,parachitest,paranpttest),paste(name,".fs.Test.pairwise.csv",sep=""))

		if(levels(group)>2){
			aov <- sapply(c(1:ncol(data)),function(x){summary(aov(data[,x]~group))[[1]][["Pr(>F)"]][1]})
#			fisher <- sapply(c(1:ncol(data)),function(x){fisher.test(table(data[,x]==0,group),alternative=TRUE)$p.value})
			fisher <- sapply(c(1:ncol(data)),function(x){
					t<-table(data[,x]==0,group);
					p<-1;
					if(nrow(t)>1){p<-fisher.test(t,alternative=TRUE)$p.value}
					return(p);})
			tukey <- sapply(c(1:ncol(data)),function(x){TukeyHSD(aov(data[,x]~group))$group[,4]})
			colnames(tukey) <- colnames(data)
			nptukey <- para_Mann(data,group,conf.level=0.95)
			write.csv(t(rbind(fisher,aov,tukey,t(nptukey))),paste(name,".fs.Test.group.csv",sep=""))
		}
	detach("package:coin")

	################  Random Forest Selection ############
	library(randomForest)
	library(Boruta)
		set.seed(125)
		boruta <- Boruta(group~., data=cbind(group,data), doTrace=2, maxRuns=12, ntree=500) #  Default maxRuns=4
		borutaplot(boruta,paste(name, ".fs.rf.png",sep=""))
		tmp <- TentativeRoughFix(boruta)
		write.csv(data.frame(boruta$finalDecision,tmp$finalDecision),paste(name, ".fs.rf.csv",sep=""))
#		getConfirmedFormula(clean.Boruta)
#		getNonRejectedFormula(clean.Boruta)
	detach("package:Boruta")
	detach("package:randomForest")

	################  Indicator Species Analysis  #############
	library(labdsv)
		isa <- indval(data, group, numitr=1000)
#		summary(isa)
		summary(isa, p=0.05, digits=2, show=p)
		isa.core <- t(rbind(isa$indcls,isa$maxcls,isa$pval))[isa$pval<=0.05,]
		write.csv(isa.core, paste(name,".fs.isa.csv",sep=""))
	detach("package:labdsv")

	################   ELASTIC NET   ####################
	library(glmnet)
		elastic_mod <- paraelastic(data, group, family="multinomial")
		ENfsplot(elastic_mod,paste(name, ".fs.EN.png",sep=""))
		write.csv(as.matrix(elastic_mod[[5]]),paste(name,".fs.EN.coef.csv",sep=""))
		write.csv(elastic_mod[[6]],paste(name,".fs.EN.par.csv",sep=""))
	detach("package:glmnet")	

	################      HHSVM      ####################




	################  SVM-RFE selection  ################
	library(e1071)
	library(kernlab)
		ftable <- para_svmrfe(data, group, ID, kernel="linear")
		write.csv(ftable,paste(name, ".fs.svmrfe.",kernel,".csv",sep=""))
	detach("package:e1071")
	detach("package:kernlab")

	################   Return part of the model   #############
	return(elastic_mod)
}


CL_uni <- function(physeq, cat, name, kernel="linear"){
	data <- data.frame(otuTable(physeq))
	samData <- data.frame(samData(physeq))
	tree <- tre(physeq)
	taxonomy <- data.frame(taxTab(physeq))

	ID <- sample.names(physeq)
	feature <- species.names(physeq)
	group <- samData[,cat]
	CL_basic(data, group, ID, name, kernel=kernel)
}

CL_basic <- function(data, group, ID, name, kernel="linear"){
	##########  SVM  ##########
	library(e1071)
		svm <- svm(group~., data, cross=nrow(data), kernel=kernel)
		svm.hat <- predict(svm,decision.values=TRUE)
		tabTrain <- confusion(group, svm.hat)
		print(summary(svm))
	detach("package:e1071")
	##########  RF  ##########
	library(randomForest)
		rf <- randomForest(group~., data=data, proximity=TRUE)
		print(rf)
		rf.cmd <- cmdscale(1-rf$proximity)
		col=c("blue","red","green","purple")
		plot(rf.cmd, col=c("blue","red","green","purple")[as.numeric(group)], ann=FALSE)
	detach("package:randomForest")
	################   Logestic Regression ##############
		glm.r <- glm(group~.,cbind(group,data), family=binomial())

}

################   Basic Linear Regression   ###############
LR_basic <- function(data, y, ID, name){
}

################   Shuffle Data   ###############
shuffleData <- function(Data, depth, permu=1000){
	plotdata <- matrix(0,nrow(Data),ncol(Data))
	Bacname <- rownames(Data)
	for(i in 1:ncol(Data)){	plotdata[,i] <- as.numeric(as.vector(Data[,i]))	}
	rownames(plotdata) <- rownames(Data)
	colnames(plotdata) <- colnames(Data)

	result <- list()
	for(i in 1:ncol(plotdata)){
		weight <- as.vector(plotdata[,i])
		ss <- table(as.factor(plotdata[,i]))
		if(sum(weight)>depth){
			ss <- replicate(permu, sample(rep(Bacname[weight>0],weight[weight>0]),depth))
#			ss <- replicate(permu, sample(Bacname[weight>0],depth,prob=weight)
			ss <- factor(as.vector(ss),levels=Bacname)
			result[[i]] <- table(ss)
		}
		else{ result[[i]] <- plotdata[,i]*permu	}
	}

	final <- do.call(cbind,result)/permu
	colnames(final) <- colnames(Data)
	final <- final[rowSums(final)>0,]
	return(final)
}

#######   Para-T-test   ######
para_t.test <- function(Data, cat, mu=0, paired=FALSE, alternative="two.sided", var.equal=FALSE, conf.level=0.95){
	col <- combn(levels(cat),2)
	ncat <- length(levels(cat))
	if(length(paired)==1){if(paired==FALSE){
		paired <- rep(FALSE,ncol(col))
	}}
	result <- matrix(0,ncol(Data),(ncat+ncol(col)))
	colnm <- paste(levels(cat),"mean",sep=" ")
	for(i in 1:ncol(col)){
		Tdata1 <- Data[cat==col[1,i],]
		Tdata2 <- Data[cat==col[2,i],]
		colnm <- c(colnm, paste(col[1,i],col[2,i],sep="<->"))
		for(j in 1:ncol(Data)){
			if(i < ncat){
				result[j,1] <- mean(Tdata1[,j])
				result[j,i+1] <- mean(Tdata2[,j])
			}
			if(	(length(Tdata1[,j])-1)*(length(Tdata2[,j])-1)	){
				Tresult <- t.test(Tdata1[,j],Tdata2[,j],mu=mu,alternative=alternative,paired=paired[i],var.equal=var.equal,conf.level=0.95)
				result[j,(ncat+i)] <- Tresult$p.value
			}
			else{	result[j,(ncat+i)] <- NA	}
		}
	}
	rownames(result) <- colnames(Data)
	colnames(result) <- colnm
	return(result)
}

#######   Para-Chisquare-test   ######
para_chi.test <- function(Data, cat, paired=FALSE, alternative="two.sided", conf.level=0.95){
	col <- combn(levels(cat),2)
	ncat <- length(levels(cat))
	if(length(paired)==1){if(paired==FALSE){
		paired <- rep(FALSE,ncol(col))
	}}
	result <- matrix(0,ncol(Data),(ncat+ncol(col)))
	colnm <- paste(levels(cat),"(existance)",sep="")
	for(i in 1:ncol(col)){
		Cdata1 <- Data[cat==col[1,i],]
		Cdata2 <- Data[cat==col[2,i],]
		colnm <- c(colnm, paste(col[1,i],col[2,i],sep="<->"))
		for(j in 1:ncol(Data)){
			r1 <- Cdata1[,j]
			r2 <- Cdata2[,j]
			a <- length(r1[r1!=0])
			b <- length(r2[r2!=0])
			c <- length(r1)-a
			d <- length(r2)-b
			if(i < ncat){
				result[j,1] <- paste(a,length(r1),sep=":")
				result[j,i+1] <- paste(b,length(r2),sep=":")
			}
			if(paired[i]==TRUE){
				result[j,(ncat+i)] <- mcnemar.test(matrix(c(a,b,c,d),2,2),correct=TRUE)$p.value
			}
			else{
				result[j,(ncat+i)] <- fisher.test(matrix(c(a,b,c,d),2,2),alternative=alternative)$p.value
			}
		}
	}
	rownames(result) <- colnames(Data)
	colnames(result) <- colnm
	return(result)
}

#######   Para-Non-parametric-T-test   ######
para_npt.test <- function(Data, cat, mu=0, paired=FALSE, alternative="two.sided", conf.level=0.95){
	col <- combn(levels(cat),2)
	ncat <- length(levels(cat))
	if(length(paired)==1){if(paired==FALSE){
		paired <- rep(FALSE,ncol(col))
	}}
	result <- matrix(0,ncol(Data),ncol(col))
	colnm <- c()
	for(i in 1:ncol(col)){
		NPTdata1 <- Data[cat==col[1,i],]
		NPTdata2 <- Data[cat==col[2,i],]
		colnm <- c(colnm, paste(col[1,i],col[2,i],sep="<->"))
		for(j in 1:ncol(Data)){
			NPTresult <- wilcox.test(NPTdata1[,j],NPTdata2[,j],mu=mu,alternative=alternative,paired=paired[i],conf.level=0.95,correct=TRUE)
			result[j,i] <- NPTresult$p.value
		}
	}
	rownames(result) <- colnames(Data)
	colnames(result) <- colnm
	return(result)
}

#######   Para-Kolmogorov-test   ######
para_ks.test <- function(Data, cat, mu=0, paired=FALSE, alternative="two.sided", conf.level=0.95){
}

#######   Para-Mann's NP multiple range test   #########
para_Mann <- function(alldata, group, conf.level=0.95){
	cat <- levels(group)
	grc <- combn(cat,2)
	result <- matrix(0, ncol(alldata), ncol(grc)+2)
	rownames(result) <- colnames(alldata)
	colnames(result) <- c("kruskal","NDWD",paste(grc[1,],grc[2,],sep=" - "))
	for(i in 1:ncol(alldata)){
		subdata <- data.frame(feat=alldata[,i],site=group)
		kw <- kruskal_test(feat~group, data=subdata, distribution=approximate(B=9999))
		result[i,1] <- pvalue(kw)[1]
		if(require("multcomp")){
			NDWD <- oneway_test(feat~group, data=subdata,
				ytrafo=function(data) trafo(data, numeric_trafo=rank),
				xtrafo=function(data) trafo(data, factor_trafo=function(x) model.matrix(~x-1) %*% t(contrMat(table(x),"Tukey"))),
				teststat="max", distribution=approximate(B=90000))
			### global p-value
			result[i,2] <- pvalue(NDWD)[1]
			result[i,3:ncol(result)] <- as.vector(pvalue(NDWD, method = "single-step"))
		}
	}
	return(result)
}

para_svmrfe <- function(data, group, ID, kernel="linear"){
	feature <- list()
print(kernel)
	for(n in 1:100){
		traindata <- data.frame()
		testdata <- data.frame()
		cat <- levels(group)
		for(i in 1:length(cat)){
			c_data <- data[group==cat[i],]
			c_train <- sample(1:nrow(c_data),5)
			c_test <- setdiff(1:nrow(c_data),c_train)
			traindata <- rbind(traindata, data.frame(group=rep(cat[i],5),c_data[c_train,]))
			testdata <- rbind(testdata, data.frame(group=rep(cat[i],length(c_test)),c_data[c_test,]))
		}
		g_train <- traindata[,1]
		d_train <- traindata[,-1]
		d_train <- d_train[,colSums(d_train)!=0]
		traindata <- data.frame(group=g_train,d_train)
		testdata <- testdata[,colnames(traindata)]
		levels(testdata[,1]) <- levels(traindata[,1])

		svmrfe <- svmRFE(group~., traindata, testdata, kernel=kernel)
#		svmRFEplot(svmrfe,filename="svmrfe.1.png")
#		write.csv(svmrfe$featureRank,"svmrfe.1.csv")
		frank <- svmrfe$featureRank
#		index <- sort((frank[,3]+frank[,4]),index.return=TRUE)$ix
		index <- sort((frank[,2]+frank[,4]),index.return=TRUE)$ix
#		index <- sort((frank[,2]+frank[,3]+frank[,4]),index.return=TRUE)$ix

		feature[[n]] <- rownames(frank)[1:index[length(index)]]
	}
	return(table(unlist(feature)))
}

#######   svmRFE feature selection  #######
svmRFE <- function(formula, train, test, kernel="linear"){
	fstr <- toString(formula)
	type <- strsplit(fstr,"\\, ")[[1]][2]

	featureRankedList <- matrix(0,(ncol(train)-1),7)
	rownames(featureRankedList) <- c(1:(ncol(train)-1))
	colnames(featureRankedList) <- c("F_rank","Classification_Acc","Train_Acc","Test_Acc","Opt_Cost","Opt_Gamma","nSV")
	featuresort <- rep(0,(ncol(train)-1))

	rankedFeatureIndex <- ncol(train)-1
	survivingFeatures <- colnames(train)
	survivingFeatures <- c(survivingFeatures[survivingFeatures!=type],type)
	result <- list()

	while(length(survivingFeatures)>1){
		# Initialize Data
		train_data <- train[,survivingFeatures]
		test_data <- test[,survivingFeatures]
		if(kernel == "radial"){
			# Grid Search find a coarse Grid
			obj <- gridsearch(formula, train=train_data, test=test_data, cross=nrow(train_data), 
					ranges=list(gamma=2^(-15:3),cost=2^(-5:15)),kernel="radial")
			coarse_Gamma <- log2(obj$best.model$gamma)
			coarse_C <- log2(obj$best.model$cost)
			range_Gamma <- seq(coarse_Gamma-2,coarse_Gamma+2,by=0.25)
			range_C <- seq(coarse_C-2,coarse_C+2,by=0.25)

			# Grid Search find best parameters
			obj_loo <- gridsearch(formula, train=train_data, test=test_data, cross=nrow(train_data), 
					ranges=list(gamma=2^range_Gamma,cost=2^range_C), kernel="radial")

			Gamma_loo <- obj_loo$best.model$gamma
			C_loo <- obj_loo$best.model$cost
print(obj_loo$best.performance)
			svmModel_loo <- obj_loo$best.model
			svmModel_loo <- svm(formula, data=train_data, cross=nrow(train_data), type="C-classification", 
					kernel=kernel, cost=C_loo, gamma=Gamma_loo)
		}

		if(kernel == "linear"){
			svmModel_loo <- svm(formula, train_data, cross=nrow(train_data), kernel=kernel)
#cachesize=500, 
		}

		# compute ranking criteria
		rankingCriteria <- svmweights(svmModel_loo)
		ranking <- sort(rankingCriteria, index.return=TRUE)$ix
#		result$model[[rankedFeatureIndex]] <- svmModel_loo
#		result$rankC[[rankedFeatureIndex]] <- rankingCriteria

		svmModel_loo_hat <- predict(svmModel_loo, train_data, decision.values=TRUE)
		svmModel_test_hat <- predict(svmModel_loo, test_data, decision.values=TRUE)
		tabTrain_loo <- confusion(train_data[,type], svmModel_loo_hat, printit=FALSE)
		tabTrain_test <- confusion(test_data[,type], svmModel_test_hat, printit=FALSE)

		# update feature ranked list
		featureRankedList[rankedFeatureIndex,1] <- rankedFeatureIndex
		featureRankedList[rankedFeatureIndex,2] <- svmModel_loo$tot.acc/100
#		featureRankedList[rankedFeatureIndex,2] <- obj_loo$best.performance
		featureRankedList[rankedFeatureIndex,3] <- tabTrain_loo$acc
		featureRankedList[rankedFeatureIndex,4] <- tabTrain_test$acc
		featureRankedList[rankedFeatureIndex,5] <- svmModel_loo$cost
		featureRankedList[rankedFeatureIndex,6] <- svmModel_loo$gamma
		featureRankedList[rankedFeatureIndex,7] <- sum(svmModel_loo$nSV)
		featuresort[rankedFeatureIndex] <- survivingFeatures[ranking[1]]

		rankedFeatureIndex <- rankedFeatureIndex-1
		# eliminate the feature with smallest ranking criterion
		survivingFeatures <- survivingFeatures[-ranking[1]]
	}
	rownames(featureRankedList) <- featuresort
	result$featureRank <- featureRankedList
	
	return(result)
}

################################################
# weights and Criteria of the hiperplane
################################################
svmweights <- function(model){
	rankingCriteria <- rep(0,ncol(model$SV))

#	rbf <- function(u, v, gamma){
#		exp(-gamma*sum((u-v)^2))
#	}
#	class(rbf) <- "kernel"
#	rbf <- rbfdot(sigma=model$gamma)

	if(model$nclasses==2){
		if(model$kernel==0){
			w <- t(model$coefs) %*% model$SV
			rankingCriteria <- w * w
		}
		if(model$kernel==1){
		}
		if(model$kernel==2){
			for(f in 1:ncol(model$SV)){
				KMat <- (model$coefs %*% t(model$coefs)) * kernelMatrix(rbfdot(sigma=model$gamma),model$SV[,-f],model$SV[,-f])
				rankingCriteria[f] <- -sum(KMat)
			}
		}
	}

	else{
		start <- c(1, cumsum(model$nSV)+1)
		start <- start[-length(start)]

		W <- matrix(0,ncol(model$SV),choose(model$nclasses,2))
		count <- 1
		for(i in 1:(model$nclasses-1)){
			for(j in (i+1):model$nclasses){
				## ranges for class i and j:
				ri <- start[i]:(start[i] + model$nSV[i] - 1)
				rj <- start[j]:(start[j] + model$nSV[j] - 1)
				## coefs and SV for (i,j):
				coefs <- c(model$coefs[ri, j-1], model$coefs[rj, i])
				SV <- data.matrix(model$SV[c(ri,rj),])
				if(model$kernel==0){
					w <- t(coefs) %*% SV
					W[,count] <- w * w
				}
				if(model$kernel==1){
				}
				if(model$kernel==2){
					for(nf in 1:ncol(model$SV)){
						KMat <- (coefs %*% t(coefs)) * kernelMatrix(rbfdot(sigma=model$gamma),SV[,-nf],SV[,-nf])
						W[nf,count] <- -sum(KMat)
					}
				}
				count <- count+1
			}
		}
		rankingCriteria <- rowMeans(W)
	}
	return(rankingCriteria)
}

gridsearch <- function(formula, train, test, cross=10, kernel="radial", ranges=list(gamma=2^(-15:3),cost=2^(-5:15))){
	fstr <- toString(formula)
	type <- strsplit(fstr,"\\, ")[[1]][2]

	para_grid <- expand.grid(ranges)
	para_list <- lapply(seq_len(nrow(para_grid)), function(i){para_grid[i,]})

	para_svm <- lapply(para_list, function(x){
			svm(formula, data=train, cross=cross, gamma=x$gamma, cost=x$cost, 
			cachesize=500, type="C-classification", kernel="radial")})

	para_train_predict <- lapply(para_svm, function(x){predict(x, train, decision.values=TRUE)})
	para_train_acc <- lapply(para_train_predict, function(x){confusion(train[,type], x, printit=FALSE)})

	para_test_predict <- lapply(para_svm, function(x){predict(x, test, decision.values=TRUE)})
	para_test_acc <- lapply(para_test_predict, function(x){confusion(test[,type], x, printit=FALSE)})

	performance <- data.frame(gamma=unlist(lapply(para_svm,function(x){x$gamma})),
				cost=unlist(lapply(para_svm,function(x){x$cost})),
				TrainError=unlist(lapply(para_svm,function(x){1-mean(x$acc)/100})),
				TrainDispersion=unlist(lapply(para_svm,function(x){sd(x$acc)/100})),
				TrainPredict=unlist(lapply(para_train_acc,function(x){1-x$acc})),
				TestPredict=unlist(lapply(para_test_acc,function(x){1-x$acc}))
			)

	performance$SumError <- performance$TrainPredict+performance$TestPredict
#	rank <- order(performance$SumError, performance$cost, performance$gamma)
	rank <- order(performance$TrainError, performance$cost, performance$gamma)
	result <- list(best.model=para_svm[[rank[1]]], best.performance=performance[rank[1],], performance=performance)
	return(result)
}

###### A function that calculates the confusion matrix and overall accuracy #####
confusion <- function(actual, predicted, names=NULL, printit=TRUE, prior=NULL){
	if(is.null(names)){	names <- levels(actual)	}
	result <- list()
	tab <- table(actual, predicted)
	acctab <- t(apply(tab, 1, function(x)x/sum(x)))
	dimnames(acctab) <- list(Actual=names,"Predicted (cv)"=names)
	result$tab <- acctab
	if(is.null(prior)){
		relnum <- table(actual)
		prior <- relnum/sum(relnum)
		acc <- sum(tab[row(tab)==col(tab)])/sum(tab)
		result$acc <- acc
	}
	else{
		acc <- sum(prior*diag(acctab))
		names(prior) <- names
		result$acc <- acc
	}
	if(printit){
		print(round(c("Overall accuracy"=acc,"Prior frequency"=prior),4))
	}
	if(printit){
		cat("\nConfusion matrix","\n")
		print(round(acctab,4))
	}
	return(result)
}

paraelastic <- function(data, group, family="multinomial"){
	alpha <- seq(0,1,0.01)
	paracvfit1 <- lapply(alpha, function(x){cv.glmnet(data.matrix(data),group,nfolds=5,alpha=x,family=family,type.measure="class")})
	paracvfit2 <- lapply(alpha, function(x){cv.glmnet(data.matrix(data),group,nfolds=5,alpha=x,family=family,type.measure="mse")})
	paracvfit1[[1]] <- NULL
	paracvfit2[[1]] <- NULL

	cvm1 <- matrix(unlist(lapply(paracvfit1,function(x){x$cvm})),nrow=(length(alpha)-1),byrow=TRUE)
	cvup1 <- matrix(unlist(lapply(paracvfit1,function(x){x$cvup})),nrow=(length(alpha)-1),byrow=TRUE)
	#cvsd1 <- matrix(unlist(lapply(paracvfit1,function(x){x$cvsd})),nrow=(length(alpha)-1),byrow=TRUE)
	#nzero1 <- matrix(unlist(lapply(paracvfit1,function(x){x$nzero})),nrow=(length(alpha)-1),byrow=TRUE)
	lamda1 <- matrix(unlist(lapply(paracvfit1,function(x){x$lambda})),nrow=(length(alpha)-1),byrow=TRUE)

	cvm2 <- matrix(unlist(lapply(paracvfit2,function(x){x$cvm})),nrow=(length(alpha)-1),byrow=TRUE)
	cvup2 <- matrix(unlist(lapply(paracvfit2,function(x){x$cvup})),nrow=(length(alpha)-1),byrow=TRUE)
	#cvsd2 <- matrix(unlist(lapply(paracvfit2,function(x){x$cvsd})),nrow=(length(alpha)-1),byrow=TRUE)
	#nzero2 <- matrix(unlist(lapply(paracvfit2,function(x){x$nzero})),nrow=(length(alpha)-1),byrow=TRUE)
	lamda2 <- matrix(unlist(lapply(paracvfit2,function(x){x$lambda})),nrow=(length(alpha)-1),byrow=TRUE)

	result1.cvm <- which(cvm1==min(cvm1), arr.ind=TRUE)
	result1.cvup<- which(cvup1==min(cvup1), arr.ind=TRUE)
	result2.cvm <- which(cvm2==min(cvm2), arr.ind=TRUE)
	result2.cvup<- which(cvup2==min(cvup2), arr.ind=TRUE)

	lam1.min.cvm <- unlist(lapply(c(1:nrow(result1.cvm)),function(x){lamda1[result1.cvm[x,1],result1.cvm[x,2]]}))
	lam1.index.cvm <- which(lam1.min.cvm==max(lam1.min.cvm), arr.ind=TRUE)
	mod1.cvm <- paracvfit1[[result1.cvm[lam1.index.cvm,1]]]

	lam1.min.cvup <- unlist(lapply(c(1:nrow(result1.cvup)),function(x){lamda1[result1.cvup[x,1],result1.cvup[x,2]]}))
	lam1.index.cvup <- which(lam1.min.cvup==max(lam1.min.cvup), arr.ind=TRUE)
	mod1.cvup <- paracvfit1[[result1.cvup[lam1.index.cvup,1]]]

	lam2.min.cvm <- unlist(lapply(c(1:nrow(result2.cvm)),function(x){lamda2[result2.cvm[x,1],result2.cvm[x,2]]}))
	lam2.index.cvm <- which(lam2.min.cvm==max(lam2.min.cvm), arr.ind=TRUE)
	mod2.cvm <- paracvfit2[[result2.cvm[lam2.index.cvm,1]]]

	lam2.min.cvup <- unlist(lapply(c(1:nrow(result2.cvup)),function(x){lamda2[result2.cvup[x,1],result2.cvup[x,2]]}))
	lam2.index.cvup <- which(lam2.min.cvup==max(lam2.min.cvup), arr.ind=TRUE)
	mod2.cvup <- paracvfit2[[result2.cvup[lam2.index.cvup,1]]]

	mod1.f.cvm.1se  <- do.call("cBind",coef(mod1.cvm, s="lambda.1se"))
	mod1.f.cvm.min  <- do.call("cBind",coef(mod1.cvm, s="lambda.min"))
	mod1.f.cvup.1se <- do.call("cBind",coef(mod1.cvup,s="lambda.1se"))
	mod1.f.cvup.min <- do.call("cBind",coef(mod1.cvup,s="lambda.min"))
	mod2.f.cvm.1se  <- do.call("cBind",coef(mod2.cvm, s="lambda.1se"))
	mod2.f.cvm.min  <- do.call("cBind",coef(mod2.cvm, s="lambda.min"))
	mod2.f.cvup.1se <- do.call("cBind",coef(mod2.cvup,s="lambda.1se"))
	mod2.f.cvup.min <- do.call("cBind",coef(mod2.cvup,s="lambda.min"))

	mod.coef <- cBind(mod1.f.cvm.1se,mod1.f.cvm.min,mod1.f.cvup.1se,mod1.f.cvup.min,
				mod2.f.cvm.1se,mod2.f.cvm.min,mod2.f.cvup.1se,mod2.f.cvup.min)
	colnames(mod.coef) <- rep(levels(group),8)

	mod.par <- matrix(c(result1.cvm[lam1.index.cvm,1]*0.01,mod1.cvm$lambda.1se,mod1.cvm$lambda.min,
			result1.cvup[lam1.index.cvup,1]*0.01,mod1.cvup$lambda.1se,mod1.cvup$lambda.min,
			result2.cvm[lam2.index.cvm,1]*0.01,mod2.cvm$lambda.1se,mod2.cvm$lambda.min,
			result2.cvup[lam2.index.cvup,1]*0.01,mod2.cvup$lambda.1se,mod2.cvup$lambda.min),
			4,3,byrow=TRUE)
	colnames(mod.par) <- c("alpha","lambda.1se","lambda.min")
	rownames(mod.par) <- c("cla.cvm","cla.cvup","mse.cvm","mse.cvup")

	result <- list()
	result[[1]] <- mod1.cvm
	result[[2]] <- mod1.cvup
	result[[3]] <- mod2.cvm
	result[[4]] <- mod2.cvup
	result[[5]] <- mod.coef
	result[[6]] <- mod.par
	return(result)
}


