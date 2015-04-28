################   Basic Multivariate Analysis   ###############
SA_uni <- function(physeq, samtype, spetype=NULL, col.sam=NULL, col.spe=NULL, shape=NULL, label=NULL, facet_formula=facet_formula, keepOnlyTheseTaxa=NULL, keepOnlyTheseSample=NULL, 
		distlist=c("bray","jaccard"), max.dist=0.3, threshold=0, axes=c(1,2), seq.depth=100, heatscale=heatscale, trans=trans, name, flip=FALSE, td.axes=NULL){
	#######    phyloseq detach   ########
	data <- data.frame(otuTable(physeq))
	if(speciesAreRows(physeq)==TRUE){	data <- t(data)	}
	tmp1 <- data[rowSums(data)!=0,]
	tmp2 <- data[rowSums(data)==0,]
	if(nrow(tmp2)){
		for(i in 1:nrow(tmp2)){
			n <- sample(2:10)
			if(ncol(tmp2)<10){n <- sample(2:ncol(tmp2))}
			tmp2[i,sample(1:ncol(tmp2))[1:n[1]]] <- 0.1/10^n[c(2:n[1],1)]
		}
		tmp2[is.na(tmp2)] <- 0
		data <- rbind(tmp1,tmp2)
	}
	otuTable(physeq) <- otuTable(data, speciesAreRows=FALSE)
	samData <- data.frame(samData(physeq))
	tree <- tre(physeq)
	taxonomy <- data.frame(taxTab(physeq))
	ID <- sample.names(physeq)
	feature <- species.names(physeq)
	group <- factor(samData[,samtype])
	type <- factor(taxonomy[,spetype])

	#######    Taxanomy Analysis   #######
	taxa.info <- taxa_analysis(physeq, samtype=samtype, color=col.sam, facet_formula=facet_formula, threshold=threshold, normalize=seq.depth, 
					keepOnlyTheseTaxa=keepOnlyTheseTaxa, keepOnlyTheseSample=keepOnlyTheseSample, OTUpoints=FALSE, labelOTUs=FALSE)
		#####   Abundane plot   #####
		png(file=paste(name, "tax.boxplot", "png", sep="."), width=6000, height=2000, res=300)
			print(taxa.info$tax.boxplot)
		dev.off()
		png(file=paste(name, "tax.barplot", "png", sep="."), width=6000, height=2000, res=300)
			print(taxa.info$tax.barplot)
		dev.off()
		png(file=paste(name, "tax.richness.plot", "png", sep="."), width=3200, height=2000, res=300)
			print(taxa.info$tax.richness.plot)
		dev.off()

	#####   Plot Group based raw heatmap   #####
	taxa.info[["raw.heatmap.ori"]] <- my.plot_heatmap(physeq, ps.ord=NULL, method=NULL, label=label, col.type=c(samtype, spetype), color=c(col.sam, col.spe), 
			heatscale=heatscale, keepOnlyTheseTaxa=keepOnlyTheseTaxa, keepOnlyTheseSample=keepOnlyTheseSample, normalize=seq.depth, trans=NULL)
	taxa.info[[paste("raw.heatmap.",trans$name,sep="")]] <- my.plot_heatmap(physeq, ps.ord=NULL, method=NULL, label=label, col.type=c(samtype, spetype), color=c(col.sam, col.spe), 
			heatscale=heatscale, keepOnlyTheseTaxa=keepOnlyTheseTaxa, keepOnlyTheseSample=keepOnlyTheseSample, normalize=seq.depth, trans=trans)
		#####   Define scale for heatmap  #####
		hei <- nspecies(physeq)
		wid <- nsamples(physeq)
		if(!is.null(keepOnlyTheseTaxa)){	hei <- sum(keepOnlyTheseTaxa%in%species.names(physeq))	}
		if(!is.null(keepOnlyTheseSample)){	hei <- sum(keepOnlyTheseSample%in%sample.names(physeq))	}

		png(file=paste(name, "tax.rawheatmap.ori", "png", sep="."), width=wid*50+500, height=hei*50+200, res=300)
			print(taxa.info$raw.heatmap.ori)
		dev.off()
		png(file=paste(name, "tax.rawheatmap", trans$name, "png", sep="."), width=wid*50+500, height=hei*50+200, res=300)
			print(taxa.info[[paste("raw.heatmap.",trans$name,sep="")]])
		dev.off()

	write.csv(taxa.info$tax.diversity, paste(name,".diversity.csv",sep=""))
	write.csv(taxa.info$tax.diver.test, paste(name,".diver.test.csv",sep=""))

	#######    Sample Distance Analysis   #######
		#######    Sample Distance Matrix    #######
		dist.all <- list()
		if(sum(getslots.phyloseq(physeq)%in%"tre")){
			registerDoParallel(cores=8)
			dist.all[["unifrac.w.nn"]][["distmatrix"]] <- distance(physeq, method="unifrac", weighted=TRUE, normalized=FALSE, parallel=TRUE, fast=TRUE)
			dist.all[["unifrac.w.nor"]][["distmatrix"]] <- distance(physeq, method="unifrac", weighted=TRUE, normalized=TRUE, parallel=TRUE, fast=TRUE)
			dist.all[["unifrac.nw"]][["distmatrix"]] <- distance(physeq, method="unifrac", weighted=FALSE, parallel=TRUE, fast=TRUE)
			dist.all[["dpcoa"]][["distmatrix"]] <- distance(physeq, "dpcoa")
			dist.all[["jsd"]][["distmatrix"]] <- distance(physeq, "jsd", parallel=TRUE)
		}
		for(i in 1:length(distlist)){	dist.all[[distlist[i]]][["distmatrix"]] <- distance(physeq, method=distlist[i])	}

		#######   Beta diversity Mental like Test   #######
		beta.test <- matrix(0, length(dist.all), 3)
		colnames(beta.test) <- c("Anosim","MRPP","perMANOVA")
		rownames(beta.test) <- names(dist.all)

	#######   Sample Distance Plot   #######
	for(i in 1:length(dist.all)){	
		distname <- names(dist.all)[i]
		#####   Beta Diversity Test   #####
		anosim <- anosim(dist.all[[distname]][["distmatrix"]], group, permutations=1000)
		mrpp <- mrpp(dist.all[[distname]][["distmatrix"]], group, permutations=1000)
		permanova <- adonis(dist.all[[distname]][["distmatrix"]]~group, permutations=1000)
		beta.test[i,1] <- anosim$signif
		beta.test[i,2] <- mrpp$Pvalue
		beta.test[i,3] <- permanova$aov.tab$"Pr(>F)"[1]

		write.csv(as.matrix(dist.all[[distname]][["distmatrix"]]), paste(name,"dist",distname,"csv",sep="."))
		Dis_ana <- distance_analysis(physeq, distmatrix=dist.all[[distname]][["distmatrix"]], axes=axes, samtype=samtype, color=col.sam, shape=shape, 
							keepOnlyTheseTaxa=keepOnlyTheseTaxa, keepOnlyTheseSample=keepOnlyTheseSample, trans=trans, max.dist=max.dist, 
							normalize=seq.depth, heatscale=heatscale, flip=flip, td.axes=td.axes, distname=distname, name=name)
		dist.all[[distname]][["dist.analysis"]] <- Dis_ana
		#####   Output Result   #####
			write.csv(Dis_ana$Signif$p.value, paste(name, "dist", distname, "pvalue", "csv", sep="."))
			png(file=paste(name, "dist", distname, "boxplot","png", sep="."), width=3200, height=2000, res=300)
				print(Dis_ana$boxplot)
			dev.off()
			png(file=paste(name, "dist", distname, "barplot","png", sep="."), width=3200, height=2000, res=300)
				print(Dis_ana$barplot)
			dev.off()
			png(file=paste(name, "dist", distname, "pcoaplot","png", sep="."), width=3200, height=3000, res=300)
				print(Dis_ana$pcoaplot)
			dev.off()
			png(file=paste(name, "dist", distname, "nmdsplot","png", sep="."), width=3200, height=3000, res=300)
				print(Dis_ana$nmdsplot)
			dev.off()
			if(length(names(Dis_ana$igraph))){
			for(j in 1:length(names(Dis_ana$igraph))){
				queue_name <- names(Dis_ana$igraph[j])
				png(file=paste(name, "dist", distname, queue_name, "png", sep="."), width=3200, height=3000, res=300)
					print(Dis_ana$igraph[[j]])
				dev.off()
			}}
	}
	write.csv(beta.test, paste(name, "dist", "beta.test", "csv", sep="."))

	#######   Ordination Analysis   #######
	ordinate.all <- list()
	#######   DPCoA   #######
	if(sum(getslots.phyloseq(physeq)%in%"tre")){	ordinate.all[["DPCoA"]][["ordination"]] <- ordinate(physeq, method="DPCoA")	}
	#######   Correspondence Analysis   #######
	ordinate.all[["UCA"]][["ordination"]] <- ordinate(physeq, method="CCA")
	ordinate.all[["CCA"]][["ordination"]] <- ordinate(as.formula(paste("physeq", samtype, sep="~")), method="CCA")
	ordinate.all[["DCA"]][["ordination"]] <- ordinate(physeq, method="DCA")
	#######   Redundancy Analysis   #######
	ordinate.all[["PCA"]][["ordination"]] <- ordinate(physeq, method="RDA")
	ordinate.all[["RDA"]][["ordination"]] <- ordinate(as.formula(paste("physeq", samtype, sep="~")), method="RDA")
	#######   NMDS Analysis   #######
	for(i in 1:length(distlist)){	ordinate.all[[paste("NMDS",distlist[i],sep=".")]][["ordination"]] <- ordinate(physeq, method="NMDS", distance=distlist[i])	}

	#######   Ordination Plot   #######
	for(i in 1:length(ordinate.all)){	
		ordname <- names(ordinate.all)[i]
		Ord_ana <- ordinate_analysis(physeq, ordination=ordinate.all[[ordname]][["ordination"]], col.type=c(samtype, spetype), color=c(col.sam, col.spe), 
							shape=shape, label=label, normalize=seq.depth, heatscale=heatscale, trans=trans, axes=axes, td.axes=td.axes, 
							keepOnlyTheseTaxa=keepOnlyTheseTaxa, keepOnlyTheseSample=keepOnlyTheseSample, ordname=ordname, name=name)
		ordinate.all[[ordname]][["ordi.analysis"]] <- Ord_ana
		#####   Output Result   #####
		png(file=paste(name, "ordi", ordname, "splitplot", "png", sep="."), width=6400, height=3000, res=300)
			print(Ord_ana$splitplot)
		dev.off()
		png(file=paste(name, "ordi", ordname, "biplot", "png", sep="."), width=3400, height=3000, res=300)
			print(Ord_ana$biplot)
		dev.off()
#		png(file=paste(name, "ordi", ordname, "supplot", "png", sep="."), width=3400, height=3000, res=300)
#			print(Ord_ana$supplot+par.axis+scale_color_manual(values=c(col.sam, col.spe)))
#		dev.off()
		#####   Define scale for heatmap  #####
		hei <- nspecies(physeq)
		wid <- nsamples(physeq)
		if(!is.null(keepOnlyTheseTaxa)){	hei <- sum(keepOnlyTheseTaxa%in%species.names(physeq))	}
		if(!is.null(keepOnlyTheseSample)){	hei <- sum(keepOnlyTheseSample%in%sample.names(physeq))	}
		png(file=paste(name, "ordi", ordname, "heatmap.ori", "png", sep="."), width=wid*50+500, height=hei*50+200, res=300)
			print(Ord_ana$heatmap.ori)
		dev.off()
		png(file=paste(name, "ordi", ordname, "heatmap", trans$name, "png", sep="."), width=wid*50+500, height=hei*50+200, res=300)
			print(Ord_ana[[paste("heatmap",trans$name,sep=".")]])
		dev.off()
	}

	#######   Return Result   #######
	result <- list()
	result[["beta.test"]] <- beta.test
	result[["taxonomy.info"]] <- taxa.info
	result[["distance.info"]] <- dist.all
	result[["ordination.info"]] <- ordinate.all
	return(result)

#	SA_basic(data, group, ID, name, uni=TRUE, unidist, sleep)
}

SA_basic <- function(data, group, ID, name, color, uni=FALSE, unidist=0, sleep=20){
	bcdist <- vegdist(data,"bray")
	write.csv(as.matrix(bcdist),paste(name,".dist.bc.csv",sep=""))

	##################   PCA, (D)PCoA, COA, NMDS   ###################
	result <- list()
	test <- matrix(0,5,2)
	colnames(test) <- c("Bray-Curtis","Unifrac")
	rownames(test) <- c("Anosim","MRPP","perManova","NMDS2","NMDS3")

	library(ade4)
	library(rgl)
		pca <- dudi.pca(data, scannf=FALSE, nf=3)
		coa <- dudi.coa(data, scannf=FALSE, nf=3)
		nmds2.bc <- metaMDS(bcdist, k=2)
		nmds3.bc <- metaMDS(bcdist, k=3)

		result[[1]] <- pca$eig[1:3]
		test[4,1] <- nmds2.bc$stress
		test[5,1] <- nmds3.bc$stress

		hclusterplot(bcdist, group, color, paste(name,".hcluster.bc.png", sep=""))
		nmdsplot(nmds2.bc, group, color, paste(name,".nmds2.bc.png",sep=""), cex=1, label=ID)
		nmdsplot(nmds3.bc, group, color, paste(name,".nmds3.bc.png",sep=""), cex=1, label=ID)
		Sys.sleep(sleep)

		if(uni==TRUE){
			nmds2.uni <- metaMDS(unidist, k=2)
			nmds3.uni <- metaMDS(unidist, k=3)
			test[4,2] <- nmds2.uni$stress
			test[5,2] <- nmds3.uni$stress
			hclusterplot(unidist, group, color, paste(name,".hcluster.uni.png",sep=""))
			nmdsplot(nmds2.uni, group, color, paste(name,".nmds2.uni.png",sep=""),cex=1, label=ID)
			nmdsplot(nmds3.uni, group, color, paste(name,".nmds3.uni.png",sep=""),cex=1, label=ID)
			Sys.sleep(sleep)
		}

		pcplot(pca, group, color, paste(name,".pca.png",sep=""), cbox=0.5, label=ID)
		Sys.sleep(sleep)
		pcplot(coa, group, color, paste(name,".coa.png",sep=""), cbox=0.5, label=ID)
		Sys.sleep(sleep)

	#########################################################
	con <- file(paste(name,".SA.log",sep=""), open="wt")
	sink(file=con, type="output")

	########   Multi-test    #########
		#######  ANOSIM, MRPP, perMANOVA  #######
		anosim.bc <- anosim(bcdist, group, permutations=1000)
		mrpp.bc <- mrpp(bcdist, group, permutations=1000)
		per.bc <- adonis(bcdist~group, permutations=1000)
#		perMANOVA.bc <- adonis(data ~ group, method="bray", permutations=1000)
		test[1,1] <- anosim.bc$signif
		test[2,1] <- mrpp.bc$Pvalue
		test[3,1] <- per.bc$aov.tab$"Pr(>F)"[1]

		print(anosim.bc)
		print(mrpp.bc)
		print(per.bc)

		if(uni==TRUE){
			anosim.uni <- anosim(unidist, group, permutations=1000)
			mrpp.uni <- mrpp(unidist, group, permutations=1000)
			per.uni <- adonis(unidist~group, permutations=1000)

			test[1,2] <- anosim.uni$signif
			test[2,2] <- mrpp.uni$Pvalue
			test[3,2] <- per.uni$aov.tab$"Pr(>F)"[1]

			print(anosim.uni)
			print(mrpp.uni)
			print(per.uni)
		}
	result[[2]] <- test

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

		eco <- data.frame(t(rbind(shannon, simpson, richness)))
		eco$group <- group
		write.csv(eco, paste(name,".diver.csv",sep=""))

		eco_ttest <- para_t.test(eco[,c(1,2,3,4,6)], eco$group, paired=FALSE, alternative="two.sided")
		eco_npttest <- para_npt.test(eco[,c(1,2,3,4,6)], eco$group, paired=FALSE, alternative="two.sided")
		write.csv(t(cbind(eco_ttest,eco_npttest)), paste(name,".diver.test.csv",sep=""))

	#########################################################
	##################         DA       #####################
	library(caret)
	library(lattice)
	library(rgl)
		###########    LDA    ##########
		lda <- lda(group~., data, CV=FALSE)
		lda.hat <- predict(lda, decision.values=TRUE)
		tabTrain <- confusionMatrix(group, lda.hat$class)
		denplot(lda.hat$x, group, color, paste(name,".lda.density.png",sep=""))
		if(ncol(lda.hat$x)>1){
			ldaplot(lda.hat$x, group, color, paste(name,".lda.xyplot.png",sep=""))
		}
		Sys.sleep(sleep)

		###########    PLSDA   #########
		useBayes <- plsda(data, group, ncomp=3, probMethod="Bayes")
		useSoftmax <- plsda(data, group, ncomp=3)
		useBayes.hat <- predict(useBayes, data)
		useSoftmax.hat <- predict(useSoftmax, data)
		print(confusionMatrix(predict(useBayes, data), group))
		print(confusionMatrix(predict(useSoftmax, data), group))

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

	print(c("feature"=ncol(data),"sample"=nrow(data)))
	print(table(group))
	print(result[[1]])
	print(result[[2]])

	sink()
	close(con)

	print(c("feature"=ncol(data),"sample"=nrow(data)))
	print(table(group))
	print(result[[1]])
	print(result[[2]])
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
	feature <- sort(colnames(data))
	########################################################
	################  basic test selection  ################
	library(coin)
		parattest <- para_t.test(data,group,paired=FALSE,alternative="two.sided")
		parachitest <- para_chi.test(data,group,paired=FALSE,alternative="two.sided")
		paranpttest <- para_npt.test(data,group,paired=FALSE,alternative="two.sided")
		result1 <- cbind(parattest,parachitest,paranpttest)
		result1 <- result1[feature,]
		write.csv(result1,paste(name,".fs.Test.pairwise.csv",sep=""))

		if(length(levels(group))>2){
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

	################   SPLSDA  Selection  ################
	library(spls)
	library(nnet)
		cv.lda <- cv.splsda(data.matrix(data), as.numeric(group), classifier="lda", K=c(1:5), eta=seq(0.1,0.9,0.1), scale.x=FALSE, fold=5, n.core=4)
		cv.log <- cv.splsda(data.matrix(data), as.numeric(group), classifier="logistic", K=c(1:5), eta=seq(0.1,0.9,0.1), scale.x=FALSE, fold=5, n.core=4)

		splsda.lda <- splsda(data.matrix(data), as.numeric(group), classifier="lda", eta=cv.lda$eta.opt, K=cv.lda$K.opt, scale.x=FALSE)
		splsda.log <- splsda(data.matrix(data), as.numeric(group), classifier="logistic", eta=cv.log$eta.opt, K=cv.log$K.opt, scale.x=FALSE)
		print(splsda.lda)
		print(splsda.log)

		par <- matrix(c(splsda.lda$eta,splsda.lda$K,splsda.lda$kappa,splsda.log$eta,splsda.log$K,splsda.lda$kappa),3,2)
		colnames(par) <- c(splsda.lda$classifier, splsda.log$classifier)
		rownames(par) <- c("eta","K","kappa")
		write.csv(splsda.lda$W, paste(name, ".fs.splsda.lda.csv", sep=""))
		write.csv(splsda.log$W, paste(name, ".fs.splsda.log.csv", sep=""))
		write.csv(par, paste(name, ".fs.splsda.par.csv", sep=""))
	detach("package:spls")
	detach("package:nnet")

	################  Indicator Species Analysis  #############
	library(labdsv)
		isa <- indval(data, group, numitr=1000)
#		summary(isa)
		summary(isa, p=0.05, digits=2, show=p)
		result2 <- data.frame(isa$indcls,isa$maxcls,isa$pval)
		result2 <- result2[feature,]
		isa.core <- result2[isa$pval<=0.05,]
		write.csv(isa.core, paste(name,".fs.isa.csv",sep=""))
	detach("package:labdsv")

	################  Random Forest Selection ############
	library(randomForest)
	library(Boruta)
		set.seed(125)
		boruta <- Boruta(group~., data=cbind(group,data), doTrace=2, maxRuns=12, ntree=500) #  Default maxRuns=4
		borutaplot(boruta,paste(name, ".fs.rf.png",sep=""))
		tmp <- TentativeRoughFix(boruta)
		result3 <- cbind(boruta$finalDecision,tmp$finalDecision)
		result3 <- result3[feature,]
		write.csv(data.frame(boruta$finalDecision,tmp$finalDecision), paste(name, ".fs.rf.csv",sep=""))
#		getConfirmedFormula(clean.Boruta)
#		getNonRejectedFormula(clean.Boruta)
	detach("package:Boruta")
	detach("package:randomForest")

	################   ELASTIC NET   ####################
	library(glmnet)
		elastic_mod <- paraelastic(data, group, family="multinomial")
		ENfsplot(elastic_mod,paste(name, ".fs.EN.png",sep=""))
		result4 <- as.matrix(elastic_mod[[5]])[feature,]
		write.csv(as.matrix(elastic_mod[[5]]),paste(name,".fs.EN.coef.csv",sep=""))
		write.csv(elastic_mod[[6]],paste(name,".fs.EN.par.csv",sep=""))
	detach("package:glmnet")	

	################      HHSVM      ####################

	################  SVM-RFE selection  ################
	library(e1071)
	library(kernlab)
		ftable <- para_svmrfe(data, group, ID, kernel="linear")
		result5 <- ftable
		write.csv(ftable,paste(name, ".fs.svmrfe.",kernel,".csv",sep=""))
	detach("package:e1071")
	detach("package:kernlab")

	################    Combine Result Score  #################
	score <- matrix(0,length(feature),6)
		rownames(score) <- feature
		colnames(score) <- c("ENET","ISA","RF","SVMRFE","Test","Score")
		score[,5] <- rowSums(cbind(as.numeric(result1[,3]),as.numeric(result1[,6]),as.numeric(result1[,7]))<=0.05)
		score[score>=2] <- 1	
		score[,2] <- result2[,3]<=0.05
		score[names(result5),4] <- result5
		result3[result3==3] <- 0
		score[,3] <- rowSums(result3/2)
		score[,1] <- rowSums(abs(result4)>0)/2
		score[,6] <- score[,1]/2+score[,2]/0.5+score[,3]+score[,4]/20+score[,5]/0.5
	write.csv(score, paste(name, ".fs.final.csv", sep=""))

	##################   Tune Score cutoff   ##################
		cre <- seq(0,max(score[,6]),by=0.1)
		flist <- lapply(cre, function(x){feature[score[,6]>=x]})

		accuracy <- sapply(flist,function(x){CL_basic(data[x], group, ID, "all", kernel="linear", plot=FALSE)})
		unlink("all.CL.log")
		rownames(accuracy) <- c("PLSDAlda", "PLSDAbayes", "SVM", "RandomForest")
		write.csv(t(rbind(cre, accuracy)), paste(name, ".fs.tune.csv", sep=""))

		png(file=paste(name,".fs.tune.png",sep=""), width=2000, height=1000, res=300)
			par(mar=c(2.5,2.5,2,1)+0.1, lwd=0.5, font.axis=1, cex.axis=1, pch=21, cex=1, mgp=c(2,0,0), tck=-0.01)
			plot(cre,accuracy[2,],main="",ylim=c(0,100),col="green",bg="green",ann=FALSE)
#			points(cre,accuracy[2,],col="cyan",bg="cyan")
			points(cre,accuracy[3,],col="red",bg="red")
			points(cre,accuracy[4,],col="blue",bg="blue")
			mtext(side=1, text="score", font=2, line=1)
			mtext(side=2, text="Accuracy", font=2, line=1)
			legend("bottomright",legend=c("PLSDA","SVM","RandomForest"),col=c("green","red","blue"),pt.bg=c("green","red","blue"),pch=21,cex=1)
		dev.off()
		##################################

##############################################

	################   Return part of the model   #############
#	return(elastic_mod)
	return(score)
}


CL_uni <- function(physeq, cat, name, kernel="linear", plot=TRUE){
	data <- data.frame(otuTable(physeq))
	samData <- data.frame(samData(physeq))
	tree <- tre(physeq)
	taxonomy <- data.frame(taxTab(physeq))

	ID <- sample.names(physeq)
	feature <- species.names(physeq)
	group <- samData[,cat]
	result <- CL_basic(data, group, ID, name, kernel=kernel, plot=plot)
	return(result)
}

CL_basic <- function(data, group, ID, name, kernel="linear", plot=TRUE){
	result <- rep(0,4)

	con <- file(paste(name,".CL.log",sep=""), open="wt")
	sink(file=con, type="output")

	library(caret)
	##########  PDA  ##########
	if(ncol(data)>2){
		###########    PLSDA   #########
		useBayes <- plsda(data, group, ncomp=3, probMethod="Bayes")
		useSoftmax <- plsda(data, group, ncomp=3)
		useBayes.hat <- predict(useBayes, data)
		useSoftmax.hat <- predict(useSoftmax, data)
		cof1 <- confusionMatrix(predict(useBayes, data), group)
		cof2 <- confusionMatrix(predict(useSoftmax, data), group)
		print(cof1)
		print(cof2)
		result[1] <- cof1$overall[1]*100
		result[2] <- cof2$overall[1]*100
		detach("package:klaR")
		detach("package:pls")
	}
	##########  SVM  ##########
	library(e1071)
		svm <- svm(group~., data, cross=nrow(data), kernel=kernel)
		svm.hat <- predict(svm,decision.values=TRUE)
		tabTrain <- confusionMatrix(group, svm.hat)
		cof3 <- summary(svm)
		print(cof3)
		result[3] <- cof3$tot.accuracy
	detach("package:e1071")
	##########  RF  ##########
	library(randomForest)
		count <- 20
		while(count > 0){
			trf <- randomForest(group~., data=data, proximity=TRUE)
			acc <- 100-trf$err.rate[nrow(trf$err.rate),1]*100
			if(result[4] < acc){
				result[4] <- acc;
				count <- 20;
				rf <- trf;
			}
			count <- count-1;
		}
		print(rf)
		if(plot==TRUE){
			rf.cmd <- cmdscale(1-rf$proximity)
			gn <- nlevels(group)
			color <- c("blue","red","green","purple")[1:gn]
			png(file=paste(name,".rf.mds.png",sep=""), width=1000, height=1000, res=300)
				par(mar=c(1.5,1.5,1,1)+0.1, font.axis=1, mgp=c(2,0.15,0), cex=0.8, tck=-0.005)
				plot(rf.cmd, pch=21, col=color[as.numeric(group)], bg=color[as.numeric(group)], ann=FALSE)
				setoff <- max(rf.cmd[,2])/30-min(rf.cmd[,2])/30
				s.label(cbind(rf.cmd[,1], rf.cmd[,2]+setoff), lab=ID, clabel=0.5, add.plot=TRUE, boxes=F, grid=F)
				legend("topright",legend=levels(group),col=color,pt.bg=color,pch=21,cex=1)
			dev.off()
		}
	detach("package:randomForest")
	################   Logestic Regression ##############
		if(ncol(data)<nrow(data)){
			glm.r <- glm(group~.,cbind(group,data), family=binomial(logit))
			print(summary(glm.r))
			print(predict(glm.r,type="response"))
		}
	detach("package:caret")

	sink()
	close(con)

	return(result)
}

################   Basic Linear Regression   ###############
LR_uni <- function(physeq, cat, name, family=gaussian()){
	data <- data.frame(otuTable(physeq))
	samData <- data.frame(samData(physeq))
	tree <- tre(physeq)
	taxonomy <- data.frame(taxTab(physeq))

	ID <- sample.names(physeq)
	feature <- species.names(physeq)
	RV <- as.vector(samData[,cat])
	LR_basic(data, RV, ID, name=name, family=family)
}

LR_basic <- function(data, RV, ID, name, family=gaussian()){
	fstr <- family$family
	linkstr <- family$link

	##############    ELASTIC NET   #################
	if(fstr=="gaussian"|fstr=="bionomial"|fstr=="poisson"){
		library(glmnet)
		elastic_mod <- paraelastic(data, RV, family=fstr)
		ENfsplot(elastic_mod,paste(name, ".LR.EN.png",sep=""))
		write.csv(as.matrix(elastic_mod[[5]]),paste(name,".LR.EN.coef.csv",sep=""))
		write.csv(elastic_mod[[6]],paste(name,".LR.EN.par.csv",sep=""))
		detach("package:glmnet")
	}

	#############   AIC for Combined Model   ############
	if(ncol(data)<nrow(data)){
			m0 <- lm(RV~1, data=cbind(RV,data))
			m1 <- lm(RV~., data=cbind(RV,data))
			print(summary(m1))
			print(extractAIC(m1))
			coefs <- FScoefs(m0, m1, data=cbind(RV,data))
			AICplot(coefs, paste(name,".AIC.coef.png",sep=""))
			m0 <- glm(RV~1, data=cbind(RV,data), family=family)
			m1 <- glm(RV~., data=cbind(RV,data), family=family)
			print(summary(m1))
			result <- step(m0, scope=list(lower=formula(m0),upper=formula(m1)), direction="forward", trace=FALSE)
			print(result)
			result <- step(m1, scope=list(lower=formula(m0),upper=formula(m1)), direction="backward", trace=FALSE)
			print(result)
	}

	##############     Linear Models   ###############
	if(ncol(data)>=nrow(data)){
		#############   LM for each features   ############
		lm_sep <- matrix(0,ncol(data),4)
		colnames(lm_sep) <- c("abu","pvalue","R2","R2.adj")
		rownames(lm_sep) <- colnames(data)
		dir.create(paste("./LR_",name,sep=""))
		setwd(paste("./LR_",name,sep=""))
		order <- order(colSums(data),decreasing=TRUE)
		for(j in 1:length(order)){
			i <- order[j]
			bcname <- colnames(data)[i]
			abu <- sum(data[,i])/sum(data)*100
			mod <- lm(RV~.,data=cbind(RV,data)[,c(1,i+1)])
			p <- 1-pf(summary(mod)$fstatistic[1],summary(mod)$fstatistic[2],summary(mod)$fstatistic[3])
			R2 <- summary(mod)$r.squared
			R2.adj <- summary(mod)$adj.r.squared
			lm_sep[bcname,1] <- abu
			lm_sep[bcname,2] <- p
			lm_sep[bcname,3] <- R2
			lm_sep[bcname,4] <- R2.adj
			if(p <= 0.05){
				LRplot(mod, paste(j,round(abu,2),bcname,"lm.png",sep="_"))
			}
		}
		setwd("../")
		write.csv(lm_sep,paste(name,".LR.lm.sep.csv",sep=""))
	}
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

################   AIC_FS_coefs for Regression   ###############
FScoefs <- function(m0, m1, data, trace=FALSE, direction="forward") {
	keepCoef <- function(m, aic){
		all <- names(coef(m1))
		new <- names(coef(m))
		ans <- rep(0, length(all))
		ans[match(new, all)] <- coef(m)
		ans
	}
	out <- with(data, stepAIC(m0, scope=list(lower=formula(m0), upper=formula(m1)), k=0, trace=trace, keep=keepCoef, direction=direction))
	rownames(out$keep) <- names(coef(m1))
	out$keep
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
	if(family=="multinomial"|family=="binomial"){
		paracvfit1 <- lapply(alpha, function(x){cv.glmnet(data.matrix(data),group,nfolds=5,alpha=x,family=family,type.measure="class")})
	}
	if(family=="gaussian"|family=="poisson"){
		paracvfit1 <- lapply(alpha, function(x){cv.glmnet(data.matrix(data),group,nfolds=5,alpha=x,family="gaussian",type.measure="deviance")})
	}
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

	mod.par <- matrix(c(result1.cvm[lam1.index.cvm,1]*0.01,mod1.cvm$lambda.1se,mod1.cvm$lambda.min,
			result1.cvup[lam1.index.cvup,1]*0.01,mod1.cvup$lambda.1se,mod1.cvup$lambda.min,
			result2.cvm[lam2.index.cvm,1]*0.01,mod2.cvm$lambda.1se,mod2.cvm$lambda.min,
			result2.cvup[lam2.index.cvup,1]*0.01,mod2.cvup$lambda.1se,mod2.cvup$lambda.min),
			4,3,byrow=TRUE)
	colnames(mod.par) <- c("alpha","lambda.1se","lambda.min")

	if(family=="multinomial"|family=="binomial"){
	rownames(mod.par) <- c("cla.cvm","cla.cvup","mse.cvm","mse.cvup")

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
	}

	if(family=="gaussian"|family=="poisson"){
	rownames(mod.par) <- c("dev.cvm","dev.cvup","mse.cvm","mse.cvup")

		mod.coef <- cBind(coef(mod1.cvm,s="lambda.1se"),coef(mod1.cvm,s="lambda.min"),coef(mod1.cvup,s="lambda.1se"),coef(mod1.cvup,s="lambda.min"),
				coef(mod2.cvm,s="lambda.1se"),coef(mod2.cvm,s="lambda.min"),coef(mod2.cvup,s="lambda.1se"),coef(mod2.cvup,s="lambda.min"))
		colnames(mod.coef) <- c("dev.cvm.1se","dev.cvm.min","dev.cvup.1se","dev.cvup.min","mse.cvm.1se","mse.cvm.min","mse.cvup.1se","mse.cvup.min")
	}

	result <- list()
	result[[1]] <- mod1.cvm
	result[[2]] <- mod1.cvup
	result[[3]] <- mod2.cvm
	result[[4]] <- mod2.cvup
	result[[5]] <- mod.coef
	result[[6]] <- mod.par
	return(result)
}


#########   Taxanomy Analysis Fundemental   #########
taxa_analysis <- function(physeq, samtype="Sample", spetype="TaxaGroup", color=NULL, keepOnlyTheseTaxa=NULL, keepOnlyTheseSample=NULL, 
				facet_formula=. ~ TaxaGroup, OTUpoints=FALSE, labelOTUs=FALSE, threshold=0, normalize=TRUE, trans=trans){
	Abu_df <- otu2df(physeq, keepOnlyTheseTaxa, keepOnlyTheseSample, threshold=threshold, normalize=TRUE)

	#######   boxplot   #######
	if(OTUpoints){
		pbox_hor <- ggplot(Abu_df, aes_string(x=spetype, y="Abu", fill=samtype)) + geom_boxplot(outlier.colour="NA") + 
				labs(y="Relative Abundance", x=spetype) + facet_grid(facet_formula, scales="free", space="free", drop=TRUE) + 
				geom_point(aes_string(color=samtype), position="jitter", size=1.5, alpha=I(1/2)) + 
				scale_fill_manual(values=color) + theme_bw() + opts(axis.text.x=theme_text(angle=90, hjust=1))
	} else {
		pbox_hor <- ggplot(Abu_df, aes_string(x=spetype, y="Abu", fill=samtype)) + geom_boxplot(outlier.size=1.5) + 
				labs(y="Relative Abundance", x=spetype) + facet_grid(facet_formula, scales="free", space="free", drop=TRUE) + 
				scale_fill_manual(values=color) + theme_bw() + opts(axis.text.x=theme_text(angle=90, hjust=1))
	}
	#####   Get facet formula   #####
	ft <- asOneFormula(facet_formula, as.formula(paste("~",samtype,"+",spetype,sep="")))
	var_form <- rev(strsplit(toString(ft),"\\, ")[[1]])[1]
	summary_formula <- as.formula(paste("Abu",var_form,sep="~"))
	Dsummary <- summaryBy(summary_formula, Abu_df, FUN=c(mean, sd, length))

	#####   Bar Plot with Error bar   #####
	pbar_hor <- ggplot(Dsummary, aes_string(x=spetype, y="Abu.mean", fill=samtype)) + geom_bar(stat="identity", position="dodge") + 
			geom_errorbar(aes(ymin=Abu.mean-Abu.sd/sqrt(Abu.length), ymax=Abu.mean+Abu.sd/sqrt(Abu.length)), width=.25, position=position_dodge(.9)) + 
			theme_bw() + facet_grid(facet_formula, scales="free", space="free", drop=TRUE) + opts(axis.text.x=theme_text(angle=90, hjust=1)) + 
			scale_fill_manual(values=color) + labs(y="Relative Abundance", x=spetype)
	if(OTUpoints){
		pbar_hor <- pbar_hor + geom_point(data=Abu_df, aes_string(x=spetype, y="Abu", color=samtype), position="jitter", size=1.5, alpha=I(1/2))
	}

	#####   Estimate Richness   #####
	unlockBinding("estimate_richness", as.environment("package:phyloseq"))
	assignInNamespace("estimate_richness", my.estimate_richness, ns="phyloseq", envir="package:phyloseq")
#	assign("estimate_richness", my.estimate_richness, "package:phyloseq")

	if(!is.null(keepOnlyTheseTaxa)){	physeq <- prune_species(keepOnlyTheseTaxa, physeq)	}
	if(!is.null(keepOnlyTheseSample)){	physeq <- prune_samples(keepOnlyTheseSample, physeq)	}

	eco <- my.estimate_richness(physeq, split=TRUE)
	group <- data.frame(samData(physeq))[,samtype]
	eco_ttest <- para_t.test(eco[,c(1,2,4,6,7)], group, paired=FALSE, alternative="two.sided")
	eco_npttest <- para_npt.test(eco[,c(1,2,4,6,7)], group, paired=FALSE, alternative="two.sided")
	eco_test <- t(cbind(eco_ttest,eco_npttest))
	rich.plot <- plot_richness_estimates(physeq, samtype, color=samtype) + scale_color_manual(values=color) + theme_bw()

	#####   Summary Result   #####
	result <- list()
	result[["tax.boxplot"]] <- pbox_hor
	result[["tax.barplot"]] <- pbar_hor
	result[["tax.diversity"]] <- eco
	result[["tax.diver.test"]] <- eco_test
	result[["tax.richness.plot"]] <- rich.plot
	return(result)
}

#########   Distance Analysis Fundemental   #########
distance_analysis <- function(physeq, distmatrix, axes=c(1,2), samtype=NULL, color=NULL, shape=NULL, normalize=TRUE, heatscale=c("white","steelblue"), max.dist=0.3, 
					keepOnlyTheseTaxa=NULL, keepOnlyTheseSample=NULL, trans=log_trans(4), title=NULL, flip=FALSE, td.axes=NULL, distname=NULL, name=NULL){
	dist <- list()
	group <- data.frame(samData(physeq))[,samtype]
	com <- rbind(levels(group),levels(group))
	com <- cbind(com, combn(levels(group),2))
	for(i in 1:ncol(com)){
		tmp <- as.vector(as.matrix(distmatrix)[group==com[1,i],group==com[2,i]])
		dist[[i]] <- tmp[tmp!=0]
	}
	distance <- unlist(dist)
	tmp <- rep(c(1:ncol(com)), sapply(dist,function(x){length(x)}))
	compare <- factor(tmp, labels=paste(com[1,],com[2,], sep=" vs. "))
	plotdata <- data.frame(compare, distance)
	Dsummary <- summaryBy(distance~compare, plotdata, FUN=c(mean, sd, length))
	signif <- pairwise.t.test(distance, compare, pool.sd=FALSE, p.adj="none")

	#######   boxplot   #######
	pbox <- ggplot(plotdata, aes(compare, distance)) + geom_boxplot() + theme_bw()
	if(flip==FALSE){
		pbox <- pbox + opts(axis.text.x=theme_text(size=8, angle=90), axis.text.y=theme_text(size=6),
				panel.grid.major=theme_blank(), panel.grid.minor=theme_line(colour='grey', linetype='dashed'), 
				legend.key = theme_blank(), axis.title.y=theme_blank(), axis.title.x=theme_blank())
	}
	if(flip==TRUE){
		pbox <- pbox + coord_flip() + opts(axis.text.x=theme_text(size=6), axis.text.y=theme_text(size=8),
				panel.grid.major=theme_blank(), panel.grid.minor=theme_line(colour='grey', linetype='dashed'), 
				legend.key = theme_blank(), axis.title.y=theme_blank(), axis.title.x=theme_blank())
	}

	#######   barplot   #######
	pbar <- ggplot(Dsummary, aes(compare, distance.mean)) + geom_bar(colour="black", fill="white") + theme_bw() + 
			geom_errorbar(aes(ymin=distance.mean-distance.sd/sqrt(distance.length), ymax=distance.mean+distance.sd/sqrt(distance.length)), width=.25, position=position_dodge(.9))

	if(flip==FALSE){
		pbar <- pbar + opts(axis.text.x=theme_text(size=8, angle=90), axis.text.y=theme_text(size=6),
				legend.key = theme_blank(), axis.title.y=theme_blank(), axis.title.x=theme_blank())
	}
	if(flip==TRUE){
		pbar <- pbar + coord_flip() + opts(axis.text.x=theme_text(size=6), axis.text.y=theme_text(size=8),
				legend.key = theme_blank(), axis.title.y=theme_blank(), axis.title.x=theme_blank())
	}

	#######   PCoA and NMDS   #######
	PCoA <- ordinate(physeq, method="PCoA", distance=distmatrix)
	NMDS <- ordinate(physeq, method="NMDS", distance=distmatrix)
		##   Notification: ordinate bug: line 239 and 243, failed to pass other par to metaMDS
		##   return(metaMDS(ps.dist)) should be: return(metaMDS(ps.dist, ...))
		##   If fixed, change code to: ordinate(physeq, method="NMDS", distance=distmatrix, k=3)
	if(!is.null(td.axes)){	NMDS.3d <- metaMDS(distmatrix, k=max(td.axes))	}

	#######   Plot Heatmap   #######
#	heat.pcoa <- my.plot_heatmap(physeq, PCoA, label="ID", species.label="Taxa")
#	heat.nmds <- my.plot_heatmap(physeq, method="NMDS", distmatrix, sample.label="ID", species.label="Taxa")

	#######   PCoA and NMDS scatter Plot   #######
	if(!is.null(keepOnlyTheseTaxa)){	physeq <- prune_species(keepOnlyTheseTaxa,   physeq)	}
	if(!is.null(keepOnlyTheseSample)){	physeq <- prune_samples(keepOnlyTheseSample, physeq)	}

	#######   2D scatter Plot   #######
	ppcoa <- plot_ordination(physeq, PCoA, type="samples", axes=axes, color=samtype, shape=shape) + geom_point(size=5) + 
			scale_color_manual(values=color) + 
#			scale_shape_manual(values=shape) + 
			theme_bw() + scale_x_continuous(name="PCo1") + scale_y_continuous(name="PCo2")
	pnmds <- plot_ordination(physeq, NMDS, type="samples", axes=axes, color=samtype, shape=shape) + geom_point(size=5) + 
			scale_color_manual(values=color) + 
#			scale_shape_manual(values=shape) + 
			theme_bw() + scale_x_continuous(name="NMDS1") + scale_y_continuous(name="NMDS2")

	#######   PCoA 3D Plot   #######
	if(!is.null(td.axes)){
		plotdata <- PCoA$vectors[,td.axes]
		if(!is.null(keepOnlyTheseSample)){
			plotdata <- plotdata[keepOnlyTheseSample,]
			group <- factor(as.vector(group[keepOnlyTheseSample]))
		}	
		cloud3d(plotdata, labels=group, filename=paste(name,"dist",distname,"pcoa3d","wrl",sep="."), type="vrml", pointstyle="s", 
			metalabels=rownames(plotdata), cols=color, scalefac=4, autoscale="independent", lab.axis=paste("PCo",td.axes,sep=""), 
			col.axis="darkblue", showaxis=TRUE, col.lab="darkgray", col.bg="white", cex.lab=1, htmlout=NULL, hwidth=1200, hheight=800, showlegend=TRUE)
		plotdata <- NMDS.3d$points[,td.axes]
		if(!is.null(keepOnlyTheseSample)){
			plotdata <- plotdata[keepOnlyTheseSample,]
			group <- factor(as.vector(group[keepOnlyTheseSample]))
		}
		cloud3d(plotdata, labels=group, filename=paste(name,"dist",distname,"nmds3d","wrl",sep="."), type="vrml", pointstyle="s", 
			metalabels=rownames(plotdata), cols=color, scalefac=4, autoscale="independent", lab.axis=paste("NMDS",td.axes,sep=""), 
			col.axis="darkblue", showaxis=TRUE, col.lab="darkgray", col.bg="white", cex.lab=1, htmlout=NULL, hwidth=1200, hheight=800, showlegend=TRUE)
	}

	#######   IGraph Network with max.dist.list  #######
	igraph <- list()
	for(i in 1:length(max.dist)){
		dist <- max.dist[i]
		ig <- make_network(physeq, type="samples", distance=distmatrix, max.dist=dist, keep.isolates=TRUE)
		if(length(ig[[3]])){
			igraph[[paste("igraph", round(dist,2), sep="_")]] <- plot_network(ig, physeq, color=samtype, shape=shape, line_weight=0.4, label=NULL)
		}
	}

	result <- list()
	result[["Dview"]] <- Dsummary
	result[["Signif"]] <- signif
	result[["PCoA"]] <- PCoA
	result[["NMDS"]] <- NMDS
	result[["boxplot"]] <- pbox
	result[["barplot"]] <- pbar
	result[["pcoaplot"]] <- ppcoa
	result[["nmdsplot"]] <- pnmds
#	result[["heatmap.pcoa"]] <- heat.pcoa
#	result[["heatmap.nmds"]] <- heat.nmds
	result[["igraph"]] <- igraph

	return(result)
}

#########   Ordination Analysis Fundemental   #########
ordinate_analysis <- function(physeq, ordination, axes=c(1,2), col.type=NULL, color=NULL, shape=NULL, label=NULL, normalize=TRUE, heatscale=c("white","steelblue"),
				keepOnlyTheseTaxa=NULL, keepOnlyTheseSample=NULL, trans=log_trans(4), title=NULL, td.axes=NULL, ordname, name){
	result <- list()

	#####   Plot Heatmap based on samtype and spetype   #####	
	result[["heatmap.ori"]] <- my.plot_heatmap(physeq, ps.ord=ordination, label=label, col.type=col.type, color=color, heatscale=heatscale, 
					keepOnlyTheseTaxa=keepOnlyTheseTaxa, keepOnlyTheseSample=keepOnlyTheseSample, normalize=normalize, trans=NULL)
	result[[paste("heatmap",trans$name,sep=".")]] <- my.plot_heatmap(physeq, ps.ord=ordination, label=label, col.type=col.type, color=color, heatscale=heatscale, 
					keepOnlyTheseTaxa=keepOnlyTheseTaxa, keepOnlyTheseSample=keepOnlyTheseSample, normalize=normalize, trans=trans)

	#####  re-Organize physeq for color/shape/label #####
	if(length(col.type)==2){
		slot1 <- col.type%in%sample.variables(physeq)
		if(sum(slot1)){	samData(physeq)$Color <- data.frame(samData(physeq))[,col.type[slot1][1]]	}
		color <- rev(color)
		slot2 <- col.type%in%rank.names(physeq)
		if(sum(slot2)){	taxTab(physeq) <- cbind(taxTab(physeq),Color=as.vector(taxTab(physeq)[,col.type[slot2][1]]))	}
		col.type <- "Color"
	}
	if(length(shape)==2){
		slot1 <- shape%in%sample.variables(physeq)
		if(sum(slot1)){	samData(physeq)$Shape <- data.frame(samData(physeq))[,shape[slot1][1]]	}
		shape <- rev(shape)
		slot2 <- shape%in%rank.names(physeq)
		if(sum(slot2)){	taxTab(physeq) <- cbind(taxTab(physeq),Shape=as.vector(taxTab(physeq)[,shape[slot2][1]]))	}
		shape <- "Shape"
	}
	if(length(label)==2){
		slot1 <- label%in%sample.variables(physeq)
		if(sum(slot1)){	samData(physeq)$Label <- data.frame(samData(physeq))[,label[slot1][1]]	}
		label <- rev(label)
		slot2 <- label%in%rank.names(physeq)
		if(sum(slot2)){	taxTab(physeq) <- cbind(taxTab(physeq),Label=as.vector(taxTab(physeq)[,label[slot2][1]]))	}
		label <- "Label"
	}

	#####   Subset Sample and Species  #####
	if(!is.null(keepOnlyTheseTaxa)){	physeq <- prune_species(keepOnlyTheseTaxa,   physeq)	}
	if(!is.null(keepOnlyTheseSample)){	physeq <- prune_samples(keepOnlyTheseSample, physeq)	}
	#####   Build Raw Ordination Plot DF   #####
	DF <- plot_ordination(physeq, ordination, axes, type="split", color=col.type, shape=shape, justDF=TRUE)
#	DF.3d <- plot_ordination(physeq, ordination, td.axes, type="split", color=col.type, shape=shape, justDF=TRUE)
	#####   Split Plot   #####
	result[["splitplot"]] <- ggplot(DF, aes_string(x=names(DF)[1], y=names(DF)[2], color=col.type, shape=shape, na.rm=TRUE)) + theme_bw() + 
					geom_point(aes(size=id.type), na.rm=TRUE) + scale_color_manual(values=color) + scale_size_manual("type", values=c(samples=5, species=2)) + 
					geom_text(size=4, aes_string(label=label), vjust=2, na.rm=TRUE) + facet_wrap(~id.type, nrow=1)
	#####  Biplot with Arrow  #####
	result[["biplot"]] <- ggplot(DF, aes_string(x=names(DF)[1], y=names(DF)[2], color=col.type, shape=shape, na.rm=TRUE)) + theme_bw() + 
					geom_point(subset = .(id.type=="samples"), aes(size=id.type), na.rm=TRUE) + 
					geom_segment(subset = .(id.type=="species"), aes_string(x=0, y=0, xend=names(DF)[1], yend=names(DF)[2], colour=col.type), size=1, arrow=arrow(length=unit(0.2,"cm")), alpha=0.75) + 
					scale_color_manual(values=color) + scale_size_manual("type", values=c(samples=5, species=1)) + 
					geom_text(size=4, aes_string(label=label), vjust=2, na.rm=TRUE) + geom_hline(aes(0), size=.2) + geom_vline(aes(0), size=.2)
	#####   Stress/Eigen Plot   #####
	if(sum(class(ordination)%in%c("metaMDS","monoMDS"))){
		#####  NMDS Stress Info  #####
		k <- ncol(ordination$points)
		gof <- goodness(ordination)
		if(sum(gof)){cexp <- gof/mean(gof)*0.8} else{cexp <- 0.8}
		result[["stress"]] <- ordination$stress
		#####  Stress Plot  #####
#		stressplot(ordination, p.col="blue", l.col="red", lwd=2, cex=0.8, ann=FALSE)
#		mtext(side=1, text="Observed Dissimilarity", font=2, cex=0.9, line=1)
#		mtext(side=2, text="Ordination Distance", font=2, cex=0.9, line=1)
#		result[["stressplot"]] <- stressplot
		#####  Biplot  #####
		result[["biplot"]] <- ggplot(DF, aes_string(x=names(DF)[1], y=names(DF)[2], color=col.type, shape=shape, na.rm=TRUE)) + theme_bw() + 
					geom_point(aes(size=id.type), na.rm=TRUE) + scale_size_manual("type", values=c(samples=5, species=2)) + 
					geom_text(size=4, aes_string(label=label), vjust=2, na.rm=TRUE) + scale_color_manual(values=color)
	} else if(sum(class(ordination)%in%"cca")){
		#####  Eigen Value  #####
		if(sum(ordination$CCA$eig)){eig <- ordination$CCA$eig} else{eig <- ordination$CA$eig}
		result[["eigen"]] <- eig
#		result[["eigenplot"]] <- plot_eigen(vector=100*eig/sum(eig))
	} else if(class(ordination)=="dpcoa"){
		#####  Eigen Value  #####
		eig <- ordination$eig
		result[["eigen"]] <- eig
#		result[["eigenplot"]] <- plot_eigen(vector=100*eig/sum(eig))
	} else{
		eigen <- ordination$evals
		decorana <- ordination$evals.decorana
		origin <- ordination$origin
		result[["eigen"]] <- eigen
#		result[["eigenplot"]] <- plot_eigen
	}
	return(result)
}

#######   Modified Version of phyloseq:::otu2df   #######
otu2df <- function(physeq, keepOnlyTheseTaxa=NULL, keepOnlyTheseSample=NULL, threshold=0, normalize=TRUE){
#taxMap=rank.names(physeq), samMap=sample.variables(physeq), 
	#######  Define Order  #######
	speorder <- species.names(physeq)
	samorder <- sample.names(physeq)
	if(speciesAreRows(physeq)){	physeq <- t(physeq)	}
	dfdata <- data.frame(otuTable(physeq))[samorder, speorder]
	# Normalization Part
	if(is.numeric(normalize)){
		dfdata <- dfdata/normalize
	} else if(normalize){
		dfdata <- dfdata/rowSums(dfdata)
	} else {	dfdata <- dfdata	}
	samples <- data.frame(samData(physeq))[samorder, ]
	species <- data.frame(taxTab(physeq))[speorder, ]
	#######  Organize Tables  #######
	if(!is.null(keepOnlyTheseTaxa)){
		dfdata <- dfdata[,speorder %in% keepOnlyTheseTaxa]
		species <- subset(species, speorder %in% keepOnlyTheseTaxa)
		species <- droplevels(species)
	}
	if(!is.null(keepOnlyTheseSample)){
		dfdata <- dfdata[samorder %in% keepOnlyTheseSample,]
		samples <- subset(samples, samorder %in% keepOnlyTheseSample)
		samples <- droplevels(samples)
	}
	#######  Reshape Abundance Table  #######
	tmp.Abu <- reshape(dfdata, direction="long", varying=colnames(dfdata), times=colnames(dfdata), v.names="Abu", timevar="TaxaGroup", idvar="Sample", ids=rownames(dfdata))
	tmp.sam <- samples[rep(1:nrow(samples),times=nrow(species)),]
	tmp.spe <- species[rep(1:nrow(species), each=nrow(samples)),]
	df <- data.frame(tmp.Abu, tmp.sam, tmp.spe)
	#######  Apply Restriction to df  ########
	df <- df[df$Abu>=threshold,]
	return(df)
}

my.estimate_richness <- function(physeq, split=TRUE){
	# Check for singletons, and then warning if they are missing.
	# These metrics only really meaningful if singletons are included.
	if(!any(otuTable(physeq)==1)){
		warning("The experiment object you have provided does not have\n",
		"any singletons. This is highly suspicious. Results of richness\n",
		"estimates are probably unreliable, or wrong, if you have already\n",
		"trimmed low-abundance taxa from the data.\n",
		"It is recommended that you find the un-trimmed data and retry.")
	}
	# If we are not splitting sample-wise, sum the species. Else, enforce orientation.
	if(!split){
		OTU <- speciesSums(physeq)
	} else if(split){
		OTU <- as(otuTable(physeq), "matrix")
		if( speciesAreRows(physeq) ){ OTU <- t(OTU) }
	}
	# Some standard richness parameters
	shannon <- diversity(OTU, index="shannon", MARGIN=1, base=exp(1))
	simpson <- diversity(OTU, index="simpson", MARGIN=1, base=exp(1))
	if(sum(ceiling(OTU)) > sum(OTU)){
		warning("Ceiling Data were used\n",
		"Results of richness estimates are probably unreliable.\n")
	}
	richness <- estimateR(ceiling(otuTable(physeq)))
#		Srar <- rarefy(Data_C, ceiling(min(rowSums(Data_C))))
#		sac <- specaccum(Data)
#		plot(sac, ci.type="polygon", ci.col="yellow")
	return(data.frame(t(rbind(richness, shannon, simpson))))
}

my.plot_network <- function(g, physeq=NULL, type="samples", color=NULL, shape=NULL, point_size=4, alpha=1, 
        label="value", hjust=1.35, line_weight=0.5, line_color=color, line_alpha=0.4, layout.method=layout.fruchterman.reingold){
	# Make the edge-coordinates data.frame
	edgeDF    <- data.frame(get.edgelist(g))
	edgeDF$id <- 1:length(edgeDF[, 1])
	# Make the vertices-coordinates data.frame
	vertDF    <- layout.method(g)
	colnames(vertDF) <- c("x", "y")
	vertDF    <- data.frame(value=g[[9]][[3]][["name"]], vertDF)

	# If phyloseq object provided,
	# AND it has the relevant additional data
	# THEN add it to vertDF
	if( !is.null(physeq) ){
		extraData <- NULL
		if( type == "samples" & !is.null(sampleData(physeq, FALSE)) ){
			extraData <- sampleData(physeq)[as.character(vertDF$value), ]
		} else if( type == "species" & !is.null(taxTab(physeq, FALSE)) ){
			extraData <- taxTab(physeq)[as.character(vertDF$value), ]
		}
		# Only mod vertDF if extraData exists
		if( !is.null(extraData) ){
			vertDF <- data.frame(vertDF, extraData)
		}
	}

	# Combine vertex and edge coordinate data.frames
	graphDF   <- merge(melt(edgeDF, id="id"), vertDF, by = "value")
	# Initialize the ggplot
	p <- ggplot(vertDF, aes(x, y))
	# Strip all the typical annotations from the plot, leave the legend
	p <- p + theme_bw() + opts(panel.grid.major=theme_blank(), panel.grid.minor=theme_blank(), axis.ticks=theme_blank(), panel.border=theme_blank(), 
					axis.text.x=theme_blank(), axis.text.y=theme_blank(), axis.title.x=theme_blank(), axis.title.y=theme_blank())
	# Add the graph vertices as points
	p <- p + geom_point(aes_string(color=color, shape=shape), size=point_size) + scale_color_manual(values=colval) + scale_shape_manual(values=shaval)
	# Add the text labels
	if( !is.null(label) ){	p <- p + geom_text(aes_string(label=label), size = 2, hjust=hjust)}
	# Add the edges:
	p <- p + geom_line(aes_string(group="id", color=line_color), graphDF, size=line_weight, alpha=line_alpha)
	return(p)
}

my.plot_heatmap <- function(physeq, ps.ord=NULL, method="NMDS", distance="bray", label=NULL, col.type=NULL, color=NULL, normalize=TRUE, heatscale=c("#000033", "#66CCFF"), 
				keepOnlyTheseTaxa=NULL, keepOnlyTheseSample=NULL, na.value=NULL, trans=log_trans(4), breaks=c(1,2,4,8,16,32,60,80,100), max.label=250, ...){
##   Modified version of plot_heatmap function in phyloseq
##   1, Can take ordination-class as input if no method and distance are provided
##   2, Flexible Color/Label defination for axis-label by samtype and/or spetype with Color Legend
##   3, Customarize keepOnlyTheseTaxa, keepOnlyTheseSample, Normalization.

	## Define Label and Color by Group
	sample.label <- species.label <- NULL
	sample.label <- label[label%in%sample.variables(physeq)][1]
	label <- rev(label)
	species.label <- label[label%in%rank.names(physeq)][1]

	sample.color <- species.color <- NULL
	samtype <- col.type[col.type%in%sample.variables(physeq)][1]
	col.type <- rev(col.type)
	spetype <- col.type[col.type%in%rank.names(physeq)][1]

	# Copy the approach from NeatMap and Phyloseq by doing ordination on samples, can take ordination subject for analysis
	# Get ps.ord ordination-class if ps.ord is not provided	
	if(is.null(ps.ord) & !is.null(method)){
		# Capture the NMDS iterations cat() output with capture.output
		junk <- capture.output( ps.ord <- ordinate(physeq, method, distance, ...), file=NULL)
	}

	# Enforce orientation and transfer otu table into matrix class
	if(!speciesAreRows(physeq)){	physeq <- t(physeq)	}
	mot <- as(otuTable(physeq), "matrix")
	# Keep selected info and Define normalization or not. 
	if(is.numeric(normalize)){
		mot <- mot/normalize*100
	} else if(normalize){
		mot <- sapply(1:ncol(mot),function(x){mot[,x]/sum(mot[,x])}*100)
	} else {	mot <- mot	}
	colnames(mot) <- sample.names(physeq)

	# Group label by samtype and spetype and then trim Taxa and Sample
	if(!is.null(keepOnlyTheseTaxa)){	mot <- mot[keepOnlyTheseTaxa,]	}
	if(!is.null(keepOnlyTheseSample)){	mot <- mot[,keepOnlyTheseSample]	}
	if(!is.na(samtype)){
		tmp <- data.frame(samData(physeq)[colnames(mot), samtype])
		mot <- mot[, rownames(tmp)[order(tmp)]]
	}
	if(!is.na(spetype)){
		tmp <- data.frame(taxTab(physeq)[rownames(mot), spetype])
		mot <- mot[rownames(tmp)[order(tmp)], ]
	}

	# Initialize sample and species order vectors as NULL
	# Reorder by the angle in radial coordinates on the 2-axis plane.
	species.order <- rownames(mot)
	sample.order <- colnames(mot)
	if(!is.null(ps.ord)){
		reduction.result <- scores(ps.ord, choices=c(1, 2), display="sites")
		sample.order <- sample.names(physeq)[order(RadialCoords(reduction.result)[, "theta"])]
		sample.order <- sample.order[sample.order%in%colnames(mot)]
		mot <- mot[, sample.order]
		# Re-order species if possible
		test <- try(scores(ps.ord, choices=c(1, 2), display="species"), TRUE)
		if( class(test) != "try-error" & !is.null(test) ){
			species.reduct <- scores(ps.ord, choices=c(1, 2), display="species")
			species.order  <- species.names(physeq)[order(RadialCoords(species.reduct)[, "theta"])]
			species.order <- species.order[species.order%in%rownames(mot)]
			mot <- mot[species.order, ]
		}
	}

	# Alternative to melt that won't re-order your carefully-re-ordered stuff...
	adf <- data.frame(expand.grid(y=rownames(mot), x=colnames(mot)), value=as(mot, "vector"))

	# Now the plotting part
	p <- ggplot(adf, aes(x=as.factor(x), y=as.factor(y), fill=value)) + geom_tile() + opts(axis.text.x=theme_blank(),axis.text.y=theme_blank())
	p <- update_labels(p, list(fill="Abundance", y="", x=""))

	## Define Color by Group
	if(!is.null(sample.label)){
		sample.color <- rep("black", length(sample.order))
		if(!is.null(samtype)){if(!is.na(samtype)){
			sam.group <- data.frame(samData(physeq))[sample.order, samtype]
			sample.color <- color[as.vector(sam.group)]
		}}
	}
	if(!is.null(species.label)){
		species.color <- rep("black", length(species.order))
		if(!is.null(spetype)){if(!is.na(spetype)){
			spe.group <- data.frame(taxTab(physeq))[species.order, spetype]
			species.color <- color[as.vector(spe.group)]
		}}
	}

	## Axis Relabeling (Skipped if more than max.label):
	# Re-write sample-labels to some sample variable...
	if(!is.null(sample.label) & length(sample.order)<=max.label){if(!is.na(sample.label)){
		p <- p + opts(axis.text.x=theme_text(size=2*treetextsize(0.1*length(sample.order)), colour=sample.color, angle = -90, hjust = 0))
		# Make a sample-named vector of the values for sample.label
		labvec <- as(getVariable(physeq, sample.label), "vector")
		names(labvec) <- sample.names(physeq)
		if(!is.null(sample.order)){
			# Re-order according to sample.order
			labvec <- labvec[sample.order]
		}
		# Add the sample.label re-labeling layer
		p <- p + scale_x_discrete(sample.label, labels=labvec)
	}}

	if(!is.null(species.label) & length(species.order)<=max.label){if(!is.na(species.label)){
		# Make a species-named vector of the values for species.label
		p <- p + opts(axis.text.y=theme_text(size=2*treetextsize(0.1*length(species.order)), colour=species.color, hjust=1))
		labvec <- as(taxTab(physeq)[, species.label], "vector")
		names(labvec) <- species.names(physeq)
		if( !is.null(species.order) ){
			# Re-order according to species.order
			labvec <- labvec[species.order]
		}
		# Add the species.label re-labeling layer
		p <- p + scale_y_discrete(species.label, labels=labvec)
	}}

	## Define na.value color
	if(is.null(na.value)){
		col.va <- round(col2rgb(heatscale[1])*1.2-col2rgb(heatscale[2])*0.2)
		col.va[col.va<0] <- 0
		col.va[col.va>255] <- 255
		na.value <- rgb(t(col.va), alpha=255, max=255)
	}

	# Color scale transformations
	if(!is.null(trans)){
		p <- p + scale_fill_gradient(low=heatscale[1], high=heatscale[2], trans=trans, breaks=breaks, na.value=na.value)
	} else {
		p <- p + scale_fill_gradient(low=heatscale[1], high=heatscale[2], na.value=na.value)
	}
	return(p)
}

####   Attach RadialCoords function  ####
RadialCoords <- function(pos){
	pos <- as.matrix(pos)
	n.points <- nrow(pos)
	radial <- matrix(0, n.points, 2)
	xc <- mean(pos[,1])
	yc <- mean(pos[,2])
	for(i in 1:n.points){
		radial[i,1] <- sqrt((pos[i,2]-yc)^2+(pos[i,1]-xc)^2)
		radial[i,2] <- atan2(pos[i,2]-yc, pos[i,1]-xc)
	}
	rownames(radial) <- rownames(pos)
	colnames(radial) <- c("r", "theta")
	return(radial)
}

####   Attach treetextsize function  ####
treetextsize <- function(n){
	# empirically chosen size-value calculator.
	s <- 6 * exp(-n/100)
	# enforce a floor.
	s <- ifelse(s > 0.5, s, 0.5)
	# enforce a max
	s <- ifelse(s < 4, s, 4)
	return(s)
}

####   Attach Plot3d from corrplot.r   ####
rgl_plot3d <- function(plotdata, group=NULL, color, lab=NULL, type="s", size=1){
library(rgl)
	if(!is.null(group)){	color <- color[as.vector(group)]	}
	if(length(type)>1){	type <- type[as.vector(group)]	}
	if(length(size)>1){	size <- size[as.vector(group)]	}
	plot3d(plotdata, col=color, type=type, size=size, xlab=lab[1], ylab=lab[2], zlab=lab[3])
detach("package:rgl")
}




