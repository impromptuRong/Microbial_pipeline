########  plot for PCA, COA
pcplot <- function(pca, cat, filename, cbox=1, label=FALSE){
	kip <- round(100*pca$eig/sum(pca$eig),digits=1)
#	cumsum(kip)

	png(file=filename, width=2000, height=2000, res=250)
#	scatter(pca, grid=F)
	catlog <- factor(cat)
	indic <- levels(catlog)
	color <- c("blue", "red", "green", "purple")
	par(mfrow=c(2,2))
	s.class(dfxy=pca$li, fac=catlog, grid=F, col=color, clabel=cbox, xax=1, yax=2, possub="bottomright", csub = 1.2, 
		sub=paste("x: PC1","(",kip[1],"% of inertia)","\ny: PC2","(", kip[2],"% of inertia)"), addaxes=TRUE)
	if(length(label)>1){	s.label(cbind(pca$li[,1], pca$li[,2]+0.1), lab=label, clabel=0.5, add.plot=TRUE, boxes=F, grid=F)	}

	s.class(dfxy=pca$li, fac=catlog, grid=F, col=color, clabel=cbox, xax=1, yax=3, possub="bottomright", csub = 1.2, 
		sub=paste("x: PC1","(",kip[1],"% of inertia)","\ny: PC3","(", kip[3],"% of inertia)"), addaxes=TRUE)
	if(length(label)>1){	s.label(cbind(pca$li[,1], pca$li[,3]+0.1), lab=label, clabel=0.5, add.plot=TRUE, boxes=F, grid=F)	}

	s.class(dfxy=pca$li, fac=catlog, grid=F, col=color, clabel=cbox, xax=2, yax=3, possub="bottomright", csub = 1.2, 
		sub=paste("x: PC2","(",kip[2],"% of inertia)","\ny: PC3","(", kip[3],"% of inertia)"), addaxes=TRUE)
	if(length(label)>1){	s.label(cbind(pca$li[,2], pca$li[,3]+0.1), lab=label, clabel=0.5, add.plot=TRUE, boxes=F, grid=F)	}

	label <- paste("x: PC1(",kip[1],"%)\ny: PC2(",kip[2],"%)\ny: PC3(",kip[3],"%)\n")
	barplot(kip, main="Eigen Values",sub=label,ylim=c(0,100))
	legend("bottomright",legend=indic,col=color,pt.bg=color,pch=21)
	lines(c(1:length(kip)),cumsum(kip), type='l', lwd=2, col="red",ann=FALSE,ylim = c(0,100))
	#s.corcircle(pca$co,sub="Correlation Circle", grid=F, csub=1.2)##Plot of the factorial maps of a correlation circle

	dev.off()

	tp2col <- as.vector(cat)
	for(i in 1:length(indic)){
		tp2col[tp2col==indic[i]] <- color[i]
	}
	plot3d(pca$li,col=tp2col,xlab="PC1",ylab="PC2",zlab="PC3",type="s",size=1)
}

########  plot for NMDS
nmdsplot <- function(nmds, cat, filename, cex=0.8, label=FALSE){
k <- ncol(nmds$points)
ud <- floor(min(nmds$points))
up <- ceiling(max(nmds$points))
gof <- goodness(nmds)
if(sum(gof)){	cexp <- (gof/mean(gof))*cex	}
else{	cexp <- cex	}
color <- c("blue", "red", "green", "purple")
color <- factor(cat,label=color[1:length(levels(cat))])

if(k==2){
	png(file=filename, width=2000, height=1000, res=250)
	par(mfrow=c(1,2), mar=c(2.5,3,0.5,1)+0.1, lwd=0.5, font.axis=2,cex.axis=cex, mgp=c(2,0.25,0), tck=-0.01)
	
#	fig <- ordiplot(nmds, xlab=NULL, ylab=NULL, ann=FALSE, choices=c(1,2), type="none", cex=cex, display="sites")
#	type="t" for text
#	points(fig$sites,pch=21,cex=cexp,col=as.vector(color),bg=as.vector(color))
	plot(nmds$points,ann=FALSE,pch=21,cex=cexp,col=as.vector(color),bg=as.vector(color))
	mtext(side=1, text="NMDS1", font=2, cex=cex+0.1, line=1)
	mtext(side=2, text="NMDS2", font=2, cex=cex+0.1, line=1.5)
#xlim=c(ud,up),ylim=c(ud,up),

#	if(label){	s.label(cbind(fig$sites[,1],fig$sites[,2]+0.1), lab=label, clabel=0.5, add.plot=TRUE, boxes=F, grid=F)	}
	if(length(label)>1){	s.label(cbind(nmds$points[,1],nmds$points[,2]+0.1), lab=label, clabel=0.3, add.plot=TRUE, boxes=F, grid=F)	}

	legend("bottomright",legend=levels(cat),col=levels(color),pt.bg=levels(color),pch=21,cex=cex)
	legend("topleft",legend=paste("Stress=",round(nmds$stress,digits=3),sep=""),bty="n",cex=cex)

	stressplot(nmds, p.col="blue", l.col="red", lwd=2, cex=cex, ann=FALSE)
	mtext(side=1, text="Observed Dissimilarity", font=2, cex=cex+0.1, line=1)
	mtext(side=2, text="Ordination Distance", font=2, cex=cex+0.1, line=1.5)

	dev.off()
}

if(k==3){
	png(file=filename, width=2000, height=2000, res=250)
	par(mfrow=c(2,2), mar=c(2.5,3,0.5,1)+0.1, lwd=0.5, font.axis=2,cex.axis=cex, mgp=c(2,0.25,0), tck=-0.01)

	plot(nmds$points[,-3],ann=FALSE,pch=21,cex=cexp,col=as.vector(color),bg=as.vector(color))
	mtext(side=1, text="NMDS1", font=2, cex=cex+0.1, line=1.5)
	mtext(side=2, text="NMDS2", font=2, cex=cex+0.1, line=1.5)
	if(length(label)>1){	s.label(cbind(nmds$points[,1],nmds$points[,2]+0.1), lab=label, clabel=0.5, add.plot=TRUE, boxes=F, grid=F)	}
	legend("bottomright",legend=levels(cat),col=levels(color),pt.bg=levels(color),pch=21,cex=cex)
#xlim=c(ud,up),ylim=c(ud,up),

	plot(nmds$points[,-2],xlim=c(ud,up),ylim=c(ud,up),ann=FALSE,pch=21,cex=cexp,col=as.vector(color),bg=as.vector(color))
	mtext(side=1, text="NMDS1", font=2, cex=cex+0.1, line=1.5)
	mtext(side=2, text="NMDS3", font=2, cex=cex+0.1, line=1.5)
	if(length(label)>1){	s.label(cbind(nmds$points[,1],nmds$points[,3]+0.1), lab=label, clabel=0.3, add.plot=TRUE, boxes=F, grid=F)	}
	legend("bottomleft",legend=paste("Stress=",round(nmds$stress,digits=3),sep=""),bty="n",cex=cex)

	plot(nmds$points[,-1],xlim=c(ud,up),ylim=c(ud,up),ann=FALSE,pch=21,cex=cexp,col=as.vector(color),bg=as.vector(color))
	mtext(side=1, text="NMDS2", font=2, cex=cex+0.1, line=1.5)
	mtext(side=2, text="NMDS3", font=2, cex=cex+0.1, line=1.5)
	if(length(label)>1){	s.label(cbind(nmds$points[,2],nmds$points[,3]+0.1), lab=label, clabel=0.3, add.plot=TRUE, boxes=F, grid=F)	}
#	legend("topright",legend=levels(cat),col=levels(color),pt.bg=levels(color),pch=21,cex=cex)
#	legend("topleft",legend=paste("Stress=",round(nmds$stress,digits=3),sep=""),bty="n",cex=cex)

	stressplot(nmds, p.col="blue", l.col="red", lwd=2, cex=cex, ann=FALSE)
	mtext(side=1, text="Observed Dissimilarity", font=2, cex=cex+0.1, line=1.5)
	mtext(side=2, text="Ordination Distance", font=2, cex=cex+0.1, line=1.5)

	dev.off()

#	tp2col <- as.vector(cat)
#	for(i in 1:length(indic)){
#		tp2col[tp2col==indic[i]] <- color[i]
#	}
	plot3d(nmds$points,col=as.vector(color),xlab="NMDS1",ylab="NMDS2",zlab="NMDS3",type="s",size=1)
}

}

denplot <- function(model, group, filename){
	plotdata <- model$x
	npic <- ncol(plotdata)
	png(file=filename, width=1000*npic, height=1000, res=250)
	indic <- levels(group)
	color <- c("blue", "red", "green", "purple")
	for(i in 1:npic){
		plot <- densityplot(~plotdata[,i], groups=group, col=color, xlab=paste("LD",i,sep=""), ylab="Density")
		print(plot, position=c((i-1)/npic,0,i/npic,1), split=c(1,1,1,1), more=TRUE)
	}
#	legend("topright",legend=indic,col=color,pch=21,cex=0.5)
	dev.off()
}

ldaplot <- function(model, group, filename){
	plotdata <- model$x
	cb <- combn(ncol(plotdata),2)
	npic <- ncol(cb)

	png(file=filename, width=1000*npic, height=1000, res=250)
	indic <- levels(group)
	color <- c("blue", "red", "green", "purple")

	for(i in 1:npic){
		plot <- xyplot(plotdata[,cb[2,i]]~plotdata[,cb[1,i]], groups=group, col=color, pch=16, 
			xlab=paste("LD",cb[1,i],sep=""), ylab=paste("LD",cb[2,i],sep=""))
		print(plot, position=c((i-1)/npic,0,i/npic,1), split=c(1,1,1,1), more=TRUE)	
	}
#	legend("topright",legend=indic,col=color,pch=16,cex=1,xjust=.5, yjust=.5)
#	label <- paste("x: LD1(",69,"%)\ny: LD2(",19.7,"%)\ny: LD3(",11.3,"%)\n")
#	plot4 <- NULL
#	plot4 <- pairs(model$x, main="My Title ", pch=21, bg=color,sub=label)
#	legend("topright",legend=indic,col=color,pch=21,cex=0.5)

	tp2col <- as.vector(group)
	for(i in 1:length(indic)){
		tp2col[tp2col==indic[i]] <- color[i]
	}
	plot3d(plotdata,col=tp2col,xlab="LD1",ylab="LD2",zlab="LD3",type="s",size=1)
#	decorate3d(xlim=c(-100,100), ylim=c(-100,100), zlim=c(-100,100)) 

	dev.off()
}

borutaplot <- function(model, filename){
	png(file=filename, width=2000, height=2000, res=250)
	par(mfrow=c(2,1), mar=c(2,2.5,1,1)+0.1, lwd=0.5, font.axis=2,cex.axis=0.8, mgp=c(2,0.25,0), tck=-0.025)
	plot(model, xlab=NULL, ylab=NULL, ann=FALSE, xaxt='n', outline=FALSE)
	mtext(side=1, text="Attributes", font=2, cex=0.9, line=1)
	mtext(side=2, text="Z-Scores", font=2, cex=0.9, line=1.5)
	plotZHistory(model, xlab=NULL, ann=FALSE)
	mtext(side=1, text="Random Forest Run", font=2, cex=0.9, line=1)
	mtext(side=2, text="Z-Scores", font=2, cex=0.9, line=1.5)
	dev.off()
}

svmRFEplot <- function(svmrfe, filename, label=FALSE){
	featureRank <- svmrfe$featureRank
	kertype <- svmrfe$model[[1]]$kernel
	if(kertype == 0){	title <- "SVM-RFE(linear kernel)"	}
	if(kertype == 1){	title <- "SVM-RFE(polynomial kernel)"	}
	if(kertype == 2){	title <- "SVM-RFE(rbf kernel)"	}

	cex <- 1-log2(nrow(featureRank)/10)/10
	if(nrow(featureRank)<10){	cex <- 1	}

	png(file=filename, width=2000, height=1000, res=250)

		par(mar=c(2.5,2.5,2,1)+0.1, lwd=0.5, font.axis=2, cex.axis=cex, pch=21, cex=cex, mgp=c(3,0.4,0), tck=-0.01)
		plot(rownames(featureRank),as.vector(featureRank[,2]),main=title,ylim=c(0,1),col="blue",bg="blue",ann=FALSE,xaxt='none')
		points(rownames(featureRank),as.vector(featureRank[,3]),col="red",bg="red")
		points(rownames(featureRank),as.vector(featureRank[,4]),col="green",bg="green")

		title(main=title)
		if(label==TRUE){	axis(side=1,at=rownames(featureRank),labels=as.vector(featureRank[,1]),las=2)	}
		mtext(side=1, text="featureRank", font=2, cex=cex+0.1, line=1)
		mtext(side=2, text="Accuracy", font=2, cex=cex+0.1, line=1.5)
		legend("bottomright",legend=c("Cla_Acc","Train_Acc","Test_Acc"),col=c("blue","red","green"),pt.bg=c("blue","red","green"),pch=21,cex=1)

	dev.off()
}

heatmaplot <- function(data, group, ID, method, filename){
	sc <- hclust(dist(data, method=method), "ward")
	ddc <- as.dendrogram(sc)
	colInd <- order.dendrogram(ddc)
	gc <- hclust(dist(t(data), method=method), "ward")
	ddr <- as.dendrogram(gc)
	rowInd <- order.dendrogram(ddr)

###########  define color
	col.set <- c("purple","blue","cyan","yellow")
	png(file=filename, width=3000, height=3000, res=250)
	margins = c(5,18,13)
	colt <- redgreen(75)
	col1 <- col.set
	keysize <- 4.5
	lmat<-rbind(c(0,4,4),c(0,1,1),c(3,2,2),c(0,5,6))

	lwid <- c(keysize,6,6)
	lhei <- c(keysize-1,0.5,40,2.0)

	layout(lmat, widths=lwid, heights=lhei, respect=FALSE)

		########## color bar NormTumor ######################
		par(mar = c(0.5, 0, 0, margins[2]))
		image(matrix(as.numeric(group)[colInd], ncol=1), col=col1, axes=FALSE, xaxt='n', yaxt='n', xlab='', ylab='')
		mtext(side=2, "group", line=0.5, cex=1.3, las=1)

		############# main image ################################
		par(mar = c(margins[1],0, 0, margins[2]))
		image(1:nrow(data), 1:ncol(data), data.matrix(data[colInd,rowInd]), xlim=0.5+c(0,nrow(data)),ylim=0.5+c(0,ncol(data)), axes=FALSE, xlab="", ylab="", col=colt)
		axis(1,at=(1:nrow(data))+0.3, rownames(data)[colInd], las=2, cex.axis=1.0)
		axis(4,at=(1:ncol(data))+0.3, colnames(data)[rowInd], las=1 , cex.axis=1.0)

		########## color bar DMRest ######################
		par(mar=c(margins[1],0,0,0))
		plot(ddr, horiz=TRUE, axes=FALSE, yaxs="i", leaflab="none")
		#############  dendrogram ##########################
		par(mar=c(0,0,0,margins[2]))
		plot(ddc, axes=FALSE, xaxs="i", leaflab="none")
		#############  dendrogram ##########################
		par(mar=c(0,0,0.8,1.5))
		plot(c(-5,5),c(-5,5), axes=FALSE, xaxt="n", yaxt="n", main="", xlab="", ylab="", type="n")
#		rect(-1,-4.5,2,4)
		legend("bottomleft", legend=levels(group), col=col1, pch=15, cex=1.0, bty="n")
		#############  color ##########################
		par(mar=c(2.0, 0, 0, 0.3), cex=1)
		dummy.x <- seq(min(data), max(data), length=length(colt))
		dummy.z <- matrix(dummy.x, ncol=1)
		image(x=dummy.x, y=1, z=dummy.z, xlab="", ylab="", yaxt="n", col=colt)

	dev.off()
}

hclusterplot <- function(dist, group, filename){
	sc <- hclust(as.dist(dist),"ward")
	ddc <- as.dendrogram(sc)
	colInd <- order.dendrogram(ddc)
	col.set <- c("blue","red","green","purple")

	png(filename, width=3000, height=1500, res=250)
	layout(matrix(c(1,2),2,1), widths=c(1), heights=c(6,1), respect=FALSE)
	par(mar=c(0,5,5,1))
	plot(sc,main="H-Clustering",axes=TRUE)
	legend("topright", legend=levels(group), col=col.set, pch=15, cex=2.0, bty="n")
	par(mar=c(3,6,0,2))
	image(matrix(as.numeric(group)[colInd], ncol=1), col=col.set, axes=FALSE, xaxt='n', yaxt='n', xlab='', ylab='')
	mtext(side=2, "group", line=0.5, cex=1.3, las=1)

	dev.off()
}

ENfsplot <- function(elastic_mod, filename){
	png(file=filename, width=2000, height=2000, res=250)
	par(mfrow=c(2,2), mar=c(2,2.5,1,1)+0.1, lwd=0.5, font.axis=2,cex.axis=0.8, mgp=c(1,0.25,0), tck=-0.025)
	plot(elastic_mod[[1]])
	plot(elastic_mod[[2]])
	plot(elastic_mod[[3]])
	plot(elastic_mod[[4]])
	dev.off()
}





