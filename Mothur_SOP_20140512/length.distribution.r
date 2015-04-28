args <- commandArgs(TRUE)
pro <- "Urith"

##############   Original Distribution   ###############
Align <- read.table(paste(pro,".shhh.trim.unique.align.summary",sep=""),header=TRUE)
Align <- read.table(paste(pro,".shhh.trim.unique.summary",sep=""),header=TRUE)
align.start <- rep(Align[1,2],times=Align[1,7])
align.end <- rep(Align[1,3],times=Align[1,7])
align.nbase <- rep(Align[1,4],times=Align[1,7])

for(i in 2:nrow(Align))
{
	tmp <- rep(Align[i,2],times=Align[i,7])
	align.start <- c(align.start,tmp)
	tmp <- rep(Align[i,3],times=Align[i,7])
	align.end <- c(align.end,tmp)
	tmp <- rep(Align[i,4],times=Align[i,7])
	align.nbase <- c(align.nbase,tmp)
}
align.final <- cbind(align.start,align.end,align.nbase)
align.final <- align.final[order(align.final[,3]),]

########   Distribution after Screen   #########
Screen <- read.table(paste(pro,".shhh.trim.unique.good.align.summary",sep=""),header=TRUE)
Screen <- read.table(paste(pro,".shhh.trim.unique.good.summary",sep=""),header=TRUE)
screen.start <- rep(Screen[1,2],times=Screen[1,7])
screen.end <- rep(Screen[1,3],times=Screen[1,7])
screen.nbase <- rep(Screen[1,4],times=Screen[1,7])

for(i in 2:nrow(Screen))
{
	tmp <- rep(Screen[i,2],times=Screen[i,7])
	screen.start <- c(screen.start,tmp)
	tmp <- rep(Screen[i,3],times=Screen[i,7])
	screen.end <- c(screen.end,tmp)
	tmp <- rep(Screen[i,4],times=Screen[i,7])
	screen.nbase <- c(screen.nbase,tmp)
}

screen.final <- cbind(screen.start,screen.end,screen.nbase)
screen.final <- screen.final[order(screen.final[,3]),]

############   Plot Distribution   ##############
png(file=paste(pro,".Seq.summary.png",sep=""), width=2000, height=2000, res=250)
par(mfrow=c(2,2))
hist(align.final[,3],main="Sequence average length before screen",xlab="Nbase",ylab="Frequency",breaks=100,xlim=c(0, 550))
plot(align.final[,3],align.final[,1],main=paste("Total = ",nrow(align.final),sep=""),xlab="Nbase",ylab="Position",
	xlim=c(0,550),ylim=c(0,(max(align.final)+100)),col="blue",pch=21,cex=0.5,bg="blue",lwd=0)
points(align.final[,3],align.final[,2],pch=21,cex=0.5,col="red",bg="red",lwd=0)
legend("topright",legend=c("Start","End"),col=c("blue","red"),pch=21,cex=0.5)

hist(screen.final[,3],main="Sequence average length after screen",xlab="Nbase",ylab="Frequency",breaks=100,xlim=c(0, 550))
plot(screen.final[,3],screen.final[,1],main=paste("Total = ",nrow(screen.final),sep=""),xlab="Nbase",ylab="Position",
	xlim=c(0,550),ylim=c(0,(max(screen.final)+100)),pch=21,cex=0.5,col="blue",bg="blue",lwd=0)
points(screen.final[,3],screen.final[,2],pch=21,cex=0.5,col="red",bg="red",lwd=0)
legend("topright",legend=c("Start","End"),col=c("blue","red"),pch=21,cex=0.5)

dev.off()

