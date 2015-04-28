Oread <- read.csv("coverage.check.csv",header=FALSE)
m <- round(mean(Oread[,4]),3)
s <- round(sd(Oread[,4]),3)

png(file="coverage.check.png", width=2000, height=2000, res=250)
hist(Oread[,4],main="Survived Sequence coverage distribution",xlab="Coverage",ylab="Frequency",breaks=100,xlim=c(0.50, 1))
legend("topright",legend=c(paste("mean = ",m,sep=""),paste("SD = ",s,sep="")),pch=21,cex=1)
dev.off()

