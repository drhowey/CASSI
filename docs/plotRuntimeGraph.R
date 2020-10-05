

setwd("K:/code/cassi/docs")
 
##graphs for CASSI website
run.times<-c(127, 127, 136, 157, 163, 351, 8914, 27859, 77520)
timePerMil<-(run.times/4498500)*10^6
   timePerMil
   
labels<-c("WZ", "AFE", "AWU", "JE", "PLINK FE", "LR", "PLINK LR", "LR+1 cov", "LR+4 covs")

colours<-c(rainbow(6)[1], rainbow(6)[2], rainbow(6)[3], rainbow(6)[4], "grey", rainbow(6)[5], "grey", "#FF66FFFF", rainbow(6)[6])


png("runtimesCASSI.png", height=600, width=1000, bg="white")
par(mfrow=c(1,2), mar=c(3.1, 4.9, 2.1, 2.1))

barplot(timePerMil[1:6], names=labels[1:6], col=colours[1:6], ylab="seconds per 1 million SNP pairs", ylim=c(0,100), cex.names=1.3, cex.axis=1.6, cex.lab=1.6)

barplot(timePerMil[6:9], names=labels[6:9], col=colours[6:9], ylab="seconds per 1 million SNP pairs", cex.names=1.3,  cex.axis=1.6, cex.lab=1.6)

dev.off()

######################################################################################################
####plot update version with linear regression

setwd("K:/code/cassi/docs")
 
##graphs for CASSI website
run.times<-c(127, 127, 136, 157, 163, 351, 1128, 1810, 4403, 8914, 27859, 77520)
timePerMil<-(run.times/4498500)*10^6
   timePerMil
   
labels<-c("WZ", "AFE", "AWU", "JE", "PLINK FE", "LR", "LIN", "LIN+1 cov", "LIN+4 covs", "PLINK LR", "LR+1 cov", "LR+4 covs")

colours<-c(rainbow(8)[1], rainbow(8)[2], rainbow(8)[3], rainbow(8)[4], "grey", rainbow(8)[5], rainbow(8)[6], "#8044FFFF", rainbow(8)[7], "grey", "#FF44BFFF", rainbow(8)[8])


png("runtimesCASSI2.png", height=600, width=1000, bg="white")
par(mfrow=c(1,2), mar=c(8.1, 4.9, 2.1, 2.1))

barplot(timePerMil[1:7], names.arg=labels[1:7], col=colours[1:7], ylab="seconds per 1 million SNP pairs", ylim=c(0,255), cex.names=1.5, cex.axis=1.6, cex.lab=1.6, las=3)

barplot(timePerMil[6:12], names.arg=labels[6:12], col=colours[6:12], ylab="seconds per 1 million SNP pairs", cex.names=1.5,  cex.axis=1.6, cex.lab=1.6, las=3)

dev.off()

