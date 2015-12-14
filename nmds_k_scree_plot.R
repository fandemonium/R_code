library(vegan)

args<-commandArgs(TRUE)
data<-read.delim(args[1], sep="\t", header=T)
data$Classification.Number<-NULL
data<-data[rowSums(is.na(data[, -1:-8]))==0,] ##remove rows with NA's
data.trans<-decostand(data[, -1:-8], "total")
#data.mds<-metaMDS(data.trans, k=3, autotransform=F)

NMDS.scree<-function(x) { #where x is the name of the data frame variable
plot(rep(1,10),replicate(10,metaMDS(x,autotransform=F,k=1)$stress/100),xlim=c(1,nrow(x)),ylim=c(0,0.5),xlab="# of Dimensions",ylab="Stress",main="NMDS stress plot")
 
for (i in 1:(nrow(x)-2)) {
points(rep(i+1,10),replicate(10,metaMDS(x,autotransform=F,k=i+1)$stress/100))
}
}

pdf("test.pdf")
NMDS.scree(data.trans)
dev.off()

