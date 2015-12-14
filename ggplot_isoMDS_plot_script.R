## script is written for pitfoaming arisa data currently

library(vegan)
library(MASS)
library(ggplot2)
library(RColorBrewer)

## RyanJW's nmds plot function:
##NMDS for 16S
ggplot.NMDS<-function(XX,ZZ,COLORS){
	library(ggplot2)
MDS1<-data.frame(scores(XX))$Dim1
MDS2<-data.frame(scores(XX))$Dim2
Treatment<-ZZ

NMDS<-data.frame(MDS1,MDS2,Treatment)

NMDS.mean=aggregate(NMDS[,1:2],list(group=Treatment),mean)

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
  {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
  }

  df_ell <- data.frame()
  for(g in levels(NMDS$Treatment)){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$Treatment==g,],
                    veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                    ,group=g))
  }

X1<-ggplot(data = NMDS, aes(MDS1, MDS2)) + geom_point(aes(color = Treatment, shape = Treatment),size=3,alpha=0.75) +
    geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=group), size=2, linetype=5)+theme_bw()+theme(aspect.ratio=1)+scale_color_manual(values=COLORS)+theme(axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))+theme(legend.title=element_text(size=15),legend.text=element_text(size=15))
X1    
}


args<-commandArgs(TRUE)

## Angela's arisa data,  read in needs to have samples as rows, otu's as columns. 
data<-read.delim(args[1], sep="\t", header=T)
data$Classification.Number<-NULL
data<-data[rowSums(is.na(data[, -1:-8]))==0,] ##remove rows with NA's
data.trans<-decostand(data[, -1:-8], "total")
data.dis<-vegdist(data.trans)
data.mds<-isoMDS(data.dis, k=3)

for (i in names(data[, 3:8])){
	print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
	colorCount = length(unique(data[, i]))
	getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
	colors = getPalette(colorCount)
	pdf(paste("arisa_nmds_",i,".pdf", sep=""))
	ggplot.NMDS(data.mds, data[, i], colors)
	dev.off()
	adonis(data.trans~data[,i])
} 
