##plot NMDS using ggplot modified from RJW's code
ggplot.NMDS.ellipse<-function(XX, df, COLORS){
	## df consists of a dataframe with 2 columns (hull factors and ellipse factors)
        library(ggplot2)
MDS1<-data.frame(scores(XX))$NMDS1
MDS2<-data.frame(scores(XX))$NMDS2

## new NMDS dataframe:
#col1: MDS1
#col2: MDS2
#col3: hull factor
#col4: eclipse factor
NMDS<-data.frame(MDS1,MDS2,df)

## Hulls
hull.data<-data.frame()
for (i in levels(NMDS[, 3])){
	temp<-subset(NMDS, NMDS[, 3] == i)
	temp.hull<-temp[chull(temp[, c(1:2)]),]
	hull.data<-rbind(hull.data, temp.hull)
} 

## Ellipses
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
  {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
  }

  df_ell <- data.frame()
  for(g in levels(NMDS[, 4])){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS[, 4]==g,],
                    veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                    ,group=g))
  }

X1<-ggplot() + 
#    geom_point(aes(color = Treatment),size=3,alpha=0.75) +
    geom_polygon(data=hull.data, aes_string(x="MDS1", y="MDS2", fill=colnames(hull.data[, 3]), color = "grey75"), alpha=0.3) +
    geom_path(data=df_ell, aes_string(x="MDS1", y="MDS2",colour=colnames(df_ell[, 4])), size=2, linetype=5)+theme_bw()+theme(aspect.ratio=1)+scale_color_manual(values=COLORS)+theme(axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))+theme(legend.title=element_text(size=15),legend.text=element_text(size=15))
X1    
}

