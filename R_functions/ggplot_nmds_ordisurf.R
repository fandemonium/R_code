## input: nmds, factors to group by, ordisurf grid data

ggplot.NMDS.ordisurf<-function(XX, YY, ZZ,COLORS){
        ## XX: metaMDS output
        ## YY: extracted nmds ordisurf data.frame (ordi.sf function)
        ## ZZ: factor for sample grouping
        ## COLORS: custome colors for the grouping

        library(ggplot2)
	
MDS1<-data.frame(scores(XX))$NMDS1
MDS2<-data.frame(scores(XX))$NMDS2
Treatment<-ZZ
Ordi<-YY

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

X1<-ggplot() +
geom_raster(data=Ordi, aes(x=x, y=y, fill = z))+
geom_contour(data=Ordi, aes(x=x, y=y, z=z), color="white")+
scale_fill_gradient(high = "darkgrey", low = "white", breaks=c(min(Ordi$z), (min(Ordi$z)+max(Ordi$z))/2, max(Ordi$z)), labels=c(round(min(Ordi$z),2), round(((min(Ordi$z)+max(Ordi$z))/2), 2), round(max(Ordi$z),2))) +
geom_point(data = NMDS, aes(x=MDS1, y=MDS2, color = Treatment) ,size=3,alpha=0.75) + 
#scale_shape_manual(values=c(21:25)) +
#scale_fill_manual(values=COLORS)+
geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=group), size=2, linetype=5)+
scale_color_manual(values=COLORS, labels=c("No-foam", "Crust", "Foam")) +
theme_classic()+ theme(axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),axis.line.y = element_line(colour = 'black', size=1, linetype='solid'))+
theme(aspect.ratio=1) +
labs(fill = "MPR", x = "MDS1", y= "MDS2") +
theme(axis.text.x=element_text(size=26, face = "bold", color = "black"),axis.text.y=element_text(size=26, face = "bold", color= "black"),axis.title.x=element_text(size=30, face= "bold"),axis.title.y=element_text(size=30, face="bold"))+
guides(color=guide_legend(title=NULL)) +
theme(legend.title = element_text(size = 18), legend.text=element_text(size=16), legend.position = "top", legend.justification=c(0.5,0.5), legend.background = element_rect(fill=(alpha = 0)))


X1    
}
