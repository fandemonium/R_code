##plot NMDS using ggplot modified from RJW's code
ggplot.NMDS.ellipse.arrow<-function(XX,YY, ZZ,COLORS){
	## XX: metaMDS output
	## YY: extracted nmds envfit data.frame
	## ZZ: factor for sample grouping
	## COLORS: custome colors for the grouping

        library(ggplot2)
MDS1<-data.frame(scores(XX))$NMDS1
MDS2<-data.frame(scores(XX))$NMDS2
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

X1<-ggplot(data = NMDS, aes(MDS1, MDS2)) + geom_point(aes(color = Treatment),size=1.5,alpha=0.75) +
    geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=group), size=2, linetype=5)+
    theme_classic()+ theme(axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),axis.line.y = element_line(colour = 'black', size=1, linetype='solid')) +
    theme(aspect.ratio=1)+
    geom_segment(data = YY, aes(x=0, xend=MDS1, y=0, yend=MDS2), arrow = arrow(length=unit(0.5, "cm")), color = "black", size = 1) +
    geom_text(data=YY, aes(x=MDS1, y=MDS2, label= row.names(YY)), nudge_x = 0.02, nudge_y = 0.01)+
    scale_color_manual(values=COLORS, labels=c("No-foam", "Crust", "Foam")) +
    theme(axis.text.x=element_text(size=26, face = "bold"),axis.text.y=element_text(size=26, face = "bold"),axis.title.x=element_text(size=30, face= "bold"),axis.title.y=element_text(size=30, face="bold"))+theme(legend.title=element_blank(), legend.text=element_text(size=22), legend.position = "top", legend.justification=c(1,0), legend.background = element_rect(fill=(alpha = 0)))
X1    
}

