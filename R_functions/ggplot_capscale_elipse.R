##plot capscale output using ggplot modified from RJW's code
## this capscale output with ellipse function will take experiment factors containing NA's 
## and plot them using one big metaMDS output. 

ggplot.capscale.ellipse<-function(CA, FACTOR, COLORS){
	library(ggplot2)
	library(vegan)

mds <- data.frame(scores(CA)$sites)[, 1:2]

mds$group <- FACTOR

mds.narm <- subset(mds, !is.na(group))

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
  {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
  }

df_ell <- data.frame()
for(g in levels(mds.narm$group)){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(mds.narm[mds.narm$group==g,],
                    veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                    ,group=g))
}

X1<-ggplot(data = mds.narm, aes(MDS1, MDS2)) + geom_point(aes(color = group),size=1.5,alpha=0.75) +
    geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=group), size=2, linetype=5)+
    theme_classic()+ theme(axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),axis.line.y = element_line(colour = 'black', size=1, linetype='solid')) +
    theme(aspect.ratio=1)+
    scale_color_manual(values=COLORS) + #, labels=c("No-foam", "Crust", "Foam")) +
    theme(axis.text.x=element_text(size=26, face = "bold"),axis.text.y=element_text(size=26, face = "bold"),axis.title.x=element_text(size=30, face= "bold"),axis.title.y=element_text(size=30, face="bold"))+theme(legend.title=element_blank(), legend.text=element_text(size=22), legend.position = "top", legend.justification=c(1,0), legend.background = element_rect(fill=(alpha = 0)))
X1    
}

