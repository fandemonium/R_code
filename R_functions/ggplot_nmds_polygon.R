##plot NMDS using ggplot modified from RJW's code
ggplot.NMDS.poly<-function(XX, YY, COLORS){
	## YY is the vector with hull factors
	library(ggplot2)
	library(RColorBrewer)
	library(vegan)
	
	MDS1<-data.frame(scores(XX))$NMDS1
	MDS2<-data.frame(scores(XX))$NMDS2
	NMDS<-data.frame(MDS1,MDS2,YY)
	
	## Hulls
	find_hull <- function(df){
		df[chull(df[, 1], df[, 2]), ]
	}
	hull.data <- ddply(NMDS, colnames(NMDS)[3], find_hull)	
	
	X1<-ggplot(NMDS, aes_string(x="MDS1", y="MDS2", color = colnames(NMDS)[3], fill = colnames(NMDS)[3])) + 
	    geom_point(size=3) +
	    geom_polygon(data=hull.data, alpha=0.3) +
	    theme_bw() +
            theme(aspect.ratio=1) +
            scale_color_manual(values=COLORS) +
            scale_fill_manual(values=COLORS) +
            theme(axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20)) +
            theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position = c(0.3, 0.1), legend.justification=c(1, 0))
	X1    
}

