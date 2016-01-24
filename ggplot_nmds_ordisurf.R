## input: nmds, factors to group by, ordisurf grid data

ggplot.NMDS.ordisurf<-function(XX,ZZ,COLORS, YY){
        library(ggplot2)
MDS1<-data.frame(scores(XX))$NMDS1
MDS2<-data.frame(scores(XX))$NMDS2
Treatment<-ZZ
Ordi<-YY

NMDS<-data.frame(MDS1,MDS2,Treatment)

#X1<-ggplot(data = NMDS, aes(x=MDS1, y=MDS2)) + 
#stat_contour(data = Ordi, aes(x = x, y = y, z = z, color = rev(..level..), binwidth = 2)) +
#geom_point(aes(shape = Treatment),size=3,alpha=0.75) +
#theme_bw() +
#theme(aspect.ratio=1) +
#labs(color = "Substrate Density") +
#scale_color_gradient(high = "darkgreen", low = "darkolivegreen1") +
#theme(axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20)) +
#theme(legend.title=element_text(size=15),legend.text=element_text(size=15))
#
X1<-ggplot() +
geom_point(data = NMDS, aes(x=MDS1, y=MDS2, fill = Treatment, shape = Treatment) ,size=3,alpha=0.75) + 
scale_shape_manual(values=c(21:25)) +
scale_fill_manual(values=COLORS)+
#stat_contour(data = Ordi, aes(x = x, y = y, z = z, color = rev(..level..), binwidth = 2)) +
stat_contour(data = Ordi, aes(x = x, y = y, z = z, color = ..level.., binwidth = 2)) +
scale_color_gradient(high = "darkgreen", low = "darkolivegreen1") +
theme_bw() +
theme(aspect.ratio=1) +
labs(color = "Substrate Density") +
theme(axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20)) +
theme(legend.title=element_text(size=15),legend.text=element_text(size=15))

X1    
}
