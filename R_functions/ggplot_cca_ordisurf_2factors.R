## input: cca, factors to group by, ordisurf grid data

ggplot.cca.ordisurf.2f<-function(CCA, DF_2F, COLORS, ORDI, GRADIENT){
scrs <- data.frame(scores(CCA, display="sites"))
CAP1<-scrs$CAP1
CAP2<-scrs$CAP2

cca.df<-data.frame(CAP1,CAP2,DF_2F)
print(names(cca.df))

X1<-ggplot() +
geom_point(data = cca.df, aes_string(x="CAP1", y="CAP2", fill = colnames(cca.df)[3], shape = colnames(cca.df)[4]) ,size=3,alpha=0.75) + 
scale_shape_manual(values=c(21:25)) +
scale_fill_manual(values=c(COLORS), guide = guide_legend(override.aes = list(shape = 23)))+
stat_contour(data = ORDI, aes(x = x, y = y, z = z, color = ..level.., binwidth = 2)) +
scale_color_gradient(high = "darkgreen", low = "darkolivegreen1") +
theme_bw() +
theme(aspect.ratio=1) +
labs(color = GRADIENT) +
theme(axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20)) +
theme(legend.title=element_text(size=15),legend.text=element_text(size=15))

X1    
}
