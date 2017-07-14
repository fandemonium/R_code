## ggplot envfit arrows for full nmds (with the first 2 dimensions), with missing data for envfit


mds.envfit.arrows <- function(XX, YY, SAMPLE_ID){
	## XX: full data metaMDS result
	## YY: sample information (data.frame)
	## SAMPLE_ID: colname for the sample ID that could be used to merge
	library(vegan)
	library(ggplot2)
	NMDS<-data.frame(XX$points)[, 1:2]
	NMDS$SAMPLES<-row.names(NMDS)
	NMDS.si<-merge(NMDS, YY, by.x="SAMPLES", by.y=SAMPLE_ID)
	nmds.envfit.df<-data.frame()
	for (i in names(YY)[! names(YY) %in% SAMPLE_ID]){
		NMDS.temp<-NMDS.si[, c("MDS1", "MDS2", i)]
		NMDS.temp<-subset(NMDS.temp, !is.na(NMDS.temp[, i]))
		NMDS.envfit<-envfit(NMDS.temp[, c("MDS1", "MDS2")] ~ NMDS.temp[, i], perm = 999)
		print(NMDS.envfit)
		temp<-data.frame(NMDS.envfit$vectors$arrows * sqrt(NMDS.envfit$vector$r), NMDS.envfit$vectors$pval)
		row.names(temp)<-i
		nmds.envfit.df<-rbind(nmds.envfit.df, temp)
	}
colnames(nmds.envfit.df)[3]<-"pval"
nmds.envfit.df
}
