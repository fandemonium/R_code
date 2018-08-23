# base functions of calculated spearman and hoeffding corelations
# works on numberic dataframe
# p value will be adjusted by fdr
# final results will be subsetted for spearman significant or hoeffding significance

cc_spearman_hoeff <- function(df){
	library(Hmisc)
	library(plyr)
	library(reshape2)
	
	temp <- df[ complete.cases(df), ]
	results_sp<-rcorr(as.matrix(temp),type="spearman")
	results_hd<-hoeffd(as.matrix(temp))
	
	#make two seperate objects for p-value and correlation coefficients
	rhos<-results_sp$r
	sp_ps<-results_sp$P
	ds<-results_hd$D
	ds_ps<-results_hd$P
	
	# going to melt these objects to 'long form' where the first two columns make up the pairs of OTUs, I am also removing NA's as they are self-comparisons, not enough data, other bad stuff
	sp_melt<-na.omit(melt(sp_ps))
	ds_melt<-na.omit(melt(ds_ps))
	
	#creating a qvalue (adjusted pvalue) based on FDR
	sp_melt$spearman_qval<-p.adjust(sp_melt$value, "fdr")
	ds_melt$hoeffding_qval<-p.adjust(ds_melt$value, "fdr")
	
	#making column names more relevant
	names(sp_melt)[3]<-"spearman_pval"
	names(ds_melt)[3]<-"hoeffding_pval"
	
	# now melting the rhos, note the similarity between ps_melt and rhos_melt
	rho.melt<-na.omit(melt(rhos))
	d.melt<-na.omit(melt(ds))
	
	names(rho.melt)[3]<-"rho"
	names(d.melt)[3]<-"D"
	
	#merging together  then subset
	sp_merged<-merge(sp_melt,rho.melt,by=c("Var1","Var2"))
	ds_merged<-merge(ds_melt, d.melt,by=c("Var1","Var2"))
	merged<-merge(sp_merged, ds_merged, by=c("Var1", "Var2"))
	merged<-subset(merged, spearman_qval < 0.05 | hoeffding_qval < 0.05)
}


	

