library(Hmisc)
library(plyr)
library(reshape2)
library(igraph)
library(fdrtool)

### input args: 1. phyloseq RDS object; 
args<-commandArgs(TRUE)
input_path<-unlist(strsplit(as.character(args[1]), "/", fixed=T))
full_name<-input_path[length(input_path)]
input_name<-strsplit(as.character(full_name), ".", fixed=T)

# phyloseq RDS object 
physeq<-readRDS(args[1])
print(physeq)

# at genus level:
#physeq<-tax_glom(physeq, "genus")

combined_barn_cc<-data.frame()
for (i in unique(data.frame(sample_data(physeq))$foam.type)){
	tryCatch({
	physeq_sub<-subset_samples(physeq, foam.type==i)
	physeq_sub<-prune_taxa(taxa_sums(physeq_sub)>0, physeq_sub)
	
	otu<-data.frame(otu_table(physeq_sub))
	si<-data.frame(sample_data(physeq_sub))
	totu<-data.frame(t(otu))
	
	# merging sample information and otu table:
	dataset<-merge(si, totu, by.x="SAMPLES", by.y="row.names")
	print(dim(dataset))
	
	to_exclude<-readLines(args[2]) #args[2]: the column names to be excluded
	temp<-dataset[, !names(dataset) %in% to_exclude]
	
	# making an object that has all the results in it (both rho and P values)
	results_sp<-rcorr(as.matrix(temp),type="spearman")
	results_hd<-hoeffd(as.matrix(temp))
	
	#make two seperate objects for p-value and correlation coefficients
	rhos<-results_sp$r
	sp_ps<-results_sp$P
	ds<-results_hd$D
	ds_ps<-results_hd$P
	
	# going to melt these objects to 'long form' where the first two columns make up the pairs of OTUs, I am also removing NA's as they are self-comparisons, not enough data, other bad stuff
	sp_melt<-melt(sp_ps)
	sp_melt<-subset(sp_melt, !is.na(value) & Var1 != Var2 | is.na(value) & Var1 == Var2)
	ds_melt<-melt(ds_ps)
	ds_melt<-subset(ds_melt, !is.na(value) & Var1 != Var2 | is.na(value) & Var1 == Var2)
	
	#creating a qvalue (adjusted pvalue) based on FDR
	sp_melt$spearman_qval<-p.adjust(sp_melt$value, "fdr")
	ds_melt$hoeffding_qval<-p.adjust(ds_melt$value, "fdr")
	#	sp_melt$spearman_qval<-fdrtool(sp_melt$value, statistic="pvalue", plot=F,verbose=F)$qval
	#	ds_melt$hoeffding_qval<-fdrtool(ds_melt$value, statistic="pvalue", plot=F,verbose=F)$qval
	
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
	merged<-subset(merged, spearman_qval < 0.05 | hoeffding_qval < 0.05 | is.na(spearman_qval) | is.na(hoeffding_qval))

	merged$foam.type<-i
	combined_barn_cc<-rbind(combined_barn_cc, merged)
	print(paste("finished ",i,sep=""))
	}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
	
# you can write the results out into a flat tab delimited table
write.table(combined_barn_cc, paste(unlist(input_name)[1], "_core_RawAbun_cc_results_otu.txt", sep=""), sep="\t", row.names=F, quote=F)

