library(phyloseq)
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
physeq<-subset_taxa(physeq, domain!="Archaea" & domain!="unclassified_Root")

combined_barn_cc<-data.frame()
for (i in unique(data.frame(sample_data(physeq))$id)){
	tryCatch({
	physeq_sub<-subset_samples(physeq, id==i)
	physeq_sub<-prune_taxa(taxa_sums(physeq_sub)>0, physeq_sub)

	otu<-data.frame(otu_table(physeq_sub))
	si<-data.frame(sample_data(physeq_sub))
#	tax<-data.frame(tax_table(physeq_sub))
#	row.names(otu)<-tax$genus
	totu<-data.frame(t(otu))
	
	# merging sample information and otu table:
	dataset<-merge(si, totu, by.x="SAMPLES", by.y="row.names")
	print(dim(dataset))
	
	temp<-dataset[,-c(1:21, 108:120, 122)]
	# making an object that has all the results in it (both rho and P values)
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
	#	sp_melt$spearman_qval<-fdrtool(sp_melt$value, statistic="pvalue", plot=F,verbose=F)$qval
	#	ds_melt$hoeffding_qval<-fdrtool(ds_melt$value, statistic="pvalue", plot=F,verbose=F)$qval
	
	#making column names more relevant
	names(sp_melt)[3]<-"spearman_pval"
	names(ds_melt)[3]<-"hoeffding_pval"
	
	# if you are of the opinion that it is a good idea to subset your network based on adjusted P-values (qval in this case), you can then subset here
	sp_sub<-subset(sp_melt, spearman_qval < 0.05)
	ds_sub<-subset(ds_melt, hoeffding_qval < 0.05)
	
	# now melting the rhos, note the similarity between ps_melt and rhos_melt
	rhos_melt<-na.omit(melt(rhos))
	ds_melt<-na.omit(melt(ds))
	
	names(rhos_melt)[3]<-"rho"
	names(ds_melt)[3]<-"D"
	
	#merging together 
	sp_merged<-merge(sp_sub,rhos_melt,by=c("Var1","Var2"))
	ds_merged<-merge(ds_sub, ds_melt,by=c("Var1","Var2"))
	merged<-merge(sp_merged, ds_merged, by=c("Var1", "Var2"))
	
	merged$id<-i
	combined_barn_cc<-rbind(combined_barn_cc, merged)
	print(paste("finished ",i,sep=""))
	}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
	
# you can write the results out into a flat tab delimited table
write.table(combined_barn_cc, paste(unlist(input_name)[1], "_combined_barn_RawAbun_cc_results_otu.txt", sep=""), sep="\t", row.names=F, quote=F)

### now we can calculate stats for the network
##final_stats<-data.frame()
##for(i in 1:length(unique(final_results$trt))){
##	temp<-subset(final_results, trt==as.vector(unique(final_results$trt))[i])
##	temp.graph<-(graph.edgelist(as.matrix(temp[,c(1,2)]),directed=FALSE))
##	E(temp.graph)$weight<-temp$rho
##	temp.graph<-simplify(temp.graph)
##	stats<-data.frame(row.names((as.matrix(igraph::degree(temp.graph,normalized=TRUE)))),(as.matrix(igraph::degree(temp.graph,normalized=TRUE))),(as.matrix(igraph::betweenness(temp.graph))))
##	names(stats)<-c("otus","norm_degree","betweenness")
##	stats$trt<-as.vector(unique(final_results$trt))[i]
##	stats$clustering_coeff<-igraph::transitivity(temp.graph,type="global")
##	stats$clustering_coeff_rand<-igraph::transitivity(igraph::erdos.renyi.game(length(V(temp.graph)),length(E(temp.graph)),type="gnm"))
##	stats$cluster_ratio<-stats$clustering_coeff/stats$clustering_coeff_rand
##	final_stats<-rbind(final_stats,stats)
##	print(paste("finished ",as.vector(unique(final_results$trt))[i],sep=""))
##}
### you can write the results out into a flat tab delimited table
##write.table(final_stats, paste(unlist(input_name)[1], "_final_stats.txt", sep=""), sep="\t", row.names=F, quote=F)
#
###combined_cc<-read.delim("test/all_valid_samples_min_taxasums_5_min_seq_10k_physeq_combined_barn_cc_results.txt", sep="\t", header=T)
##combined_cc<-read.delim(args[1], sep="\t", header=T)
#
##meta<-read.delim("foaming_status_cc/meta_w_genus_information.txt", sep="\t", header=T)
#meta<-read.delim(args[2], sep="\t", header=T)
#
### separte measurement:measurement, bacteria:bacteria interactions
#temp<-merge(combined_barn_cc, meta[, c("genus", "domain")], by.x="Var1", by.y="genus")
#temp<-merge(temp, meta[, c("genus", "domain")], by.x="Var2", by.y="genus")
### bacteria to bacteria
#bac.bac<-subset(temp, domain.x=="Bacteria" & domain.y=="Bacteria")[, 1:9]
### bacteria to measurements
#bac.m<-rbind(subset(temp, domain.x=="Bacteria" & domain.y=="measurements")[, 1:9], subset(temp, domain.y=="Bacteria" & domain.x=="measurements")[, 1:9])
#
#write.table(bac.bac, paste(unlist(input_name)[1], "_combined_RelaAbun_cc_bac_bac.txt", sep=""), sep="\t", row.names=F, quote=F)
#write.table(bac.m, paste(unlist(input_name)[1], "_combined_RelaAbun_cc_bac_measurements.txt", sep=""), sep="\t", row.names=F, quote=F)
#
