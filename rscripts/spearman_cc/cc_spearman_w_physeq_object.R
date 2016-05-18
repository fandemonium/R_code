library(phyloseq)
library(Hmisc)
library(plyr)
library(reshape2)
library(igraph)
library(fdrtool)

### input args: 1. phyloseq RDS object; 2. the last column that is a factor; 
args<-commandArgs(TRUE)
input_path<-unlist(strsplit(as.character(args[1]), "/", fixed=T))
full_name<-input_path[length(input_path)]
input_name<-strsplit(as.character(full_name), ".", fixed=T)

# phyloseq RDS object 
physeq<-readRDS(args[1])
print(physeq)

# at otu level:
physeq<-subset_taxa(physeq, domain!="Archaea" & domain!="unclassified_Root")

otu<-data.frame(otu_table(physeq))
si<-data.frame(sample_data(physeq))
tax<-data.frame(tax_table(physeq))
totu<-data.frame(t(otu))

# merging sample information and otu table:
dataset<-merge(si, totu, by.x="SAMPLES", by.y="row.names")
print(dim(dataset))

# co-occurrence by Foaming.Status
#subset the data for a particular treatment
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
sp_melt$spearman_qval<-fdrtool(sp_melt$value, statistic="pvalue", plot=F,verbose=F)$qval
ds_melt$hoeffding_qval<-fdrtool(ds_melt$value, statistic="pvalue", plot=F,verbose=F)$qval
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


#merging together and remove negative rhos
sp_merged<-merge(sp_sub,rhos_melt,by=c("Var1","Var2"))
ds_merged<-merge(ds_sub, ds_melt,by=c("Var1","Var2"))
merged<-merge(sp_merged, ds_merged, by=c("Var1", "Var2"))

# you can write the results out into a flat tab delimited table
write.table(merged, paste(unlist(input_name)[1], "_final_results.txt", sep=""), sep="\t", row.names=F, quote=F)

