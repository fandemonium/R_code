library(reshape2)
library(plyr)

args<-commandArgs(TRUE)

## /mnt/scratch/yangfan1/pitfoaming/Processed_16s/3_cdhit_clustering/R/all_16s_donotuseRemoved_no_si_psmelt_table.txt
data<-read.delim(args[1], sep="\t", header=T, fill=T)
## /mnt/scratch/yangfan1/pitfoaming/Processed_16s/3_cdhit_clustering/R/fixed_inputs/sample_information_donotuseRemoved_rearranged.txt
si<-read.delim(args[2], sep="\t", header=T, fill=T)

temp<-ddply(data, .(Sample, domain, phylum, class, order, family, genus), summarize, abundance=sum(Abundance))
colnames(temp)[1]<-"SAMPLES"
genus_out<-dcast(temp, SAMPLES ~ genus, value.var="abundance")

data_full<-merge(si, genus_out, "SAMPLES")
print(dim(data_full))
write.table(data_full, "all_493sample_cc_input.txt", sep="\t", row.names=F, quote=F)

ntrt<-subset(data_full, Treated.Status=="0")
write.table(ntrt, "all_sample_notrt_cc_input.txt", sep="\t", row.names=F, quote=F)

cfVcnf<-subset(data_full, grepl("CF|CNF", Name))
write.table(cfVcnf, "cfVcnf_cc_input.txt", sep="\t", row.names=F, quote=F)

cfVcnf_ntrt<-subset(cfVcnf, Treated.Status=="0")
write.table(cfVcnf_ntrt, "cfVcnf_notrt_cc_input.txt", sep="\t", row.names=F, quote=F)
 

