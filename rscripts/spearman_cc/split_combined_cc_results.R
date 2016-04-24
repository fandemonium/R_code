args<-commandArgs(TRUE)
input_path<-unlist(strsplit(as.character(args[1]), "/", fixed=T))
full_name<-input_path[length(input_path)]
input_name<-strsplit(as.character(full_name), ".", fixed=T)

#combined_cc<-read.delim("test/all_valid_samples_min_taxasums_5_min_seq_10k_physeq_combined_barn_cc_results.txt", sep="\t", header=T)
combined_cc<-read.delim(args[1], sep="\t", header=T)

#meta<-read.delim("foaming_status_cc/meta_w_genus_information.txt", sep="\t", header=T)
meta<-read.delim(args[2], sep="\t", header=T)

## separte measurement:measurement, bacteria:bacteria interactions
temp<-merge(combined_cc, meta[, c("genus", "domain")], by.x="Var1", by.y="genus")
temp<-merge(temp, meta[, c("genus", "domain")], by.x="Var2", by.y="genus")
## bacteria to bacteria
bac.bac<-subset(temp, domain.x=="Bacteria" & domain.y=="Bacteria")[, 1:9]
## bacteria to measurements
bac.m<-rbind(subset(temp, domain.x=="Bacteria" & domain.y=="measurements")[, 1:9], subset(temp, domain.y=="Bacteria" & domain.x=="measurements")[, 1:9])

write.table(bac.bac, paste(unlist(input_name)[1], "_combined_cc_bac_bac.txt", sep=""), sep="\t", row.names=F, quote=F)
write.table(bac.m, paste(unlist(input_name)[1], "_combined_cc_bac_measurements.txt", sep=""), sep="\t", row.names=F, quote=F)


