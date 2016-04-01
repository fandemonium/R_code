args<-commandArgs(TRUE)
input_path<-unlist(strsplit(as.character(args[1]), "/", fixed=T))
full_name<-input_path[length(input_path)]
input_name<-strsplit(as.character(full_name), ".", fixed=T)

#final_results<-read.delim("all_493sample_cc_input_final_results.txt", sep="\t", header=T)
final_results<-read.delim(args[1], sep="\t", header=T)

strong_results<-subset(final_results, rho >= args[2])

pdf(paste(unlist(input_name)[1], "_rho_", args[2], "_histo.pdf", sep=""))
h<-hist(strong_results$rho)
h$density<-h$counts/sum(h$counts)*100
print(plot(h, freq=FALSE, ylab='Percentage'))
dev.off()

