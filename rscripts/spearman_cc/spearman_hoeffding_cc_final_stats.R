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

final_results<-read.delim(args[1], sep="\t", header=T)

#barn_foaming.rate<-read.delim("barn_foaming_rate.txt", sep="\t", header=T)
#barn_foaming.rate<-read.delim(args[2], sep="\t", header=T)
#final_results<-merge(final_results, barn_foaming.rate[, c("id", "category")], "id")
#final_results$Var1<-paste(final_results$Var1, final_results$category, sep="::")
#final_results$Var2<-paste(final_results$Var2, final_results$category, sep="::")
#final_results<-data.frame(final_results[, c("Var1", "Var2", "id")], final_results[, 4:length(final_results[1,])])
final_results<-subset(final_results, D >= "0.65")
final_results$Var1<-paste(final_results$Var1, final_results$id, sep="::")
final_results$Var2<-paste(final_results$Var2, final_results$id, sep="::")
final_results<-data.frame(final_results[, c("Var1", "Var2", "id")], final_results[, 4:length(final_results[1,])])

# now we can calculate stats for the network
final_stats<-data.frame()
for(i in unique(final_results$id)){
        temp<-subset(final_results, id==i)
        temp.graph<-(graph.edgelist(as.matrix(temp[,c(1,2)]),directed=FALSE))
        E(temp.graph)$weight<-abs(temp$rho)
        temp.graph<-simplify(temp.graph)
        stats<-data.frame(row.names((as.matrix(igraph::degree(temp.graph,normalized=TRUE)))),(as.matrix(igraph::degree(temp.graph,normalized=TRUE))),(as.matrix(igraph::betweenness(temp.graph))))
        names(stats)<-c("otus","norm_degree","betweenness")
        stats$id <- i
        stats$clustering_coeff<-igraph::transitivity(temp.graph,type="global")
        stats$clustering_coeff_rand<-igraph::transitivity(igraph::erdos.renyi.game(length(V(temp.graph)),length(E(temp.graph)),type="gnm"))
        stats$cluster_ratio<-stats$clustering_coeff/stats$clustering_coeff_rand
        final_stats<-rbind(final_stats,stats)
        print(paste("finished ", i ,sep=": "))
}
# you can write the results out into a flat tab delimited table
write.table(final_stats, paste(unlist(input_name)[1], "_final_stats.txt", sep=""), sep="\t", row.names=F, quote=F) 
