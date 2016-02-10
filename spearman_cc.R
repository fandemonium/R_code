library(Hmisc)
library(plyr)
library(reshape2)
library(igraph)
library(fdrtool)

args<-commandArgs(TRUE)
input_name<-strsplit(as.character(args[1]), ".", fixed=T)

# R output df
dataset<-read.table(args[1],header=T,sep="\t",check.names=F)

# eg, Foaming.Status
treatments<-as.vector(unique(dataset$Foaming.Status))
final_results<-data.frame()
for(i in 1:length(treatments)){
	#subset the data for a particular treatment
	temp<-subset(dataset, Foaming.Status==treatments[i])
	# making an object that has all the results in it (both rho and P values)
	results<-rcorr(as.matrix(temp[,-c(1:args[2])]),type="spearman")
	#make two seperate objects for p-value and correlation coefficients
	rhos<-results$r
	ps<-results$P
	# going to melt these objects to 'long form' where the first two columns make up the pairs of OTUs, I am also removing NA's as they are self-comparisons, not enough data, other bad stuff
	ps_melt<-na.omit(melt(ps))
	#creating a qvalue based on FDR
	ps_melt$qval<-fdrtool(ps_melt$value, statistic="pvalue", plot=F,verbose=F)$qval
	#making column names more relevant
	
	names(ps_melt)[3]<-"pval"
	# if you are of the opinion that it is a good idea to subset your network based on adjusted P-values (qval in this case), you can then subset here
	ps_sub<-subset(ps_melt, qval < 0.05)
	# now melting the rhos, note the similarity between ps_melt and rhos_melt
	rhos_melt<-na.omit(melt(rhos))
	names(rhos_melt)[3]<-"rho"
	#merging together and remove negative rhos
	merged<-merge(ps_sub,subset(rhos_melt, rho > 0),by=c("Var1","Var2"))
	merged$trt<-treatments[i]
	final_results<-rbind(final_results, merged)
	print(paste("finished ",treatments[i],sep=""))
}
write.table(final_results, paste(unlist(input_name)[1], "_final_results.txt", sep=""), sep="\t", row.names=F, quote=F)

# now we can calculate stats for the network
final_stats<-data.frame()
for(i in 1:length(unique(final_results$trt))){
	temp<-subset(final_results, trt==as.vector(unique(final_results$trt))[i])
	temp.graph<-(graph.edgelist(as.matrix(temp[,c(1,2)]),directed=FALSE))
	E(temp.graph)$weight<-temp$rho
	temp.graph<-simplify(temp.graph)
	stats<-data.frame(row.names((as.matrix(igraph::degree(temp.graph,normalized=TRUE)))),(as.matrix(igraph::degree(temp.graph,normalized=TRUE))),(as.matrix(igraph::betweenness(temp.graph))))
	names(stats)<-c("otus","norm_degree","betweenness")
	stats$trt<-as.vector(unique(final_results$trt))[i]
	stats$clustering_coeff<-igraph::transitivity(temp.graph,type="global")
	stats$clustering_coeff_rand<-igraph::transitivity(igraph::erdos.renyi.game(length(V(temp.graph)),length(E(temp.graph)),type="gnm"))
	stats$cluster_ratio<-stats$clustering_coeff/stats$clustering_coeff_rand
	final_stats<-rbind(final_stats,stats)
	print(paste("finished ",as.vector(unique(final_results$trt))[i],sep=""))
}
write.table(final_stats, paste(unlist(input_name)[1], "_final_stats.txt", sep=""), sep="\t", row.names=F, quote=F)

