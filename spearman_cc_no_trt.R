library(Hmisc)
library(plyr)
library(reshape2)
library(igraph)
library(fdrtool)
library(GGally)
library(intergraph)
library(RColorBrewer)
library(ggplot2)
library(network)
library(sna)

## input in terminal: Rscript ~/Documents/repos/code/R/spearman_cc_no_trt.R strain_data.txt 16 0.94
## args[1]: tab delimited data table, located in "/Users/fanyang/Box\ Sync/Projects/Others/Jarboe/strain_data.txt"
## args[2]: last row where the data is factor
## args[3]: rho number indicating strong correlation or to subset

args<-commandArgs(TRUE)
input_name<-strsplit(as.character(args[1]), ".", fixed=T)

# R output df
df<-read.table(args[1],header=T,sep="\t",check.names=F)
dataset<-data.frame(t(df))
#names(dataset)<-dataset["strain", ]
	# making an object that has all the results in it (both rho and P values)
results<-rcorr(as.matrix(dataset[-c(1:args[2]), ]),type="spearman")
#results<-rcorr(as.matrix(dataset),type="spearman")
	#make two seperate objects for p-value and correlation coefficients
rhos<-results$r
ps<-results$P
	# going to melt these objects to 'long form' where the first two columns make up the pairs of OTUs, I am also removing NA's as they are self-comparisons, not enough data, other bad stuff
ps_melt<-na.omit(melt(ps))
	#creating a qvalue based on FDR
#ps_melt$qval<-fdrtool(ps_melt$value, statistic="pvalue", plot=F,verbose=F)$qval
	# in case of too few points, use p.adjust function instead
ps_melt$qval<-p.adjust(ps_melt$value, "fdr")
	#making column names more relevant
	
names(ps_melt)[3]<-"pval"
	# if you are of the opinion that it is a good idea to subset your network based on adjusted P-values (qval in this case), you can then subset here
ps_sub<-subset(ps_melt, qval < 0.05)
	# now melting the rhos, note the similarity between ps_melt and rhos_melt
rhos_melt<-na.omit(melt(rhos))
names(rhos_melt)[3]<-"rho"
	#merging together and remove negative rhos
final_results<-merge(ps_sub,subset(rhos_melt, rho > 0),by=c("Var1","Var2"))
write.table(final_results, paste(unlist(input_name)[1], "_final_results.txt", sep=""), sep="\t", row.names=F, quote=F)

# now we can calculate stats for the network
temp.graph<-(graph.edgelist(as.matrix(final_results[,c(1,2)]),directed=FALSE))
E(temp.graph)$weight<-final_results$rho
temp.graph<-simplify(temp.graph)
final_stats<-data.frame(row.names((as.matrix(igraph::degree(temp.graph,normalized=TRUE)))),(as.matrix(igraph::degree(temp.graph,normalized=TRUE))),(as.matrix(igraph::betweenness(temp.graph))))
names(final_stats)<-c("otus","norm_degree","betweenness")
final_stats$clustering_coeff<-igraph::transitivity(temp.graph,type="global")
final_stats$clustering_coeff_rand<-igraph::transitivity(igraph::erdos.renyi.game(length(V(temp.graph)),length(E(temp.graph)),type="gnm"))
final_stats$cluster_ratio<-final_stats$clustering_coeff/final_stats$clustering_coeff_rand
write.table(final_stats, paste(unlist(input_name)[1], "_final_stats.txt", sep=""), sep="\t", row.names=F, quote=F)

# evaluate histogram
# subset for strong cooccurring items
strong_results<-subset(final_results, rho >= args[3])

pdf(paste(unlist(input_name)[1], "_rho_", args[3], "_histo.pdf", sep=""))
h<-hist(strong_results$rho)
h$density<-h$counts/sum(h$counts)*100
print(plot(h, freq=FALSE, ylab='Percentage'))
dev.off()

temp.graph<-(graph.edgelist(as.matrix(strong_results[,c(1,2)]),directed=FALSE))
E(temp.graph)$weight<-strong_results$rho
temp.graph<-simplify(temp.graph)
gnet<-asNetwork(temp.graph)
df<-asDF(gnet)
vs<-df$vertexes
vs_phyla<-merge(vs, data.frame(t(dataset[1:16, ])), by.x="vertex.names",by.y="row.names")
vs_phyla<-arrange(vs_phyla,intergraph_id)

# creating a new group
cols<-colnames(vs_phyla)[7:19]
vs_phyla$group<-apply(vs_phyla[, cols], 1, paste, collapse="_")
# color network by each column
for (x in colnames(vs_phyla[, c(4, 7:20)])){
	colorCount = length(unique(vs_phyla[, x]))
        getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
        colors = getPalette(colorCount)
	names(colors)<-unique(vs_phyla[,x])
#	gnet %v% "x" <- lapply(vs_phyla[, x], as.character)
	set.vertex.attribute(gnet, "x", lapply(vs_phyla[, x], as.character))
        pdf(paste(args[1], "_rho_", args[3], "_", x, "_network.pdf", sep=""), height=10, width=12)

        p<-ggnet2(gnet, size=5, method="kamadakawaii", color="x", palette=colors)
        print(p)
        dev.off()
}


