###############################################
### for cc gnet on subsetted final results: ###
### eg: for i in *_results_*.txt; do Rscript ~/Documents/repos/code/R/gnet_cc_subset.R $i 0.75 ../meta_w_genus_information.txt; done ###
###############################################
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

args<-commandArgs(TRUE)
input_path<-unlist(strsplit(as.character(args[1]), "/", fixed=T))
full_name<-input_path[length(input_path)]
input_name<-strsplit(as.character(full_name), ".", fixed=T)

#final_results<-read.delim("all_493sample_cc_input_final_results.txt", sep="\t", header=T)
#final_results<-read.delim(paste(args[1], "_final_results.txt", sep=""), sep="\t", header=T)
final_results<-read.delim(args[1], sep="\t", header=T)

strong_results<-subset(final_results, rho >= args[2])

#pdf(paste(args[1], "_rho_", args[2], "_histo.pdf", sep=""))
#h<-hist(strong_results$rho)
#h$density<-h$counts/sum(h$counts)*100
#print(plot(h, freq=FALSE, ylab='Percentage'))
#dev.off()

#phyla<-read.delim("meta_w_genus_information.txt", sep="\t", header=T)
phyla<-read.delim(args[3], sep="\t", header=T)

for (i in unique(strong_results$trt)){
        temp.graph<-(graph.edgelist(as.matrix(subset(strong_results, trt==i)[,c(1,2)]),directed=FALSE))
        E(temp.graph)$weight<-subset(strong_results, trt==i)$rho
        temp.graph<-simplify(temp.graph)
        gnet<-asNetwork(temp.graph)
        df<-asDF(gnet)
        vs<-df$vertexes
        vs_phyla<-merge(vs, phyla, by.x="vertex.names",by.y="genus")
        vs_phyla<-arrange(vs_phyla,intergraph_id)
	vs_phyla<-subset(vs_phyla, !grepl("Archaea|Root", domain))
	vs_phyla$domain<-factor(vs_phyla$domain, levels=c("Bacteria", "measurements", "factors"))
	
## index the desired coloring group, in this case "vs_phyla$group1"
	group1<-data.frame(unique(vs_phyla$group1))
	group1$index<-seq(1, length(group1[, 1]))
	colnames(group1)[1]<-"group1"
	vs_phyla<-merge(vs_phyla, group1, "group1")
## change gnet vertex.name to color index
	network.vertex.names(gnet)<-vs_phyla$index
## add color index to legend labeling
	vs_phyla$index_group1<-paste(vs_phyla$index, vs_phyla$group1, sep="_")	

## plot with a for loop if need multiple columns to be plotted:
#        for (x in colnames(vs_phyla[, c(5, 9)])){
                colorCount = length(unique(vs_phyla[, "index_group1"]))
                getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
                colors = getPalette(colorCount)
		names(colors)<-unique(vs_phyla[,"index_group1"])
		gnet %v% "index_group1" <- lapply(vs_phyla[, "index_group1"], as.character)
#		set.vertex.attribute(gnet, "x", lapply(vs_phyla[, x], as.character))
                pdf(paste(unlist(input_name)[1], "_rho_", args[2], "_fs", i, "_", x, "_network.pdf", sep=""), height=10, width=12)

                p<-ggnet2(gnet, label=T, size=8, color="index_group", shape="domain", palette=colors)
		print(p)
                dev.off()
        }
}

