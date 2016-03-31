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

	## in order to have synchronized color among figures, a set of color has to be generated for everything in group1 first 
	## group2 is essentially the same as group1, so one set of colors is enough here
	for (x in colnames(phyla[, c(7, 8)])){
		colorCount = length(unique(phyla[, x]))
		getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
		colors = getPalette(colorCount)
		## index the desired coloring group, in this case "vs_phyla$group1"
		to_index<-data.frame(unique(phyla[, x]))
		to_index$index<-seq(1, length(to_index[, 1]))
		colnames(to_index)[1]<-x
		phyla_temp<-merge(phyla, to_index, x)
		phyla_temp$index_group<-paste(phyla_temp$index, phyla_temp[, x], sep="_")
		names(colors)<-unique(phyla_temp$index_group)
	
	        vs_phyla<-merge(vs, phyla_temp, by.x="vertex.names",by.y="genus")
	        vs_phyla<-arrange(vs_phyla,intergraph_id)
		vs_phyla<-subset(vs_phyla, !grepl("Archaea|Root", domain))
		vs_phyla$domain<-factor(vs_phyla$domain, levels=c("Bacteria", "measurements", "factors"))

		## change gnet vertex.name to color index
		network.vertex.names(gnet)<-vs_phyla$index
		gnet %v% "index_group" <- lapply(vs_phyla[, "index_group"], as.character)
#		set.vertex.attribute(gnet, "x", lapply(vs_phyla[, x], as.character))

                pdf(paste(unlist(input_name)[1], "_rho_", args[2], "_fs", i,"_",x, "_network.pdf", sep=""), height=10, width=12)
                p<-ggnet2(gnet, label=T, size=8, color="index_group", shape="domain", palette=colors)
		print(p)
                dev.off()
        }
}

