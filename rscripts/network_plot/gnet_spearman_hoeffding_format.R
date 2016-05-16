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

strong_results<-subset(final_results, D >= "0.65")
#barn_foaming.rate<-read.delim("barn_foaming_rate.txt", sep="\t", header=T)
barn_foaming.rate<-read.delim(args[2], sep="\t", header=T)
strong_results<-merge(strong_results, barn_foaming.rate[, c("id", "category")], "id")
strong_results$Var1<-paste(strong_results$Var1, strong_results$category, sep="::")
strong_results$Var2<-paste(strong_results$Var2, strong_results$category, sep="::")
strong_results<-data.frame(strong_results[, c("Var1", "Var2", "id")], strong_results[, 4:length(strong_results[1,])])


#phyla<-read.delim("meta_w_genus_information.txt", sep="\t", header=T)
phyla<-read.delim(args[3], sep="\t", header=T)

for (i in unique(strong_results$category)){
        temp<-subset(strong_results, category==i)
	temp.graph<-(graph.edgelist(as.matrix(temp[,c(1,2)]),directed=FALSE))
        E(temp.graph)$weight<-temp$rho
        temp.graph<-simplify(temp.graph)
        gnet<-asNetwork(temp.graph)
        df<-asDF(gnet)
        vs<-df$vertexes
	vs$vertex.names<-data.frame(do.call('rbind', strsplit(as.character(vs$vertex.names), "::", fixed=T)))[, 1]

	## in order to have synchronized color among figures, a set of color has to be generated for everything in group1 first 
	## group2 is essentially the same as group1, so one set of colors is enough here
	for (x in colnames(phyla[, c(3, 8)])){
		colorCount = length(unique(phyla[, x]))
		getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
		colors = getPalette(colorCount)
		## index the desired coloring group, in this case "vs_phyla$group1"
		to_index<-data.frame(unique(phyla[, x]))
		to_index$index<-seq(1, length(to_index[, 1]))
		colnames(to_index)[1]<-x
		phyla_temp<-merge(phyla, to_index, x)
		phyla_temp$index_group<-paste(phyla_temp$index, phyla_temp[, x], sep="_")
		phyla_temp$index_group<-reorder(phyla_temp$index_group, as.numeric(phyla_temp$index))
		names(colors)<-unique(phyla_temp$index_group)
	
	        vs_phyla<-merge(vs, phyla_temp, by.x="vertex.names",by.y="genus")
	        vs_phyla<-arrange(vs_phyla,intergraph_id)
		vs_phyla<-subset(vs_phyla, !grepl("Archaea|Root", domain))
		vs_phyla$domain<-factor(vs_phyla$domain, levels=c("Bacteria", "measurements", "factors"))

		## change gnet vertex.name to color index
		network.vertex.names(gnet)<-vs_phyla$index
		gnet %v% "index_group" <- lapply(vs_phyla[, "index_group"], as.factor)
		gnet %v% "type" <- lapply(vs_phyla[, "domain"], as.character)
#		set.vertex.attribute(gnet, "x", lapply(vs_phyla[, x], as.character))
		set.edge.attribute(gnet, "lty", ifelse(gnet %e% "weight" > 0, 1, 2))

                pdf(paste(unlist(input_name)[1], "_D_0.65", "_foaming_category_", i,"_",x, "_network.pdf", sep=""), height=10, width=18)
                p<-ggnet2(gnet, label=T, size=8, color="index_group", shape="type", shape.palette=c("Bacteria" = 16, "measurements" = 17, "factors"=15), edge.lty="lty", edge.color="black") + scale_color_manual(values=colors, guide=guide_legend(override.aes=list(size=6.5)))
		print(p)
                dev.off()
        }
}

