###############################################
### for cc gnet on categorically common final results: ###
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
library(arcdiagram)

args<-commandArgs(TRUE)
input_path<-unlist(strsplit(as.character(args[1]), "/", fixed=T))
full_name<-input_path[length(input_path)]
input_name<-strsplit(as.character(full_name), ".", fixed=T)

#final_results<-read.delim("all_493sample_cc_input_final_results.txt", sep="\t", header=T)
#final_results<-read.delim(paste(args[1], "_final_results.txt", sep=""), sep="\t", header=T)
final_results<-read.delim(args[1], sep="\t", header=T)

#phyla<-read.delim("meta_w_genus_information.txt", sep="\t", header=T)
phyla<-read.delim(args[2], sep="\t", header=T)

#strong_results<-subset(final_results, D >= "0.65")
#strong_results<-subset(final_results, abs(rho) >= "0.55")

#get rid of otu's that are not in the strong_results
#phyla<-phyla[phyla$otu %in% unique(strong_results$Var1), ]
phyla<-phyla[phyla$otu %in% unique(final_results$Var1) | phyla$otu %in% unique(final_results$Var2), ]

#temp<-strong_results
temp<-final_results
net<-(graph.edgelist(as.matrix(temp[,c(1,2)]),directed=FALSE))
#E(net)$rho<-temp$rho
E(net)$wt<-temp$wt
#E(net)$weight<-temp$rho
E(net)$weight<-temp$wt
E(net)$coeff<-temp$coeff
#net<-simplify(net, edge.attr.comb="mean")
## define: node size by degree
net_degree <- igraph::degree(net)

## generate a gradient of color for the net_degree
degree_df<-data.frame(net_degree)
colnames(degree_df)[1]<-"degree"
degree_df$index<-seq(1, length(degree_df[, 1]))
degree_df$nodes<-row.names(degree_df)
uniq_degree<-data.frame(sort(unique(net_degree)))
colnames(uniq_degree)[1]<-"degree"
colfunc <- colorRampPalette(c("grey", "black"))
uniq_degree$cols<-colfunc(length(uniq_degree$degree))
degree_df<-merge(degree_df, uniq_degree, "degree")
degree_df<-arrange(degree_df, degree_df$index)
### define groups to be plotted together
#gclus <- clusters(net)
#gclus.df<-data.frame(gclus$membership)
#degree_df<-merge(degree_df, gclus.df, by.x="nodes", by.y="row.names")
#degree_df<-degree_df[order(degree_df$gclus.membership, degree_df$index),]
##print(degree_df)

## define: positive interactions solid line; negative interactions dash line
#E(net)$lty <- ifelse(E(net)$rho > 0, 1, 2)
E(net)$lty <- ifelse(E(net)$wt < 0, 2, ifelse(E(net)$coeff == "D", 2, 1))

df<-asDF(net)
vs<-df$vertexes
edg<-df$edges
edg$index<-row.names(edg)
## define: positive interactions draw arches above
pos_l <- as.numeric(edg$index[edg$wt > 0])

vs_phyla <- merge(vs, phyla, by.x="name", by.y="otu")
vs_phyla<-arrange(vs_phyla,intergraph_id)
vs_phyla<-subset(vs_phyla, !grepl("Archaea|Root", domain))
vs_phyla$domain<-factor(vs_phyla$domain, levels=c("Bacteria", "measurements", "diet", "factors"))

## define: node shapes
vs_phyla$shape <- vs_phyla$domain
vs_phyla$shape<-gsub("\\bBacteria\\b", 16, vs_phyla$shape)
vs_phyla$shape<-gsub("\\bmeasurements\\b", 17, vs_phyla$shape)
vs_phyla$shape<-gsub("\\bdiet\\b", 15, vs_phyla$shape)
vs_phyla$shape<-gsub("\\bfactors\\b", 18, vs_phyla$shape)

V(net)$shape <- vs_phyla$shape

## rename vertices
#V(net)$name <- as.character(vs_phyla$group2)

## get nodes and edges for basic arcplot frame
net_edges = get.edgelist(net)

gclus = clusters(net)

pdf(paste(unlist(input_name)[1], "_rho_D_ietwork.pdf", sep=""), height=10, width=10)
par(mar=c(25,0,0,0)) 
#p<-arcplot(net_edges, labels=as.character(vs_phyla$group3), col.labels="black", cex.labels=log(net_degree + 1, 5), lty.arcs=E(net)$lty, lwd.arcs= E(net)$wt, col.arcs="black", above = pos_l, col.nodes="black",cex.nodes=log(net_degree + 1, 3), bg.nodes="black", pch.nodes=as.numeric(V(net)$shape), ordering=order(gclus$membership))
p<-arcplot(net_edges, labels=as.character(vs_phyla$group3_otu), col.labels=degree_df$cols, cex.labels=1.5, lty.arcs=E(net)$lty, lwd.arcs= E(net)$wt, col.arcs="black", above = pos_l, col.nodes=degree_df$cols, cex.nodes=2, bg.nodes="black", pch.nodes=as.numeric(V(net)$shape), ordering=order(gclus$membership))
print(p)
dev.off()
