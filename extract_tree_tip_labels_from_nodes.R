## usage: 
## while read node; do Rscript ~/Documents/repos/code/R/extract_tree_tip_labels_from_nodes.R ../uclust_avg_linkage/rpf_tree_uclust.0.65.nwk $node; done < ../rpf_nodes_to_group.txt

library(ape)
library(phylobase)
library(ggplot2)
library(ggtree)

args<-commandArgs(TRUE)

tree<-read.tree(args[1])
nodeid<-as.numeric(args[2])

g1<-as(tree, "phylo4")
node<-subset(g1, node.subtree=nodeid)

out_file<-paste(nodeid, ".tips.txt", sep="")
write.table(tipLabels(node), file=out_file, quote=F, row.names=F, col.names=F)
