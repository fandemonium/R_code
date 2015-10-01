library(ape)
library(phylobase)
library(ggplot2)
library(ggtree)

tree<-read.tree("uclust_avg_linkage/rpf_tree_uclust.0.65.nwk")

g1<-as(tree, "phylo4")

## when plotting the subtree, the node id will be displayed, which is different from the node label. one things to make sure you are looking at the right thing is do this after you get your tree set up: ##
nodeLabels(geospiza) <- paste("N", nodeId(geospiza, "internal"), sep="")

ggtree(subset(g1, node.subtree=2827)) #node.subtree: use node id
# can also use plot(subset(....)), but ggtree is much nicer
#ggtree(subset(g1, node.subtree=2827), branch.length='none') + geom_text(aes(label=node, size=0.5))

n2827<-subset(g1, node.subtree=2827)

## to plot the alignment of fasta with the tree, simply load the alignment fa by typing in the path
msaplot(ggtree(n2827), "../../Reference Seqs/rpfHits_jones/rpfHits_hmm.faa")

## to get tip labels:
tipLabels(n2827)
## to get node id (non single tip node):
nodeId(n2827, "internal")

## phylobase also display tree information like a table (but it's not). To capture the df like outpt:
capture.output(n2827, file="test_R_output.txt")
