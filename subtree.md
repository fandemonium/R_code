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
