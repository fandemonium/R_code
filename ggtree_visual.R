library("ape")
library("ggplot2")
library("ggtree")


args<-commandArgs(TRUE)

tree <- read.tree(args[1])
pdf(args[2])
ggtree(tree)+geom_text(aes(label=node, size=0.3))
#ggtree(tree, branch.length="none")+geom_text(aes(label=node, size=0.5))
#gzoom(tree, grep("AmpC", tree$tip.label))
dev.off()

