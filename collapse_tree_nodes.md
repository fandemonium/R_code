library("ggplot2")
library("ggtree")
tree<-read.tree("uclust_avg_linkage/rpf_tree_uclust.0.65.nwk")

# clades (node id)
#2827
#2672
#3288
#3467
#3192
#2385
#2259
#1994
#1904 
#1805 

```
p <- ggtree(tree, branch.length="none") 
#p <- ggtree(tree)
cp <- p %>% collapse(node=2827) %>% collapse(node=2672) %>% collapse(node=3288) %>% collapse(node=3467) %>% collapse(node=3192) %>% collapse(node=2385) %>% collapse(node=2259) %>% collapse(node=1994) %>% collapse(node=1904 ) %>% collapse(node=1805)

cp1<-cp + geom_point(subset=.(
node == 2827| node==2672| node==3288| node==3467| node==3192| node==2385| node==2259| node==1994| node==1904 | node==1805 
), size=5, shape=23, fill="steelblue")
cp2<-cp1 + geom_text(aes(label=node))
```

## to add hightlight
cp <- p %>% hilight(node=2827, fill="steelblue", alpha=0.6) %>% hilight(node=2672, fill="darkgreen", alpha=0.6) %>% hilight(node=3288, fill="red", alpha=0.6) %>% hilight(node=3467, fill="blue", alpha=0.6) %>% hilight(node=3192, fill="purple", alpha=0.6) %>% hilight(node=2385, fill="yellow", alpha=0.6) %>% hilight(node=2259, fill="grey", alpha=0.6) %>% hilight(node=1994, fill="#999900", alpha=0.6) %>% hilight(node=1904 , fill="#FFCCCC", alpha=0.6) %>% hilight(node=1805, fill="#009999", alpha=0.6)
