oldest.mrca<-function(tree,tips){
H<-dist.nodes(tree)
   X<-mrca(tree)
   n<-length(tips)
   nodes<-height<-vector(); k<-1
   for(i in 1:(n-1)) for(j in (i+1):n){
      nodes[k]<-X[tips[i],tips[j]]
      height[k]<-H[match(nodes[k],tree$edge[,1]),1]
      k<-k+1
   }
   z<-match(min(height),height)
   return(nodes[z])
}
