f<-function(x){
	phylum <- row.names(x)
	oid <- x[1]
	l <- rbind(as.vector(unlist(strsplit(oid, ",")), mode="list"))
	as.character(l)
}

## to use:
test2<-apply(test1, 1, f)
