library(reshape2)

args<-commandArgs(TRUE)

#df <- read.delim("strain_data.txt", sep="\t", header=T)
df <- read.delim(args[1], sep="\t", header=T)

binary<-data.frame(df$strain)
names(binary)<-"strain"
for (i in colnames(df[, 4:16])){
	df1<-data.frame(df$strain, df[, i])
	names(df1)<-c("strain", "X")
	temp<-dcast(df1, strain ~ X, fill=0, fun.aggregate=function(x) 1, value.var="X")
	l<-colnames(temp[, 2:length(temp[1,])])
	new<-apply(as.matrix(data.frame(i,l)), 1, paste, collapse=".")
	names(temp)<-c("strain", new)
	binary<-merge(binary, temp, "strain")
}
binary<-merge(binary, df[, c(1, 17:26)], "strain")
write.table(binary, args[2], sep="\t", row.names=F, quote=F)
