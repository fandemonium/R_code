library(Hmisc)
library(plyr)
library(reshape2)
library(igraph)
library(fdrtool)
library(GGally)
library(intergraph)
library(RColorBrewer)

args<-commandArgs(TRUE)

#final_results<-read.delim("all_493sample_cc_input_final_results.txt", sep="\t", header=T)
final_results<-read.delim(paste(args[1], "_final_results.txt", sep=""), sep="\t", header=T)
#final_stats<-read.delim("all_493sample_cc_input_final_stats.txt", sep="\t", header=T)
final_stats<-read.delim(paste(args[1], "_final_stats.txt", sep=""), sep="\t", header=T)

final_results$trt<-as.factor(final_results$trt)
final_stats$trt<-as.factor(final_stats$trt)

pdf(paste(args[1], "_rho_histo.pdf", sep=""))
print(hist(final_results$rho))
dev.off()

pdf(paste(args[1], "_rho_density.pdf", sep=""))
p<-ggplot(final_results)+geom_density(aes(rho,fill=trt),alpha=0.5)+theme_bw(base_size=17)+theme(aspect.ratio=1)+scale_fill_manual(name="Foaming.Status",values=c("red","black"))
print(p)
dev.off()

pdf(paste(args[1], "_Ndegree_betweenness.pdf", sep=""))
p<-ggplot(final_stats)+geom_point(aes(x=norm_degree,y=betweenness,color=trt),alpha=0.5)+scale_y_log10()+theme_bw(base_size=17)+labs(x="Normalized Degree",y="Betweenness")+theme(aspect.ratio=1)+scale_colour_manual(name="Foaming.Status",values=c("red","black"))
print(p)
dev.off()

## subset results for rho >= 0.7:
strong_results<-subset(final_results, rho >= 0.7)
for (i in colnames(strong_results[, 1:2])){
	strong_results[, i]<-gsub("/| |-", ".", strong_results[, i])
}

pdf(paste(args[1], "_rho_0.7_histo.pdf", sep=""))
print(hist(strong_results$rho))
dev.off()

phyla<-read.delim("meta_w_genus_information.txt", sep="\t", header=T)


for (i in unique(strong_results$trt)){
	temp.graph<-(graph.edgelist(as.matrix(subset(strong_results, trt==i)[,c(1,2)]),directed=FALSE))
	E(temp.graph)$weight<-subset(strong_results, trt==i)$rho
	temp.graph<-simplify(temp.graph)
	gnet<-asNetwork(temp.graph)
	df<-asDF(gnet)
	vs<-df$vertexes
	vs_phyla<-merge(vs, phyla, by.x="vertex.names",by.y="genus")
	vs_phyla<-arrange(vs_phyla,intergraph_id)
	for (x in colnames(vs_phyla[, 4:5])){
		colorCount = length(unique(vs_phyla[, x]))
		getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
		colors = getPalette(colorCount)
		pdf(paste(args[1], "_rho_0.7_fs", i, "_", x, "_network.pdf", sep=""))
		p<-ggnet(gnet, size=0, method="kamadakawaii")+geom_point(aes(colour=vs_phyla[, x]))+scale_colour_manual(values=colors)
		print(p)
		dev.off()
	}
}

