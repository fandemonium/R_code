###############################################
### for cc gnet on subsetted final results: ###
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

args<-commandArgs(TRUE)

#final_results<-read.delim("all_493sample_cc_input_final_results.txt", sep="\t", header=T)
final_results<-read.delim(paste(args[1], "_final_results.txt", sep=""), sep="\t", header=T)

strong_results<-subset(final_results, rho >= args[2])

#pdf(paste(args[1], "_rho_", args[2], "_histo.pdf", sep=""))
#h<-hist(strong_results$rho)
#h$density<-h$counts/sum(h$counts)*100
#print(plot(h, freq=FALSE, ylab='Percentage'))
#dev.off()

#phyla<-read.delim("meta_w_genus_information.txt", sep="\t", header=T)
phyla<-read.delim(args[3], sep="\t", header=T)

for (i in unique(strong_results$Foaming.Status)){
        temp.graph<-(graph.edgelist(as.matrix(subset(strong_results, Foaming.Status==i)[,c(1,2)]),directed=FALSE))
        E(temp.graph)$weight<-subset(strong_results, Foaming.Status==i)$rho
        temp.graph<-simplify(temp.graph)
        gnet<-asNetwork(temp.graph)
        df<-asDF(gnet)
        vs<-df$vertexes
        vs_phyla<-merge(vs, phyla, by.x="vertex.names",by.y="genus")
        vs_phyla<-arrange(vs_phyla,intergraph_id)
	vs_phyla<-subset(vs_phyla, !grepl("Archaea|Root", domain))
	vs_phyla$domain<-factor(vs_phyla$domain, levels=c("Bacteria", "measurements", "factors"))
        for (x in colnames(vs_phyla[, c(5, 9)])){
                colorCount = length(unique(vs_phyla[, x]))
                getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
                colors = getPalette(colorCount)
		names(colors)<-unique(vs_phyla[,x])
#		gnet %v% "x" <- lapply(vs_phyla[, x], as.character)
		set.vertex.attribute(gnet, "x", lapply(vs_phyla[, x], as.character))
                pdf(paste(args[1], "_rho_", args[2], "_fs", i, "_", x, "_network.pdf", sep=""), height=10, width=12)

                p<-ggnet2(gnet, size=5, method="kamadakawaii", color="x", palette=colors)
#			+geom_point(aes(colour=vs_phyla[,x], shape=vs_phyla$domain), size=5)+
#			scale_shape_manual(values=c(16, 8, 12))+
#			scale_colour_manual(name=x, values=colors)+
#			theme_bw()+
#			theme(aspect.ratio=1)+
#			guides(col = guide_legend(ncol = 3))
#			
                print(p)
                dev.off()
        }
}

