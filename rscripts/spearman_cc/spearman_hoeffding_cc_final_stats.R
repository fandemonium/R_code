library(phyloseq)
library(Hmisc)
library(plyr)
library(reshape2)
library(igraph)
library(fdrtool)

source("~/Documents/repos/R_code/R_functions/network_effective_size.R")

### input args: 1. phyloseq RDS object; 2. the last column that is a factor; 
args<-commandArgs(TRUE)
input_path<-unlist(strsplit(as.character(args[1]), "/", fixed=T))
full_name<-input_path[length(input_path)]
input_name<-strsplit(as.character(full_name), ".", fixed=T)

input_results<-read.delim(args[1], sep="\t", header=T)

##barn_foaming.rate<-read.delim("barn_foaming_rate.txt", sep="\t", header=T)
barn_foaming.rate<-read.delim(args[2], sep="\t", header=T)
final_results<-subset(input_results, D >= "0.65")
final_results<-data.frame(final_results[, c("Var1", "Var2", "id")], final_results[, 4:length(final_results[1,])])

# now we can calculate stats for the network
network_clustering<-data.frame()
for(i in unique(final_results$id)){
        temp<-subset(final_results, id==i)
        temp.graph<-(graph.edgelist(as.matrix(temp[,c(1,2)]),directed=FALSE))
        E(temp.graph)$weight<-abs(temp$rho)
        temp.graph<-simplify(temp.graph)
        id <- i
	
	N_nodes<-vcount(temp.graph)
	N_edges<-ecount(temp.graph)
	
	g.components <- clusters(temp.graph)
	N_clusters<-g.components$no
	Max_csize<-max(g.components$csize)
	Max_c_edges<-ecount(induced.subgraph(temp.graph, which(g.components$membership == which.max(g.components$csize))))
        
	g.density<-graph.density(temp.graph)
	g.pathlength.avg<-average.path.length(temp.graph)
	
	betcent<-centralization.betweenness(temp.graph)$centralization
	degcent<-centralization.degree(temp.graph)$centralization
	
	g.mod<-modularity(edge.betweenness.community(temp.graph))
	
	efsize_avg<-mean(effective.size(temp.graph, mode = "all"))
	
	clustering_coeff<-transitivity(temp.graph)
        
        rand<-replicate(1000, {
        	clustering_coeff_rand<-transitivity(erdos.renyi.game(length(V(temp.graph)),length(E(temp.graph)),type="gnm"))
        	})
        rand_avg<-mean(rand)
        rand_ci.975<-qt(.975, df=length(rand)-1)*(sd(rand)/sqrt(length(rand)))
                
        cluster_ratio<-clustering_coeff/rand_avg
        test<-cbind(id, N_nodes, N_edges, N_clusters, Max_csize, Max_c_edges, g.density, g.pathlength.avg, betcent, degcent, g.mod, efsize_avg, clustering_coeff, rand_avg, rand_ci.975, cluster_ratio)
        network_clustering<-rbind(network_clustering, test)
        print(paste("finished ", i ,sep=": "))
}

final_stat<-merge(barn_foaming.rate, network_clustering, "id")

# you can write the results out into a flat tab delimited table
write.table(final_stat, paste(unlist(input_name)[1], "_final_stats.txt", sep=""), sep="\t", row.names=F, quote=F) 
