library(phyloseq)
library(vegan)

args <- commandArgs(TRUE)

otu <- read.delim(args[1])
print("the dimension of the OTU table is: ")
dim(otu)

tax <- read.delim(args[2])
print("the dimension of the taxa table is: ")
dim(tax)

row.names(otu)<-otu$OTUS
otu$OTUS <- NULL

row.names(tax)<-tax$OTUS

data.phy<-phyloseq(otu_table(as.matrix(otu),taxa_are_rows=T), tax_table(as.matrix(tax)))
