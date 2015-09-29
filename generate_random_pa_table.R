#This script generates a random trait table for consentrait using an existing trait table's row names as template
## usage: Rscript ~/Documents/repos/code/R/generate_random_pa_table.R <exsiting trait table> <output random trait table>

## get row names from observed trait table ##
args <- commandArgs(TRUE)
table<-read.delim(args[1], sep="\t", header=F)

## generate an empty 10x10 matrix ##
r <- length(table$V1)
c <- 10 
m0<-matrix(0, r, c)
m1<-apply(m0, c(1, 2), function(x) sample(c(0,1), 1))
df1<-as.data.frame(m1)
class(df1)

rownames(df1)<-table$V1

write.table(df1, args[2], sep="\t", quote=F, row.names=T, col.names=F)
