library(sp)

args <- commandArgs(TRUE)

#both input tables need to have column with names "longitude" and "latitude"
given_table <- read.delim(args[1])
region_table <- read.delim(args[2]) #to overlay against, can have more than just the boundary points

coordinates(given_table) <- ~ longitude + latitude
uniq_regions <- unique(region_table[, args[3]]) #the single column name that's being used to define differen sub regions

overlayed <- data.frame()
for (i in uniq_regions){
	region <- subset(region_table, region_table[, args[3]] == i)
	ch <-chull(region[, c("longitude", "latitude")])
	region.ch <- region[c(ch, ch[1]), ]
	pol <- SpatialPolygons(list(Polygons(list(Polygon(region.ch[, c("longitude", "latitude")])), ID = i)))
	overlay.1 <- given_table[!is.na(given_table %over% pol), ]
	if (dim(overlay.1)[1] > 0){
		overlay.1$region <- i
		overlayed <- rbind(overlayed, data.frame(overlay.1))
	}
}
write.table(overlayed, args[4], sep="\t", quote=F, row.names=F)
	
