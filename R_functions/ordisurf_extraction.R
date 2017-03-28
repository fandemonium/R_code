## nmds: nmds object
## env:  2 column data frame; 1 column has header same as SAMPLE_ID, 1 column is the env variable
## SAMPLE_ID: one of the env column name for matching
ordi.sf<-function(nmds, env, SAMPLE_ID){
	library(vegan)
        NMDS<-data.frame(nmds$points)[, 1:2]
        NMDS$SAMPLES<-row.names(NMDS)
        NMDS.si<-merge(NMDS, env, by.x="SAMPLES", by.y=SAMPLE_ID)
	env.sf<-ordisurf(NMDS.si[, c("MDS1", "MDS2")]~ env[, !names(env) %in% SAMPLE_ID], plot=F)
	print(summary(env.sf))

extract.xyz <- function(obj) {
    xy <- expand.grid(x = obj$grid$x, y = obj$grid$y)
    xyz <- cbind(xy, c(obj$grid$z))
    names(xyz) <- c("x", "y", "z")
    xyz <- data.frame(na.omit(xyz))
    return(xyz)
}

	contour.vals <- extract.xyz(obj=env.sf)
}
