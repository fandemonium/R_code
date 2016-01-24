ordi.sf<-function(nmds, env){

env.sf<-ordisurf(nmds ~ env, plot=F, scaling=3)

extract.xyz <- function(obj) {
    xy <- expand.grid(x = obj$grid$x, y = obj$grid$y)
    xyz <- cbind(xy, c(obj$grid$z))
    names(xyz) <- c("x", "y", "z")
    xyz <- data.frame(na.omit(xyz))
    return(xyz)
}

contour.vals <- extract.xyz(obj=env.sf)
}
