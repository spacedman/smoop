##' Make some random data
##'
##' Completely random data
##' @title Random data
##' @param n number of points
##' @return Spatial points data frame with p, Count, and N values
##' @export
##' @author Barry Rowlingson
makeRandom <- function(n){
  d=data.frame(x=runif(n),
    y=runif(n),
    Count=as.integer(runif(n)*20)
    )
  d$p=1+runif(n,0.2,0.9)
  d$N=as.integer(0.5+d$Count/d$p)
  coordinates(d)=~x+y
  d
}
##' Make some random data from a smooth
##'
##' Use RandomFields to create a smooth surface with given parameters,
##' then create some point data samples from it
##' 
##' @title Random data
##' @param n number of points
##' @param scale scale factor for GaussRF
##' @param alpha alpha for GaussRF
##' @param mean mean of GaussRF
##' @param nugget nugget for GaussRF
##' @param variance variance of GaussRF
##' @param maxF scale factor for field
##' @return list of pts (data points) and field (raster field object)
##' @export
##' @author Barry Rowlingson
makeSmooth <- function(n,scale=0.1,alpha=1,mean=0,nugget=0,variance=1,maxF=1){
  require(RandomFields)
  require(raster)
  x <- seq(0,1,len=50)
  y <- seq(0,1,len=50)
  f <- GaussRF(x=x, y=y, model="stable", grid=TRUE,param=c(mean, variance, nugget, scale, alpha))
  f <- maxF*(f-min(f))/(max(f)-min(f))
  f <- raster(list(x=x,y=y,z=f))
  d <- data.frame(x=runif(n),y=runif(n),
    N=as.integer(runif(n,100,150)))
  coordinates(d)=~x+y
  d$f = extract(f,coordinates(d))
  d$Count=as.integer(1 + d$N*d$f)
  list(pts=d,field=f)
}
