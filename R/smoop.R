##
## (c) Copyright Lancaster University 2012,2013
##' ##' smoop - smooth over population centres
##'
##' Do a spatial smoothing where the baseline ratio is over a minimum population size.
##'
##' @param y formula for variable of interest
##' @param n formula for population
##' @param spdata spatial data frame
##' @param M minimal population size
##' @param bounds define bounds or clipping area
##' @param clip whether to clip to bounds
##' @param nx grid size x
##' @param ny grid size y
##' @return a raster stack with values and SEs and something
##' @export
smoop <- function(y,n,spdata,M,bounds=spdata,clip=FALSE,nx=64,ny=64,kernel=kernfun){
  require(raster)
  require(sp)
  require(FNN)
  yn=.getYN(y,n,spdata)
  box=bbox(bounds)
  xgrid=seq(box[1,1],box[1,2],len=nx)
  ygrid=seq(box[2,1],box[2,2],len=ny)

  gridxy = expand.grid(x=xgrid,y=ygrid)

  if(clip){
    ## FIXME check that bounds here is polygons
    xyo = overlay(SpatialPoints(gridxy,proj4string=proj4string(bounds)),bounds)
    gridxy=gridxy[!is.na(xyo),]
  }

  s = evalsmooth(gridxy,coordinates(spdata),yn$n,yn$y,M,kernel=kernel,opt=FALSE)

  s=data.frame(s)
  coordinates(s) <- gridxy
  gridded(s) <- TRUE
  s=brick(s)

  rhohat=sum(yn$n*yn$y)/sum(yn$n)
  attr(s,"rhohat")=rhohat
  return(s)
  
}
##' smoop cross-validation statistic
##'
##' 
##' @title smoop cross-validation statistic
##' @param y formula for the variable of interest
##' @param n formula for the population numbers
##' @param spdata spatial data frame
##' @param M minimum population count (vector)
##' @return vector of CV statistics corresponding to values of M
##' @author Barry Rowlingson
smoopCV <- function(y,n,spdata,M){
  yn = .getYN(y,n,spdata)
  pts=coordinates(spdata)

  laply(M,function(MM){  evalsmooth(pts,pts,yn$n,yn$y,MM,opt=TRUE)})
  
}

.getYN <- function(y,n,spdata){
  y=model.frame(y,spdata)
  if(ncol(y)!=1){stop("Incorrect model formula for y")}
  y=y[,1]
  n=model.frame(n,spdata)
  if(ncol(n)!=1){stop("Incorrect model formula for n")}
  n=n[,1]
  return(list(y=y,n=n))
}
