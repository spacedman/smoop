##' smoop - smooth over population centres
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
  y=model.frame(y,spdata)
  if(ncol(y)!=1){stop("Incorrect model formula for y")}
  y=y[,1]
  n=model.frame(n,spdata)
  if(ncol(n)!=1){stop("Incorrect model formula for n")}
  n=n[,1]
  box=bbox(bounds)
  xgrid=seq(box[1,1],box[1,2],len=nx)
  ygrid=seq(box[2,1],box[2,2],len=ny)

  gridxy = expand.grid(x=xgrid,y=ygrid)

  if(clip){
    ## FIXME check that bounds here is polygons
    xyo = overlay(SpatialPoints(gridxy,proj4string=proj4string(bounds)),bounds)
    gridxy=gridxy[!is.na(xyo),]
  }

  s = evalsmooth(gridxy,coordinates(spdata),n,y,M,kernel=kernel)

  s=data.frame(s)
  coordinates(s) <- gridxy
  gridded(s) <- TRUE
  s=brick(s)

  rhohat=sum(n*y)/sum(n)
  attr(s,"rhohat")=rhohat
  return(s)
  
}

