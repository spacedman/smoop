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

smoop1d <- function(y, n, x, M, bounds=range(x, na.rm=TRUE), nx=64, kernel=kernfun){
 xgrid=seq(min(bounds,na.rm=TRUE),max(bounds,na.rm=TRUE),len=nx)
 ygrid = 0
 gridxy = expand.grid(x=xgrid,y=ygrid)
 s = evalsmooth(gridxy,cbind(x,0),n,y,M,kernel=kernel,opt=FALSE)
 s = data.frame(s)
 s$x = xgrid
 s
}

smoopCVbiOneD <- function(N, y, n, x, M, .progress=smoopProgress()){
  spd = data.frame(y=y,n=n, xc=x)
  spd$yc = 0
  coordinates(spd) <- ~xc+yc
  smoopCVbiN(N, ~y, ~n, spd, M, .progress=.progress)
}

##' ##' smoop cross-validation statistic
##'
##' 
##' @title smoop cross-validation statistic
##' @param y formula for the variable of interest
##' @param n formula for the population numbers
##' @param spdata spatial data frame
##' @param M minimum population count (vector)
##' @return vector of CV statistics corresponding to values of M
##' @export
##' @author Barry Rowlingson
smoopCV <- function(y,n,spdata,M){
  yn = .getYN(y,n,spdata)
  pts=coordinates(spdata)

  laply(M,function(MM){  evalsmooth(pts,pts,yn$n,yn$y,MM,opt=TRUE)})
  
}
##' Cross validation by bisection
##'
##' Divide data into two parts, fit with one, compute error with the other
##' and do the same the other way round
##' @title Cross Validation
##' @param y formula for the variable of interest
##' @param n formula for the population numbers
##' @param spdata spatial data frame 
##' @param M vector of smoothing constants
##' @param .progress progress indicator
##' @return vector of sum of squared differences
##' @export
##' @author Barry Rowlingson

smoopCVbi <- function(y,n,spdata,M,.progress=smoopProgress()){
  require(plyr)
  yn = .getYN(y,n,spdata)
  pts = coordinates(spdata)

  np = nrow(spdata)
  n1 = as.integer(np/2)
  p1 = sort(sample(np,n1))
  p2 = (1:np)[-p1]

  meansqd = laply(M,
    function(MM){
      s1 = yn$y[p1]-evalsmooth(pts[p1,],pts[p2,],yn$n[p2],yn$y[p2],MM)[,1]
      s2 = yn$y[p2]-evalsmooth(pts[p2,],pts[p1,],yn$n[p1],yn$y[p1],MM)[,1]
      msqd=mean(c(s1^2,s2^2))
      msqd
    },
    .progress=.progress
    )
  meansqd
  
}

smoopCVbiN <- function(N,y,n,spdata,M,.progress=smoopProgress()){
  m = laply(1:N,function(i)smoopCVbi(y,n,spdata,M,.progress="none"),.progress=.progress)
  s = apply(m,1,mean)
  laply(1:N,function(i){m[i,]-s[i]})
}

smoopLooS <- function(y,n,spdata,M,nout=nrow(spdata),j,.progress=smoopProgress()){
  require(plyr)
  laply(M,
        function(MM){
          smoopLoo(y,n,spdata,MM,nout=nrow(spdata),j,.progress="none")$mssq
        },
        .progress=.progress
      )
}

smoopLoo <- function(y,n,spdata,M,nout=nrow(spdata),j,.progress=.smoopProgress()){
  yn = .getYN(y,n,spdata)
  pts = coordinates(spdata)
  #sx = evalsmooth(pts,pts,yn$n, yn$y, M, opt=TRUE)
  if(missing(j)){
    j = sort(sample(nrow(spdata),nout))
  }
  sxj = laply(j,
    function(jj){
      #evalsxwix(i,x=nn,n=n,pop=pop,Y=Y,M=M,kernel=kernel)$sx
      evalsmooth(pts[jj,,drop=FALSE],pts[-jj,],yn$n[-jj],yn$y[-jj],M)
    },
    .progress=.progress

    )
  Ysx=data.frame(j=j,y=yn$y[j],sx=sxj[,1])
  mssq=mean((Ysx$y-Ysx$sx)^2)
  list(Ysx=Ysx,mssq=mssq)
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
