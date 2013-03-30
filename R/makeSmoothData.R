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

makeSmooth <- function(n,scale=0.1,alpha=1,mean=0,nugget=0,variance=1,maxF=1){
  require(RandomFields)
  x <- seq(0,1,len=50)
  y <- seq(0,1,len=50)
  f <- GaussRF(x=x, y=y, model="stable", grid=TRUE,param=c(mean, variance, nugget, scale, alpha))
  f <- maxF*(f-min(f))/(max(f)-min(f))
  f <- raster(list(x=x,y=y,z=f))
  d <- data.frame(x=runif(n),y=runif(n),
    N=100)
  coordinates(d)=~x+y
  d$f = extract(f,coordinates(d))
  d$Count=as.integer(1 + d$N*d$f)
  d
}
