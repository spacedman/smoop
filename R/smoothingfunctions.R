## (c) Copyright Lancaster University 2012,2013
##
# Peter's suggested kernel function
kernfun <- function(x){
    if(is.nan(x)){
        return(1)
    }
    else if(x<=1){
       return((1-x^2)^2)
    }
    else{
        return(0)
    }
}

##' SmoothIT function
##'
##' A function to fit the smoothed surface using the best M from cross validation
##'
##' @param pts1 matrix of grid cells
##' @param pts2 matrix of GP locations
##' @param pop population vector
##' @param Y quantity of interest
##' @param kernel Smoothing  kernel function. A function that takes a single argument, d, and returns a single numeric. Default is uniform kernel.
##' @param n optional maximum k for get.knnx. Probably best to set it for optimisation purposes. See ?get.nnx
##' @param lower lower bound for optimiser 
##' @param upper upper bound for optimiser
##' @param algorithm Method for computing k-nearest neighbour. See ?get.nnx
##' @return ...
SmoothIT <- function(pts1,pts2,pop,Y,kernel=NULL,n=NULL,lower,upper,algorithm="cover_tree"){
    Mopt <- getM(pts=pts2,pop=pop,Y=Y,kernel=kernel,n=n,algorithm=algorithm,lower=lower,upper=upper)$minimum
    ans <- evalsmooth(pts1=pts1,pts2=pts2,pop=pop,Y=Y,M=Mopt,n=n,kernel=kernel,algorithm=algorithm,check=TRUE)
    attr(ans,"Mopt") <- Mopt
    return(ans)
}

##' getM function
##'
##' A function to find the best M using cross validation
##'
##' @param pts matrix of GP locations
##' @param pop population vector
##' @param Y quantity of interest
##' @param kernel Smoothing  kernel function. A function that takes a single argument, d, and returns a single numeric. Default is uniform kernel.
##' @param n optional maximum k for get.knnx. Probably best to set it for optimisation purposes. See ?get.nnx
##' @param algorithm Method for computing k-nearest neighbour. See ?get.nnx
##' @param lower lower bound for optimiser 
##' @param upper upper bound for optimiser
##' @return ...

getM <- function(pts,pop,Y,kernel,n=NULL,algorithm,lower,upper){
  force(lower)
  force(upper)
  force(algorithm)
  f <- function(m){
    evalsmooth(pts1=pts,pts2=pts,pop=pop,Y=Y,M=m,n=n,kernel=kernel,algorithm=algorithm,opt=TRUE)
  }
  return(optimise(f,lower=lower,upper=upper,tol=100))
}


##' evalsmooth function
##'
##' A function to compute Peter's variable smoothing kernel.
##'
##' @param pts1 matrix of grid cells
##' @param pts2 matrix of GP locations
##' @param pop population vector
##' @param Y quantity of interest
##' @param M smoothing threshold, the number of people
##' @param n optional maximum k for get.knnx. See ?get.nnx
##' @param kernel Smoothing  kernel function. A function that takes a single argument, d, and returns a single numeric. Default is uniform kernel.
##' @param algorithm Method for computing k-nearest neighbour. See ?get.nnx
##' @param opt Setting opt to TRUE returns the cross validation variance, sum[(Y_j-s^(-j)(x_j))^2], where s^(-j) is the smoothing kernel without data Y_j from GP at x_j 
##' @return Matrix with three columns first column is smoothed values on locations pts1, second column is variance and third column is the zscore

evalsmooth <- function(pts1,pts2,pop,Y,M,n=NULL,kernel=NULL,algorithm="kd_tree",opt=FALSE,check=FALSE){
  require(futile.logger)
  if (M<0){
    stop("M must be non-negative.")
  }
  if (M>sum(pop)){
    stop("M must be smaller than the total population")
  }
  if(is.null(kernel)){
    kernel <- function(d){return(as.numeric(d<=1))}
  }
  if(is.null(n)){
    p <- sort(pop)
    cp <- cumsum(p)
    n <- min(which(cp>=M))
    flog.debug("Number of nearest neighbours = %d",n,name="smoop")
  }
    
  nn <- get.knnx(pts2,pts1,n,algorithm=algorithm)

  if (any(nn$nn.index==-1)){
    stop("Unable to compute at least one neighbour. Try increasing n.")
  }
    
  if(opt){
    q2 <- sapply(1:dim(pts2)[1],
                  function(i){Y[i]-evalsxwix(i,x=nn,n=n,pop=pop,Y=Y,M=M,kernel=kernel)$sx})
    
    retval <- sum(q2^2)
    if(is.na(retval)){
      stop("retval is NA")
    }
    return(retval) # cross validation variance
  }

  ans <- t(sapply(1:dim(pts1)[1],
                  function(i){evalsxvx(i,x=nn,n=n,pop=pop,Y=Y,M=M,kernel=kernel)}))

  rhohat <- sum(pop*Y)/sum(pop)

  sigma2 <- .findSigma2(Y,pop)

  ans[,2] <- sigma2*ans[,2]
  ans <- cbind(ans,(ans[,1]-rhohat)/sqrt(ans[,2]))
  colnames(ans) <- c("sx","vx","zscore")
  attr(ans,"sigma2") <- sigma2
  attr(ans,"M") <- M
  if(any(is.nan(ans[,1]))){
    warning("NaNs in smoothed grid. Try increasing M?")
  }
  return(ans)      
}  

##' evalsx function
##'
##' A function to compute the smoothed value of cell with index i. 
##'
##' @param i cell index
##' @param x output from get.knnx
##' @param n optional maximum k for get.knnx. See ?get.nnx
##' @param pop population vector
##' @param Y quantity of interest
##' @param M smoothing threshold, the number of people
##' @param kernel Smoothing  kernel function. A function that takes a single argument, d, and returns a single numeric.
##' @param opt whether or not in optimisation mode see optim parameter in function evalsmooth
##' @return ...
##' 

### this is the opt=FALSE case
evalsxvx = function(i,x,n,pop,Y,M,kernel){
  sxwix=evalsxwix(i,x,n,pop,Y,M,kernel)
  vx <- sum(sxwix$wix^2/sxwix$popi)/(sum(sxwix$wix)^2)
  return(c(sxwix$sx,vx))
}

### 
evalsxwix <- function(i,x,n,pop,Y,M,kernel){
    cs <- cumsum(pop[x$nn.index[i,]])
    cd <- cumsum(x$nn.dist[i,])
    idx <- min(which(cs>=M)) #min(which(cs>=M & cd>0))
    stopifnot(idx>=1)

    nni = x$nn.index[i,1:idx]
    nnd = x$nn.dist[i,1:idx]
    popi = pop[nni]
    
    wix <- sapply(nnd/x$nn.dist[i,idx],kernel)*popi
    if(length(wix)==1 & isTRUE(all.equal(wix,rep(0,length(wix))))){
      wix=kernel(0)*popi
    }
    
    wix[nnd==0] <- kernel(0)*popi[nnd==0]
    sx <- sum(wix*Y[nni])/sum(wix)

    return(list(sx=sx,wix=wix,popi=popi))

}



.findSigma2 <- function(Y,pop){
  or <- order(pop)
  opop <- pop[or]
  oY <- Y[or]
        
  prbs <- seq(0.1,1,length.out=min(c(100,floor(length(Y)/10))))    
  qt <- quantile(opop,prbs) # percentile boundaries
  ind <- 1 # opop[1] is th 0th quantile of the sample
  qtct <- 1
  for (i in 2:length(opop)){
    if(opop[i]<=qt[qtct]){
      ind <- c(ind,ind[i-1])
    }else{
      qtct <- qtct + 1
      ind <- c(ind,qtct)
    }
  }
  sk <- sapply(1:100,function(i){var(oY[ind==i])})
  Nk <- sapply(1:100,function(i){mean(opop[ind==i])})
  logsk <- log(sk)
  logNk <- log(Nk)
  
  coe <- coefficients(lm(logsk~logNk))
  sigma2 <- exp(coe[1])
  return(sigma2)
}
