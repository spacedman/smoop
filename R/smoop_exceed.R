#' 
#' compute exceedence probs from the output of smoop
#'
#' @param sv the output from smoop
#' @param c the critical value for exceedence
#' @param rhohat the rho-hat value from smoop (usually in the smoop object)
#' 
smoop_exceed <- function(sv,c,rhohat=attr(sv,"rhohat")){
  require(raster)
  require(classInt)
   upp = raster(sv[[1]])
   low = raster(sv[[1]])
   quints = raster(sv[[1]])
   
   values(upp) <- pnorm((values(sv[[1]])-c*rhohat)/sqrt(values(sv[[2]])),lower.tail=FALSE)
   values(low) <- pnorm((values(sv[[1]])-(1/c)*rhohat)/sqrt(values(sv[[2]])))

   class = raster(sv[[1]])
   values(class)=0.5 #((1-values(upp))+values(low))/2
#   gridcol[]=(upp+low)/2
   class[upp<0.05] <- (1-upp)[upp<0.05]
   class[low<0.05] <-  low[low<0.05]

   class[is.na(sv[[1]])]=NA

   v = values(sv[[1]])
   ci = findCols(classIntervals(v,n=5))
   quints[]=ci
   
   s = stack(class,upp,low,quints)
   names(s)=c("class","Upper","Lower","Quintiles")
   s
   
 }
