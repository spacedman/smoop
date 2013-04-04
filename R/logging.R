##' logging level for smoop
##'
##' 
##' @title Logging level
##' @return nothing
##' @author Barry Rowlingson
##' @export
smoopLogLevel <- function(){
  flog.logger(name="smoop")$threshold
}
##' Test if log level for smoop is INFO or above
##'
##' 
##' @title Logging level
##' @return TRUE is smoop logging level >= INFO
##' @export
##' @author Barry Rowlingson
smoopLogInfo <- function(){
  smoopLogLevel()>=INFO
}
##' Smoop Progress indicator
##'
##' Decide on a progress indicator based on the log level
##' @title Logging level
##' @return string describing a progress bar
##' @export
##' @author Barry Rowlingson
smoopProgress <- function(){
  if(smoopLogInfo()){
    return("text")
  }else{
    return("none")
  }
}
