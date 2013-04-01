smoopLogLevel <- function(){
  flog.logger(name="smoop")$threshold
}

smoopLogInfo <- function(){
  smoopLogLevel()>=INFO
}

smoopProgress <- function(){
  if(smoopLogInfo()){
    return("text")
  }else{
    return("none")
  }
}
