
#'Find all stop codons
#'@param aaSequence AAString instance
#'@return coordinates of all X if there is any. -1 if none exists
#'@export

findX<-function(aaSequence){
  coor<-c()
  for (i in 1:nchar(aaSequence)){
    if (as.character(aaSequence[i])=='X'){
      coor<-c(coor,i)
    }
  }
  if(is.null(coor)){
    return(-1)
  }else{
    return(toString(coor))
  }
}
