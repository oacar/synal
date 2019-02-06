
#'Find all start codons
#'@param aaSequence AAString instance
#'@return coordinates of all I if there is any. -1 if none exists
#'@export

findI<-function(aaSequence){
  coor<-c()
  for (i in 1:nchar(aaSequence)){
    if (as.character(aaSequence[i])=='I'){
      coor<-c(coor,i)
    }
  }
  if(is.null(coor)){
    return(-1)
  }else{
    return(toString(coor))
  }}
