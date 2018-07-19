
#'Find all start codons
#'@param aaSequence AAString instance
#'@return coordinates of all M if there is any. -1 if none exists
#'@export

findM<-function(aaSequence){
  coor<-c()
  for (i in 1:nchar(aaSequence)){
    if (as.character(aaSequence[i])=='M'){
      coor<-c(coor,i)
    }
  }
  if(is.null(coor)){
    return(-1)
  }else{
    return(toString(coor))
  }}
