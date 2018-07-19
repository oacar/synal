#'Check if all codons are same
#'@param alignedLetter list of AAString or DNAString instances, possibly result of startCodon or stopCodon functions
#'@return TRUE if all same , FALSE otherwise

#'@export

checkAlignment<-function(alignedLetter){#Aligned letter should be list of AAString or DNAString instances
  check=c()
  for(i in 2:length(alignedLetter)){
    check<-append(check,alignedLetter[[1]]==alignedLetter[[i]])
  }
  #return(alignedLetter[[1]]==alignedLetter[[2]] && alignedLetter[[1]]==alignedLetter[[3]] && alignedLetter[[1]]==alignedLetter[[4]] && alignedLetter[[1]]==alignedLetter[[4]])
  return(all(check))
}
