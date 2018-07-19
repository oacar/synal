
#'Start Codon alignment
#'@param alignment DNAStringSet or AAStringSet object of alignment
#'@return list of first 3 nucleotide if DNAStringSet or first aminoacid if AAStringSet
#'@export

startCodon <- function(alignment) {##alphabet check is to decide if it is AA or DNA seq
  alp<-alphabet(alignment,baseOnly=TRUE)#return first 3 or first letter
  
  sequence.count<-length(alignment)-1
  if(length(alp)==4){
    list.start<-list(DNAString(turnWoGaps(alignment[[1]]))[1:3])
    for(id in 2:sequence.count){
      list.start<-append(list.start,DNAString(turnWoGaps(alignment[[id]]))[1:3])
    }
    return(list.start)
  }else{
    list.start<-list(AAString(turnWoGaps(alignment[[1]]))[1])
    for(id in 2:sequence.count){
      list.start<-append(list.start,AAString(turnWoGaps(alignment[[id]]))[1])
    }
    return(list.start)
  }
}