
#'Stop codon alignment
#'@param alignment DNAStringSet or AAStringSet object of alignment
#'@return list of last 3 nucleotide if DNAStringSet or last aminoacid if AAStringSet
#'@export

stopCodon <- function(alignment) {##alphabet check is to decide if it is AA or DNA seq
  alp<-alphabet(alignment,baseOnly=TRUE)#return last 3 or last 1 letter
  l<-length(alignment[[1]])
  sequence.count<-length(alignment)-1

  if(length(alp)==4){
    a1<-DNAString(turnWoGaps(alignment[[1]]))
    l1<-length(a1)
    list.stop<-list(a1[(l1-2):l1])
    for(id in 2:sequence.count){
      seq<-DNAString(turnWoGaps(alignment[[id]]))
      seq.length<-length(seq)

      list.stop<-append(list.stop,seq[(seq.length-2):seq.length])
    }
    return(list.stop)
  }else{
    a1<-AAString(turnWoGaps(alignment[[1]]))
    l1<-length(a1)

    list.stop<-list(a1[l1])
    for(id in 2:sequence.count){
      seq<-AAString(turnWoGaps(alignment[[id]]))
      seq.length<-length(seq)

      list.stop<-append(list.stop,seq[seq.length])
    }
    return(list.stop)
  }
}
