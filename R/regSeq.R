
#'Add gaps to sequence
#'Search sequence with any number of gaps in between using regular expressions
#'adds '[-]*' in between every letter. This helps to seacrh orf sequence on the aligned sequence since there is no way to know where the gaps are
#'@param seq character string
#'@return character string
#'@export

regSeq<-function(seq){
  regSeq<-''
  seq<-as.character(seq)
  for (i in 1:(nchar(seq)-1)){
    regSeq<-paste(regSeq, substr(seq,i,i), '[-]*', sep = '')
  }
  regSeq<-paste(regSeq, substr(seq,nchar(seq),nchar(seq)), sep = '')
  return(regSeq)
}
