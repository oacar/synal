#'Finds annotated ORF sequence from annotated sequence file downloaded from SGD
#'@param orfName annotated orf Y----- name
#'@param outputPath path if orfName sequence will be asked to written as separate file
#'@return sequence of given orf name as DNAString instance
#'@export


findYGeneSeq<- function(geneName,outputPath=NULL){

  #allnames <- scerAnnotatedSequences%>%names%>%str_split(',')%>%lapply('[',1)%>%str_split('[ ]+')%>%lapply('[',1)
  index <- which(scerAnnotatedSequences$orf_name%>%str_detect(geneName))


  if(index==0){

    warning(paste(geneName, 'cannot be found ', sep=" "))

  }else{
    Seq<-scerAnnotatedSequences[index,]$sequence%>%DNAStringSet()
    names(Seq) <- geneName

  if(is.null(outputPath)){
    return(Seq)
  }else{
    fileName<-paste0(outputPath,'/',geneName,'_sequence.fa')
    writeXStringSet(Seq, fileName)
  }

  return(Seq)
  }
}
