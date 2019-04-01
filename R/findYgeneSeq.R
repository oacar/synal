#'Finds annotated ORF sequence from annotated sequence file downloaded from SGD
#'@param orfName annotated orf Y----- name
#'@return sequence of given orf name as DNAString instance
#'@export


findYGeneSeq<- function(geneName){
  for (i in 1:length(scerAnnotatedSequences)){
    if(geneName ==  strsplit(strsplit(names(scerAnnotatedSequences)[i], ',')[[1]][1], ' ')[[1]][1]){
      index<-i
    }
  }
  if(index==0){

    warning(sprintf("%s", paste(geneName, 'cannot be found ', sep=" ")))

  }else{  Seq<-scerAnnotatedSequences[index]

  if(is.null(outputPath)){
    return(Seq)
  }else{
    fileName<-sprintf("%s", paste(outputPath,geneName, sep="/"))
    fileName<-sprintf("%s", paste(fileName,"sequence", sep="_"))
    fileName<-sprintf("%s", paste(fileName,'.fa', sep=""))
    writeXStringSet(Seq, fileName)
  }

  return(Seq)
  }
}
