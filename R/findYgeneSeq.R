#'Finds annotated ORF sequence from annotated sequence file downloaded from SGD
#'@param geneName annotated orf Y----- name
#'@param outputPath path if orfName sequence will be asked to written as separate file
#'@return sequence of given orf name as DNAString instance
#'@importFrom utils data
#'@export


findYGeneSeq<- function(geneName,outputPath=NULL){
  data(scerAnnotatedSequences,envir = environment())
  #allnames <- scerAnnotatedSequences%>%names%>%str_split(',')%>%lapply('[',1)%>%str_split('[ ]+')%>%lapply('[',1)
  seq <- tryCatch({synal::scerAnnotatedSequences[[geneName]]}, error = function(e){NULL})


  if(is.null(seq)){

    warning(paste(geneName, 'cannot be found ', sep=" "))

  }else{
    seq<-seq%>%DNAStringSet()
    names(seq) <- geneName

  if(is.null(outputPath)){
    return(seq)
  }else{
    fileName<-paste0(outputPath,'/',geneName,'_sequence.fa')
    writeXStringSet(seq, fileName)
  }

  return(seq)
  }
}
