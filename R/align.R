#'Gets DNAStringSet aligns it
#'@param mySequences DNAStringSet with all sequences in it
#'@param orfName ORF identifier name that will be used for file writing
#'@param outputDirectory path for files to be written in
#'@param ygeneSeq sequence of ORF of interest. If not given, last sequence of {mySequences} considered to be ygeneSeq
#'@param algorithm alignment algorithm to be used. Default is 'Muscle' but 'ClustalW' and 'ClustalOmega' are available
#'@return aligned DNAStringset, start and stop positions of orfName on DNAStringSet and AA alignment of the orfName sequence with other species
#'@export


align <- function(mySequences, orfName, outputDirectory=NULL,ygeneSeq=NULL,algorihm='Muscle') {


  dnaAlignmentList<-suppressWarnings(tryCatch({alignWoSmorf(mySequences,algorithm = 'Muscle',ygeneSeq)},
                                              error=function(cond){
                                                stop(paste('Alignment function gives error for',orfName, 'please make sure you are giving correct input'))
                                              }))

  DNAStr<-dnaAlignmentList[[1]]
  start<-dnaAlignmentList[[2]]
  stop<-dnaAlignmentList[[3]]
  if(start<0 || stop <0){
    stop(paste('No smorf sequence within these gene range for  ',orfName , sep = ''))
  }

  # subalignment####
  subalign<-subalignment(start,stop, DNAStr)[[1]]
  subalign<-tryCatch(DNAStringSet(muscle(subalign,quiet = T)), error=function(e) FALSE)
  if(is.logical(subalign)){
    stop(paste("One of the sequences does not have any nucleotide in the aligned region for ",orfName,sep = ''))
  }


  #Aminoacid translation####
  aa_alignment<-aaTranslation(subalign,DNAStr)
  if(is.logical(aa_alignment)){
    stop(paste("One of the sequences does not have any nucleotide in the aligned region for ",orfName,sep = ''))
  }

  if(is.null(outputDirectory)==F){
    writeXStringSet(subalign, file=paste(paste(outputDirectory,orfName, sep="/"),"subalignment.fa",sep = "_"))
    writeXStringSet(aa_alignment,file=paste(paste(outputDirectory,orfName, sep="/"),"AATranslation.fa",sep = "_"))

  }


  list(dnaAlignmentList=dnaAlignmentList,aa_alignment=aa_alignment)
}
