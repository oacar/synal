#'Gets DNAStringSet aligns it
#'@param mySequences DNAStringSet with all sequences (including the ORF to be searched) in it
#'@param orfName ORF identifier name that will be used for file writing
#'@param path path for files to be written in
#'@return aligned DNAStringset, start and stop positions of orfName on DNAStringSet and AA alignment of the orfName sequence with other species
#'@export


align <- function(mySequences, orfName, path) {


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


  writeXStringSet(subalign, file=paste(paste(path,orfName, sep="/"),"subalignment.fa",sep = "_"))

  #Aminoacid translation####
  aa_alignment<-aaTranslation(subalign,DNAStr)
  if(is.logical(aa_alignment)){
    stop(paste("One of the sequences does not have any nucleotide in the aligned region for ",orfName,sep = ''))
  }

  writeXStringSet(aa_alignment,file=paste(paste(path,orfName, sep="/"),"AATranslation.fa",sep = "_"))

  list(dnaAlignmentList=dnaAlignmentList,aa_alignment=aa_alignment)
}
