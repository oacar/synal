#'Gets DNAStringSet aligns it
#'@param mySequences DNAStringSet with all sequences (including the ORF to be searched) in it
#'@param orfName ORF identifier name that will be used for file writing
#'@param path path for files to be written in
#'@return aligned DNAStringset, start and stop positions of orfName on DNAStringSet and AA alignment of the orfName sequence with other species
#'@export


align <- function(mySequences, orfName, path) {

  dnaAlignmentList<-suppressWarnings(tryCatch({alignWoSmorf(mySequences,algorithm = 'Muscle')},
                                              error=function(cond){
                                                alignWoSmorf(mySequences,algorithm = 'ClustalW')
                                              }))

  DNAStr<-dnaAlignmentList[[1]]
  start<-dnaAlignmentList[[2]]
  stop<-dnaAlignmentList[[3]]
  if(start<0 || stop <0){
    stop(paste('No smorf sequence within these gene range for  ',orfName , sep = ''))
    #next
    # if(TRUE){
    # nPath<-paste(noseqrange,orfName, sep="/") #If regex cannot find the smorf sequence #1-Seq might be changed or #2-input files are wrong
    #file.rename(path,nPath)
    #if(TRUE){next}
    #}
  }
  #
  if(start<=20){
    sStart<-1
  }else{
    sStart<-start-20
  }
  if(abs(stop-length(DNAStr[[1]]))<=20){
    sStop<-length(DNAStr[[1]])
  }else{
    sStop<-stop+20
  }
  smallerDNAStr<-subalignment(sStart,sStop, DNAStr)[[1]]
  smallerDNAStr<-append(smallerDNAStr,mySequences[length(mySequences)])

  sDnaAlignmentList<-tryCatch({alignWoSmorf(smallerDNAStr,algorithm = 'Muscle')},
                              error=function(cond){
                                alignWoSmorf(smallerDNAStr,algorithm = 'ClustalW')
                              })

  sDNAStr<-sDnaAlignmentList[[1]]
  sStart<-sDnaAlignmentList[[2]]
  sStop<-sDnaAlignmentList[[3]]

  # subalignment####
  subalign<-subalignment(sStart,sStop, sDNAStr)[[1]]
  subalign<-append(subalign, mySequences[length(mySequences)])
  #subalign<-DNAStringSet(msa(subalign, 'Muscle', order = 'input'))
  #subalign<-chartr('N','-',subalign)warni
  subalign<-tryCatch(DNAStringSet(muscle(subalign,quiet = T)), error=function(e) FALSE)
  if(is.logical(subalign)){
    stop(paste("One of the sequences does not have any nucleotide in the aligned region for ",orfName,sep = ''))

    #newPath<-paste(noregion,orfName, sep="/")
    # if(TRUE){next}
  }

  # writeXStringSet(DNAStr, file=paste(paste(path,orfName, sep="/"),"alignment.fa",sep = "_"))

  writeXStringSet(subalign, file=paste(paste(path,orfName, sep="/"),"subalignment.fa",sep = "_"))

  #Aminoacid translation####
  aa_alignment<-aaTranslation(subalign,DNAStr)
  if(is.logical(aa_alignment)){
    stop(paste("One of the sequences does not have any nucleotide in the aligned region for ",orfName,sep = ''))
    #next
    #newPath<-paste(noregion,orfName, sep="/")
    #if(TRUE){next}
  }
  names(aa_alignment)[length(aa_alignment)] <- orfName
  writeXStringSet(aa_alignment,file=paste(paste(path,orfName, sep="/"),"AATranslation.fa",sep = "_"))

  list(dnaAlignmentList=dnaAlignmentList,aa_alignment=aa_alignment)
}
