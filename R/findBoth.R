#'Search for both start and stop codons to get the overlapping orf
#'@param DNAStr aligned DNAStringSet object
#'@param speciesIndex the index for which sequence in the DNAStringSet object to be searched for ATG and stop codons
#'@param start the start position of the S. Cerevisiae ORF over the alignment object
#'@param orfName the corresponding orf name. Used for a clear warning message
#'@param stop stop position of the S. Cerevisiae ORF over the alignment object. It is used for when the Start codon not found
#'
#'@return The list of boolean value for which will be used to determine if another subalignment and AA translation will be written and new start and stop positions
#'      
#'
#'@export


findBoth <- function(DNAStr, speciesIndex, start, stop,orfName) {
  NextStop<-findCodon(DNAStr[[speciesIndex]][start:length(DNAStr[[speciesIndex]])],c('TAA', 'TAG', 'TGA'))+start-1
  PrevStart<-findCodon(DNAStr[[speciesIndex]][1:start-1], 'ATG')
  writeFile<-FALSE
  if(PrevStart<=0 || PrevStart==Inf || NextStop<=0 || NextStop==Inf){
    return(list(writeFile=NULL, 
                PrevStart=PrevStart,
                NextStop=NextStop))
  }
  dummyStopCheck<-findCodon(DNAStr[[speciesIndex]][PrevStart:length(DNAStr[[speciesIndex]])],c('TAA', 'TAG', 'TGA'))+PrevStart-1
  if(dummyStopCheck<start){
    warning(paste('No start codon found! for ', orfName, sep = ''))
    PrevStart<-start
    NextStop<-stop
  }else{
    writeFile=TRUE
  }
  return(list(writeFile=writeFile, 
              PrevStart=PrevStart,
              NextStop=NextStop))
}
