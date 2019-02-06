
#'Search for stop codon to get the overlapping orf
#'@param DNAStr aligned DNAStringSet object
#'@param speciesIndex the index for which sequence in the DNAStringSet object to be searched for ATG and stop codons
#'@param start the start position of the S. Cerevisiae ORF over the alignment object
#'@param orfName the corresponding orf name. Used for a clear warning message
#'@param stop stop position of the S. Cerevisiae ORF over the alignment object. It is used for when the Start codon not found
#'
#'@return The list of boolean value for which will be used to determine if another subalignment and AA translation will be written and new start and stop positions
#'
#'
#'
#'@export

findNextStop <- function(DNAStr, speciesIndex, start, stop,orfName, l=-1) {
  NextStop<-findCodon(DNAStr[[speciesIndex]][start:length(DNAStr[[speciesIndex]])],c('TAA', 'TAG', 'TGA'),l)+start-1
  writeFile<-FALSE
  if(NextStop<=0 || NextStop==Inf){
    return(list(writeFile=NULL,
                PrevStart=start,
                NextStop=stop))
  }
  writeFile=TRUE
  return(list(writeFile=writeFile,
              PrevStart=start,
              NextStop=NextStop))
}
