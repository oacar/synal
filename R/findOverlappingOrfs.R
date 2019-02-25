#'Finds all overlapping ORFs
#'@param dna DNA sequence that will be searched for finding ORFs
#'@param range the range of interest. Only ORF ranges overlapping with this range will be considered
#'@return IRanges object with start and end positions determined using \code{dna} object
#'@export
findOverlappingOrfs <- function(dna,range) {
  startpositions <- str_locate_all(as.character(dna),'A[-]*T[-]*G')[[1]]
  if(length(startpositions)==0){
    return(NA)
  }
  stoppositions <- str_locate_all(as.character(dna),'T[-]*A[-]*A|T[-]*A[-]*G|T[-]*G[-]*A')[[1]]
  ranges <- IRanges()
  for(i in startpositions[,1]){
    startpos <- i
    for(j in 1:nrow(stoppositions)){
      stoppos <- stoppositions[j,2]
      if(stoppos>startpos){
        seq <- turnWoGaps(subseq(dna,start = startpos,end=stoppos))
        aatr <- suppressWarnings(translate(DNAString((seq))))

        if(alphabetFrequency(aatr)['*']>1){
          next
        }
        if(nchar((seq))%%3==0 & any(end(ranges)[start(ranges)==startpos]<stoppos)==F & any(start(ranges)[end(ranges)==stoppos]<startpos)==F){
          ranges <- append(ranges,IRanges(startpos,stoppositions[j,2]))
        }
      }
    }
  }
  ranges <- ranges[ranges%over%range]
  ranges <- ranges[width(overlapsRanges(ranges,range))>12]
  ranges
}
