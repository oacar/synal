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
  #
  eg <- expand.grid(startpositions[,1],stoppositions[,2])
  eg_list <- eg[eg[,1]<eg[,2],] #take start positions < stop positions only and create iranges object with those positions
  eg_range <- IRanges(eg_list[,1],eg_list[,2])
  eg_ranges <- eg_range[eg_range%over%range]
  ranges <- IRanges()
  for(i in 1:length(eg_ranges)){
    startpos <- start(eg_ranges[i])
    stoppos <- end(eg_ranges[i])
    seq <- turnWoGaps(subseq(dna,start = startpos,end=stoppos))
    aatr <- suppressWarnings(translate(DNAString((seq))))

    if(alphabetFrequency(aatr)['*']>1){
      next
    }
    if(nchar((seq))%%3==0 & any(end(ranges)[start(ranges)==startpos]<stoppos)==F & any(start(ranges)[end(ranges)==stoppos]<startpos)==F){
      ranges <- append(ranges,eg_ranges[i])
    }
  }


  ranges <- ranges[ranges%over%range]
  ranges <- ranges[width(overlapsRanges(ranges,range))>12]
  newranges <- IRanges()
  #from the ORFs that end at the same position, take the one started earlier.
  for(i in 1:length(unique(end(ranges)))){
    e <- unique(end(ranges))[i]
    eth <- ranges[end(ranges)==e]
    maxeth <- eth[which.max(width(eth))]
    newranges <- append(newranges,maxeth)
  }
  newranges
}
