#'Finds all overlapping ORFs
#'@param dna DNA sequence that will be searched for finding ORFs
#'@param range the range of interest. Only ORF ranges overlapping with this range will be considered
#'@return IRanges object with start and end positions determined using \code{dna} object
#'@import foreach
#'@import doParallel
#'@import IRanges
#'@export
#'

#library(foreach)
#library(doParallel)

findOverlappingOrfsPar <- function(dna,range) {
  startpositions <- str_locate_all(as.character(dna),'A[-]*T[-]*G')[[1]]
  if(length(startpositions)==0){
    return(NA)
  }
  stoppositions <- str_locate_all(as.character(dna),'T[-]*A[-]*A|T[-]*A[-]*G|T[-]*G[-]*A')[[1]]
  ranges <- IRanges()
  ranges <- foreach(i=1:length(startpositions[,1]),.packages = 'synal',.combine = 'irangesCombine') %:%
    foreach(j=1:nrow(stoppositions),.combine = 'irangesCombine') %dopar%{
      startpos <- startpositions[i,1]
      stoppos <- stoppositions[j,2]
      if(stoppos>startpos){
        seq <- turnWoGaps(subseq(dna,start = startpos,end=stoppos))
        aatr <- suppressWarnings(translate(DNAString((seq))))

        if(alphabetFrequency(aatr)['*']>1){
          return(IRanges(0,0))
        }else if(nchar((seq))%%3==0 & any(end(ranges)[start(ranges)==startpos]<stoppos)==F & any(start(ranges)[end(ranges)==stoppos]<startpos)==F){
          IRanges(startpos,stoppositions[j,2])
        }else(IRanges())
      }
    }
  #ranges <- ifelse(is.list(ranges),unlist(ranges)[[1]],ranges)
  if(is.list(ranges)){
    return(NULL)
  }
  ranges <- ranges[ranges%over%range]
  ranges <- ranges[width(overlapsRanges(ranges,range))>12]
  newranges <- IRanges()
  for(i in 1:length(unique(end(ranges)))){
    e <- unique(end(ranges))[i]
    eth <- ranges[end(ranges)==e]
    maxeth <- eth[which.max(width(eth))]
    newranges <- append(newranges,maxeth)
  }
  newranges
}

#'Combination function for IRanges to use in foreach
#'@param r1 iranges instance
#'@param r2 iranges instance
#'@return combined iranges object
#'@export
irangesCombine <- function(r1,r2){
  #d <- IRanges()
  d <- append(r1,r2)
#  d <- append(d,r2)
  d
}
