#'
#'@param DNAStr DNA alignment of syntenic block. This is the full alignment not subalignment
#'@param j this is the iterator for the ID of species of interest on DNAStr and types vector
#'@param r this is the range parameter. This specifies how far from the proto-gene region will be taken into account for ORF finding default =200
#'@param start start position of ATG of proto-gene on DNAStr
#'@param stop last nucleotide position of proto-gene on DNAStr
#'@param ygeneSeq proto-gene sequence without gaps. This is used for determining the frame of proto-gene with respect to the aligned other species ORF
#'@param types this is used for names. by default c('Spar','Smik','Skud','Sbay')
#'
#'@return list of three elements
#'@return dna = dna alignment of jth species orf with proto-gene, aa = aa alignment with full jth ORF, aaOverlap = overlapping regions of both ORFs
#'
#'@export

findBestOverlap <- function(DNAStr, j, r,start, stop, ygeneSeq, types , u=NULL) {
  u <- ifelse(is.null(u),r,u)
  subseq <- subseq(DNAStr[[j]],start=start-r,end=stop+u)
  ranges <- findOverlappingOrfs(subseq,range = IRanges(r,stop-start+r))
  if(length(ranges)==0) return(FALSE)
  score=0
  best <- NULL
  for( u in 1:length(ranges)){
    startPosRanges <- start(ranges)[u]
    stopPosRanges <- end(ranges)[u]
    sequenceRanges <-subseq(subseq,startPosRanges,stopPosRanges)
    # names(sequenceRanges) <- types[[j]]
    if(startPosRanges<r+1){
      smorfSubSeq <- subseq(DNAStr[1],start-r-1+startPosRanges,start-r-1+stopPosRanges)
      place <- str_locate(turnWoGaps(smorfSubSeq[[1]]),as.character(subseq(ygeneSeq,1,7))[[1]])
      longer=T
      if(place[1,1]%%3==2){
        smorfSubSeq <- subseq(DNAStr[1],start-r-2+startPosRanges,start-r-1+stopPosRanges)

      }else if(place[1,1]%%3==0){
        smorfSubSeq <- subseq(DNAStr[1],start-r-3+startPosRanges,start-r-1+stopPosRanges)
      }
    }else{
      smorfSubSeq <- subseq(DNAStr[1],start-r-1+startPosRanges,start-r-1+stopPosRanges)
      place <- str_locate(as.character(ygeneSeq),turnWoGaps(subseq(smorfSubSeq[[1]],1,ifelse(width(smorfSubSeq)<10,width(smorfSubSeq),10))))
      longer=F
      if(place[1,1]%%3==2){
        smorfSubSeq <- subseq(DNAStr[1],start-r-2+startPosRanges,start-r-1+stopPosRanges)

      }else if(place[1,1]%%3==0){
        smorfSubSeq <- subseq(DNAStr[1],start-r-3+startPosRanges,start-r-1+stopPosRanges)
      }
    }
    dnaSet <- append(DNAStringSet(sequenceRanges),smorfSubSeq)
    names(dnaSet)[1] <- types[j]
    aaRanges <- aaTranslation(dnaSet,DNAStr[c(j,1)])
    aaRangesSmall <- findSmallFrame(aaRanges,2,aaRanges[[2]],overlap=T,longer=longer)[c(2,1)]
    newscore <- calcIdentity(aaRangesSmall)*width(aaRanges[1])
    if(score<newscore[2,1]){
      score=newscore[2,1]
      bestAA <- aaRanges[c(2,1)]
      bestAAsmall <- aaRangesSmall
      bestDna <- dnaSet[c(2,1)]
    }
  }
  return(list(dna=bestDna,aa=bestAA, aaOverlap=bestAAsmall))
}
