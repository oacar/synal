<<<<<<< HEAD
#'Finds best overlapping ORF in terms of the number of identical aminoacids
=======
#'
>>>>>>> c1655aa13be547d7cec4b70c6fa5291cc3ca1da2
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

<<<<<<< HEAD
findBestOverlap <- function(DNAStr, j, r,start, stop, ygeneSeq, types , k=NULL) {
  k <- ifelse(is.null(k),r,k)
  subseq <- subseq(DNAStr[[j]],start=ifelse(start-r<1,1,start-r),end=ifelse(stop+k<length(DNAStr[[j]]),stop+k,length(DNAStr[[j]])))
  ranges <- findOverlappingOrfsPar(subseq,range = IRanges(r,stop-start+r))
  itr <- 0
  if(length(ranges)==0){
    return(NULL)
  }else if(all(is.na(ranges)) & itr<3 ){
    r=r-200
    subseq <- subseq(DNAStr[[j]],start=ifelse(start-r<1,1,start-r),end=ifelse(stop+k<length(DNAStr[[j]]),stop+k,length(DNAStr[[j]])))

    ranges <- findOverlappingOrfsPar(subseq,range = IRanges(r,stop-start+r))
    itr <- itr+1
  }
  score=-1
  best <- NULL
  bestAA <- NULL
  bestDna <- NULL
  bestAAsmall <- NULL
=======
findBestOverlap <- function(DNAStr, j, r,start, stop, ygeneSeq, types , u=NULL) {
  u <- ifelse(is.null(u),r,u)
  subseq <- subseq(DNAStr[[j]],start=start-r,end=stop+u)
  ranges <- findOverlappingOrfs(subseq,range = IRanges(r,stop-start+r))
  if(length(ranges)==0) return(FALSE)
  score=0
  best <- NULL
>>>>>>> c1655aa13be547d7cec4b70c6fa5291cc3ca1da2
  for( u in 1:length(ranges)){
    startPosRanges <- start(ranges)[u]
    stopPosRanges <- end(ranges)[u]
    sequenceRanges <-subseq(subseq,startPosRanges,stopPosRanges)
    # names(sequenceRanges) <- types[[j]]
    if(startPosRanges<r+1){
      smorfSubSeq <- subseq(DNAStr[1],start-r-1+startPosRanges,start-r-1+stopPosRanges)
<<<<<<< HEAD
      place <- str_locate(turnWoGaps(smorfSubSeq[[1]]),as.character(subseq(ygeneSeq,1,10))[[1]])
      if(any(is.na(place))){
        next
      }
      longer=T
      if(place[1,1]%%3==0){
        nogapscer <- turnWoGaps(DNAStr[[1]])
        placenogap <- str_locate(nogapscer,as.character(subseq(ygeneSeq,1,10))[[1]])
        smorfSubSeq <- xscat(str_sub(nogapscer,placenogap[1,1]-1,placenogap[1,1]-1),smorfSubSeq)

      }else if(place[1,1]%%3==2){
        nogapscer <- turnWoGaps(DNAStr[[1]])
        placenogap <- str_locate(nogapscer,as.character(subseq(ygeneSeq,1,10))[[1]])
        smorfSubSeq <- xscat(str_sub(nogapscer,placenogap[1,1]-2,placenogap[1,1]-1),smorfSubSeq)
      }
    }else{
      smorfSubSeq <- subseq(DNAStr[1],start-r-1+startPosRanges,start-r-1+stopPosRanges)
      nogap <- DNAString(turnWoGaps(smorfSubSeq[[1]]))
      place <- str_locate(as.character(ygeneSeq),as.character(subseq(nogap,1,ifelse(length(nogap)<10,length(nogap),10))))
      if(any(is.na(place))){
        next
      }

      longer=F
      if(place[1,1]%%3==2){
        smorfSubSeq <- xscat(str_sub(as.character(ygeneSeq),place[1,1]-1,place[1,1]-1),smorfSubSeq)
#subseq(DNAStr[1],start-r-2+startPosRanges,start-r-1+stopPosRanges)

      }else if(place[1,1]%%3==0){
        smorfSubSeq <- xscat(str_sub(as.character(ygeneSeq),place[1,1]-2,place[1,1]-1),smorfSubSeq)
      }
    }
    dnaSet <- append(DNAStringSet(sequenceRanges),smorfSubSeq)
    names(dnaSet) <- names(DNAStr)[c(j,1)]
    aaRanges <- aaTranslation(dnaSet,DNAStr[c(j,1)])
    aaRangesSmall <- findSmallFrame(aaRanges,translate(subseq(ygeneSeq,1,12))[[1]],aaRanges[[2]],overlap=T,longer=longer)[c(2,1)]
    if(is.logical(aaRangesSmall)){
      warning(paste(names(ygeneSeq),types[j],'has a problematic overlapping ORF. Nice to check'))
      next
    }
=======
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
>>>>>>> c1655aa13be547d7cec4b70c6fa5291cc3ca1da2
    newscore <- calcIdentity(aaRangesSmall)*width(aaRanges[1])
    if(score<newscore[2,1]){
      score=newscore[2,1]
      bestAA <- aaRanges[c(2,1)]
      bestAAsmall <- aaRangesSmall
      bestDna <- dnaSet[c(2,1)]
<<<<<<< HEAD
    }else if(score==newscore[2,1]){
      if(length(bestAA[[2]])<length(aaRanges[[1]])){
        score=newscore[2,1]
        bestAA <- aaRanges[c(2,1)]
        bestAAsmall <- aaRangesSmall
        bestDna <- dnaSet[c(2,1)]
      }
    }
  }
  if(is.null(bestDna) | is.null(bestAA) | is.null(bestAAsmall)){
    return(NULL)
  }else{
    return(list(dna=bestDna,aa=bestAA, aaOverlap=bestAAsmall))
  }

}
# rm(bestDna)
# rm(bestAA)
# rm(bestAAsmall)
=======
    }
  }
  return(list(dna=bestDna,aa=bestAA, aaOverlap=bestAAsmall))
}
>>>>>>> c1655aa13be547d7cec4b70c6fa5291cc3ca1da2
