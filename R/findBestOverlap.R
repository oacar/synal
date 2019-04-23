#'Finds best overlapping ORF in terms of the number of identical aminoacids
#'@param DNAStr DNA alignment of syntenic block. This is the full alignment not subalignment
#'@param j this is the iterator for the ID of species of interest on DNAStr and types vector
#'@param start start position of ATG of proto-gene on DNAStr
#'@param stop last nucleotide position of proto-gene on DNAStr
#'@param ygeneSeq proto-gene sequence without gaps. This is used for determining the frame of proto-gene with respect to the aligned other species ORF
#'@param types this is used for names. by default c('Spar','Smik','Skud','Sbay')
#'
#'@return list of three elements
#'@return dna = dna alignment of jth species orf with proto-gene, aa = aa alignment with full jth ORF, aaOverlap = overlapping regions of both ORFs
#'@import IRanges
#'@export

findBestOverlap <- function(DNAStr, j,start, stop, ygeneSeq, types , k=NULL) {
  subseq <- DNAStr[[j]]%>%as.character()
  #tic('ranges')
  map=map_alignment_sequence(subseq,turnWoGaps(subseq))

  range=IRanges(map[map==start]%>%names()%>%as.integer(),map[map==stop]%>%names()%>%as.integer())
  ranges <- findOverlappingOrfs(subseq%>%turnWoGaps(),range = range)#IRanges(start,stop))
  ranges <- IRanges(map[start(ranges)],map[end(ranges)])
  if(length(ranges)==0){
    return(NULL)
  }
 # toc()
  itr <- 0
 # tic('best')
  score=-1
  best <- NULL
  bestAA <- NULL
  bestDna <- NULL
  bestAAsmall <- NULL
  
  #print(j)
  for( u in 1:length(ranges)){
    startPosRanges <- start(ranges)[u]
    stopPosRanges <- end(ranges)[u]
    sequenceRanges <-subseq(subseq,startPosRanges,stopPosRanges)
    # names(sequenceRanges) <- types[[j]]
    if(startPosRanges<start+1){
      smorfSubSeq <- subseq(DNAStr[1],startPosRanges,stopPosRanges)
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
      smorfSubSeq <- subseq(DNAStr[1],startPosRanges,stopPosRanges)
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
    newscore <- calcIdentity(aaRangesSmall)*width(aaRanges[1])
    if(score<newscore[2,1]){
      score=newscore[2,1]
      bestAA <- aaRanges[c(2,1)]
      bestAAsmall <- aaRangesSmall
      bestDna <- dnaSet[c(2,1)]
    }else if(score==newscore[2,1]){
      if(length(bestAA[[2]])<length(aaRanges[[1]])){
        score=newscore[2,1]
        bestAA <- aaRanges[c(2,1)]
        bestAAsmall <- aaRangesSmall
        bestDna <- dnaSet[c(2,1)]
      }
    }
  }
 # toc()

  if(is.null(bestDna) | is.null(bestAA) | is.null(bestAAsmall)){
    return(NULL)
  }else{
    return(list(dna=bestDna,aa=bestAA, aaOverlap=bestAAsmall))
  }
}
# rm(bestDna)
# rm(bestAA)
# rm(bestAAsmall)
