#'Finds best overlapping ORF in terms of the number of identical aminoacids
#'@param DNAStr DNA alignment of syntenic block. This is the full alignment not subalignment
#'@param j this is the iterator for the ID of species of interest on DNAStr and types vector
#'@param start start position of ATG of proto-gene on DNAStr
#'@param stop last nucleotide position of proto-gene on DNAStr
#'@param ygeneSeq proto-gene sequence without gaps. This is used for determining the frame of proto-gene with respect to the aligned other species ORF
#'@param types this is used for names. by default c('Spar','Smik','Skud','Sbay')
#'@param map_ygene this is the map between aligned and original sequence of ORF of interest
#'@return list of three elements
#'@return dna = dna alignment of jth species orf with proto-gene, aa = aa alignment with full jth ORF, aaOverlap = overlapping regions of both ORFs(translated), dnaOverlap=overlapping regions of both ORFs
#'@import IRanges
#'@export

findBestOverlap <- function(DNAStr, j,start, stop, ygeneSeq, types , map_ygene=NULL) {

  # if(save){
  #   if(is.null(outputPath)|dir.exists(outputPath)==F){
  #     stop(paste0('You asked to save all overlapping ORFs but the output path: ',outputPath,' for save is wrong!'))
  #   }
  # }
  subseq <- DNAStr[[j]]%>%as.character()
  #tic('ranges')
  #map_ygene <- map_alignment_sequence(DNAStr[[1]]%>%as.character(),turnWoGaps(DNAStr[[1]]%>%as.character()))
  range_ygene_gapped <- IRanges(start,stop)
  map_j <- map_alignment_sequence(subseq,turnWoGaps(subseq))
  range_ygene_start<- ifelse(length(map_j[map_j==start])!=0,map_j[map_j==start]%>%names()%>%as.integer(),
                                  which.max(map_j[map_j<start])%>%names()%>%as.integer)

  range_ygene_end <- ifelse(length(map_j[map_j==stop])!=0,map_j[map_j==stop]%>%names()%>%as.integer(),
                                       which.max(map_j[map_j<stop])%>%names()%>%as.integer)
  range_ygene=IRanges(range_ygene_start,range_ygene_end)
  ranges <- findOverlappingOrfs(subseq%>%turnWoGaps(),range = range_ygene)#IRanges(start,stop))
  ranges_gapped <- IRanges(map_j[start(ranges)],map_j[end(ranges)])
  if(length(ranges)==0){
    return(NULL)
  }
  # toc()
  # tic('best')
  score=-1
  best <- NULL
  bestAA <- NULL
  bestDna <- NULL
  bestAAsmall <- NULL
  alldata <- list()
  #print(j)

  for(u in 1:length(ranges)){
    ranges_gapped_u <- ranges_gapped[u]
    ranges_u <- ranges[u]
    union_range <- union(ranges_gapped_u,range_ygene_gapped)
    intersect_range <- intersect(ranges_gapped_u,range_ygene_gapped)
    intersectStart <- start(intersect_range)
    unionStart <- start(union_range)
    intersectEnd<- end(intersect_range)
    unionEnd <- end(union_range)


    #check if intersect only has gaps:
    subseq_intersect <- subseq(DNAStr[1], intersectStart,intersectEnd)
    if(turnWoGaps(subseq_intersect)%>%nchar==0) next
    #check ygene start position
    ygene_union_int <- checkFrame(dna=DNAStr[1]%>%turnWoGaps(),map = map_ygene,unionStart,unionEnd,intersectStart,intersectEnd,start)
    other_union_int <- checkFrame(dna=DNAStr[j]%>%turnWoGaps(),map = map_j,unionStart,unionEnd,intersectStart,intersectEnd,start(ranges_gapped_u))

    #check other sequence start position
    intersectDna <- append(ygene_union_int$intersectDna,other_union_int$intersectDna)
    unionDna <-  append(ygene_union_int$unionDna,other_union_int$unionDna)#%>%muscle()%>%DNAStringSet()
    names(intersectDna) <- (types[c(1,j)])
    names(unionDna) <- (types[c(1,j)])
    unionAA <- aaTranslation(unionDna,unionDna)%>%muscle(quiet = T)%>%AAStringSet()
    intersectAA <- aaTranslation(intersectDna,intersectDna)%>%muscle(quiet = T)%>%AAStringSet()
    unionDna <-  unionDna%>%muscle(quiet = T)%>%DNAStringSet()
    intersectDna <- intersectDna%>%muscle(quiet = T)%>%DNAStringSet()
    if(is.null(unionDna) | is.null(unionAA) | is.null(intersectAA)){
      next
    }
    alldata[[u]] <- list(dna=unionDna,aa=unionAA, dnaOverlap=intersectDna,aaOverlap=intersectAA)

  if(is.logical(intersectAA)){
    warning(paste(names(ygeneSeq),types[j],'has a problematic overlapping ORF. Nice to check'))
    next
  }
  newscore <- calcIdentity(intersectAA,percent = F)/(width(ygeneSeq)/3)
  if(score<newscore[2,1]){
    score=newscore[2,1]
    bestId=u
  }else if(score==newscore[2,1]){
    if(length(bestAA[[2]])<length(unionAA[[1]])){
      score=newscore[2,1]
      bestId=u
    }
  }
  }
    return(list(seq=alldata,id=bestId))

}
# rm(bestDna)
# rm(bestAA)
# rm(bestAAsmall)
# start <- start(ranges_gapped_u)
# map <- map_ygene
# dna <- DNAStr[1]%>%turnWoGaps()
#'this function finds correct frame sequences of intersection and union of two Iranges objects
#'@param dna ungapped sequence
#'@param map mapped list of original sequence-->gapped alignment created by {map_alignment_sequence}
#'@param unionStart union range start position
#'@param unionEnd union range last position
#'@param intersectStart intersection range start position
#'@param intersectEnd intersection range end position
#'@param start start of ORF of interest, to compare with unionStart and intersectStart values
#'@return list of correct frame intersection and union dna sequences taken as subsequence of {dna}
#'
#'@export
#'
checkFrame <- function(dna, map,unionStart,unionEnd,intersectStart,intersectEnd,start) {
  start_ungapped <-map[map==start]%>%names()%>%as.integer()

  #this if else does this:
  #first checks if position x has a corresponding map key in the map list. If so, returns the key value which is the position of x
  #on the ungapped sequence. If there is no key, then it is a gap '-' character on the gapped sequence, thus returns the key for the next
  #non-gap character.
  #the goal of this is to find the position of union and intersection ranges and compare those to actual start-codons to conserve the
  #frame for translation. otherwise the actual orf frame might change due to extra gaps in the alignment
  union_start_ungapped_pos<- ifelse(length(map[map==unionStart])!=0,map[map==unionStart]%>%names()%>%as.integer(),
                                           which.min(map[map>unionStart])%>%names()%>%as.integer)

  intersect_start_ungapped_pos <- ifelse(length(map[map==intersectStart])!=0,map[map==intersectStart]%>%names()%>%as.integer(),
                                               which.min(map[map>intersectStart])%>%names()%>%as.integer)

  union_end_ungapped_pos<- ifelse(length(map[map==unionEnd])!=0,map[map==unionEnd]%>%names()%>%as.integer(),
                                    which.min(map[map>unionEnd])%>%names()%>%as.integer)

  intersect_end_ungapped_pos <- ifelse(length(map[map==intersectEnd])!=0,map[map==intersectEnd]%>%names()%>%as.integer(),
                                         which.min(map[map>intersectEnd])%>%names()%>%as.integer)


  diff_union <- (start_ungapped-union_start_ungapped_pos)

  diff_intersect <- (start_ungapped-intersect_start_ungapped_pos)

  if(diff_union%%3==1){
    union_start_ungapped_pos <- (union_start_ungapped_pos+1)
  }else if(diff_union%%3==2){
    union_start_ungapped_pos <- (union_start_ungapped_pos-1)
  }

  if(diff_intersect%%3==1){
    intersect_start_ungapped_pos <- (intersect_start_ungapped_pos-2)
  }else if(diff_intersect%%3==2){
    intersect_start_ungapped_pos <- (intersect_start_ungapped_pos-1)
  }

  intersectDna <- dna%>%subseq(intersect_start_ungapped_pos,intersect_end_ungapped_pos)%>%DNAStringSet()
  unionDna <- dna%>%subseq(union_start_ungapped_pos,union_end_ungapped_pos)%>%DNAStringSet()
  return(list(intersectDna=intersectDna,unionDna=unionDna))
}
