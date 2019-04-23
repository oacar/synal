#'Finds all overlapping ORFs
#'@param dna DNA sequence that will be searched for finding ORFs
#'@param range the range of interest. Only ORF ranges overlapping with this range will be considered
#'@return IRanges object with start and end positions determined using \code{dna} object
#'@export
findOverlappingOrfs <- function(dna,range) {
  #tic()
  startpositions <- str_locate_all(as.character(dna),'ATG')[[1]]
  if(length(startpositions)==0){
    return(NA)
  }
  stoppositions <- str_locate_all(as.character(dna),'TAA|TAG|TGA')[[1]]

  orfpos <- str_locate_all(as.character(dna),'ATG(?:[ATGC]{3})*?(?:TAA|TAG|TGA)')[[1]]
  ranges <- IRanges()
  #
  eg <- expand.grid(startpositions[,1],stoppositions[,2])
  eg_list <- eg[(eg[,1]<eg[,2]) & (eg[,2]-eg[,1]+1)%%3==0,] #take start positions < stop positions only and create iranges object with those positions
  eg_range <- IRanges(start=eg_list[,1],end=eg_list[,2])
  eg_ranges <- eg_range[eg_range%over%range]
  newranges <- IRanges()

  for(i in 1:length(unique(start(eg_ranges)))){
    e <- unique(start(eg_ranges))[i]
    eth <- eg_ranges[start(eg_ranges)==e]
    maxeth <- eth[which.min(width(eth))]
    newranges <- append(newranges,maxeth)
  }
  newranges2 <- IRanges()
  #from the ORFs that end at the same position, take the one started earlier.
  for(i in 1:length(unique(end(newranges)))){
    e <- unique(end(newranges))[i]
    eth <- newranges[end(newranges)==e]
    maxeth <- eth[which.max(width(eth))]
    newranges2 <- append(newranges2,maxeth)
  }
  #toc()
  newranges2 <- newranges2[width(overlapsRanges(newranges2,range))>12]

  newranges2
}
