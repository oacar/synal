

#'find Codon
#'Searches for start or stop codons
#'For start codon search the input sequence should end right before the position where it alignes with start codon of proto-gene ORF
#'For stop codon search the input sequence should start from the position where it alignes with start codon of proto-gene ORF so for actual position on the DNA alignment (that start - 1 ) should be added the returned value
#'@param sequence DNAString object
#'@param codonSet ('ATG') for start codon or ('TAA', 'TAG', 'TGA') for stop codon search
#'@return index of the A for start codon or index of last letter of stop codon for stop codon search
#'@seealso \code{\link{findBoth}} for usage
#'@export

findCodon<-function(sequence, codonSet){

  if(codonSet!='ATG' && any(codonSet!=c('TAA', 'TAG', 'TGA') )){
    stop('codonSet should be either \'ATG\' or c(\'TAA\', \'TAG\', \'TGA\')')
  }
  min_i<-Inf
  min_<-Inf
  min_lengthGaps<-Inf
  i<-1
  seqWoGaps<-turnWoGaps(sequence)

  if(any(codonSet=='ATG')){
    allCoordinates<-gregexpr(regSeq(codonSet),seqWoGaps)[[1]]
    allCoordinatesWGaps<-gregexpr(regSeq(codonSet),sequence)[[1]]
    while(i<=length(allCoordinates)){
      diff_<-nchar(seqWoGaps)+1-allCoordinates[i]
      if((diff_)%%3==0 && diff_<min_){
        min_<-diff_
        min_i<-allCoordinatesWGaps[i]

      }
      i=i+1
    }
    return(min_i)

  } else {
    for( codon_ in codonSet){
      allCoordinates<-gregexpr(regSeq(codon_),seqWoGaps)[[1]]
      allCoordinates_withGaps<-gregexpr(regSeq(codon_),sequence)[[1]]
      i<-length(allCoordinates)
      while(i>0){
        diff_<-allCoordinates[i]-1
        if((diff_)%%3==0 && diff_<min_ && diff_>length(seqWoGaps)){
          min_<-diff_
          min_index<-i
          min_i<-allCoordinates_withGaps[i]
          min_lengthGaps<-attr(allCoordinates_withGaps,'match.length')[min_index]
        }
        i=i-1
      }
    }
    return(min_i+min_lengthGaps-1)
  }}
