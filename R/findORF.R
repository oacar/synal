#'Find orf start and stop coordinates
#'@param AAString instance
#'@return coordinates of M---X sequence if exists and longer than I---X, returns I---X if it exists and Longer than M--X
#'@export
findORF<- function(aa_seq){
  ##input: aminoacid sequence and specific aminoacid character to be searched for (M for start or X for stop codon)
  ##returns M---X sequence if exists and longer than I---X, returns I---X if it exists and Longer than M--X
  ## warning: does not checks for coordinates of the character or multiple occurences
  #use regexpr to find M----X sequence
  #attr(gregexpr("M[^X]+X",aa_alignment[[3]]))
  check=FALSE
  mx=regexpr("M[^X]*X",aa_seq)[[1]][1]
  mLength<-attr(regexpr("M[^X]*X",aa_seq), "match.length")[1]
  ix=regexpr("I[^X]*X",aa_seq)[[1]][1]
  iLength<-attr(regexpr("I[^X]*X",aa_seq), "match.length")[1]
  if(mx!=-1){
    if(ix!=-1 && mLength<iLength){
      return(list(ix, ix+iLength-1))
    }else{
      return(list(mx, mx+mLength-1))
    }
  }else{
    return(list(ix, ix+iLength-1))
  }

}
