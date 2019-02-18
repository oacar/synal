#'Find orf start and stop coordinates
#'@param AAString instance or character vector
#'@return coordinates of M---X or I---X. It calcualtes longest orf by considering all M to X pairs and first finds all non-overlapping ORFs and then returns the longest one.
#'@export
findORF<- function(aa_seq,findI=F){
  ##input: aminoacid sequence and specific aminoacid character to be searched for (M for start or X for stop codon)
  ##returns M---X sequence if exists and longer than I---X, returns I---X if it exists and Longer than M--X
  ## warning: does not checks for coordinates of the character or multiple occurences
  #use regexpr to find M----X sequence
  #attr(gregexpr("M[^X]+X",aa_alignment[[3]]))
  check=FALSE
  # mx=regexpr("M[^X]*X",aa_seq)[[1]][1]
  # mLength<-attr(regexpr("M[^X]*X",aa_seq), "match.length")[1]
  # ix=regexpr("I[^X]*X",aa_seq)[[1]][1]
  # iLength<-attr(regexpr("I[^X]*X",aa_seq), "match.length")[1]
  if(findI){
    mxs <- str_locate_all(aa_seq,'M|I')[[1]]

  }else{
    mxs <- str_locate_all(aa_seq,'M')[[1]]
  }
  xs <- str_locate_all(aa_seq,'X')[[1]]

  orfRanges <- IRanges()
  for(i in mxs[,1]){
    for(j in xs[,1]){
      if(i<j){
          if(str_count(subseq(aa_seq,start = i,end=j),'X')==1 & any(end(orfRanges)[start(orfRanges)==i]<j)==F  & any(start(orfRanges)[end(orfRanges)==j]<i)==F){
           orfRanges <- append(orfRanges,IRanges(i,j))
        }
      }
     }
  }
  max=0
  for(i in 1:length(orfRanges)){
    aatr <- turnWoGaps(subseq(aa_seq,start = start(orfRanges)[i],end=end(orfRanges)[i]))
    if(nchar(aatr)>max){
      max=nchar(aatr)
      longest <- orfRanges[i]
    }

  }
  return(list(start(longest),end(longest)))
  # if(mx!=-1){
  #   if(ix!=-1 && mLength<iLength){
  #     return(list(ix, ix+iLength-1))
  #   }else{
  #     return(list(mx, mx+mLength-1))
  #   }
  # }else{
  #   return(list(ix, ix+iLength-1))
  # }

}
