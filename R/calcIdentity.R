
#'Calculate %ID between sequences of alignments
#'Only %ID between Scer and other sequences are calculated
#'Gaps are treated as difference if one sequence has gap and other does not. If both sequence have gap in the same position it will not be count as similarity
#'@param alignment DNAStringSet or AAStringSet object of alignments
#'@return nx1 matrix with with %IDs

#'@export

calcIdentity<-function(alignment){
  #this identity does not exclude gaps !!!!!!!
  #calculate identity matrix with percentage
  counts<-c()
  
  for(i in 2:length(alignment)){
    c<-0
    for(j in 1:length(alignment[[1]])){
      if((alignment[[1]][j]==alignment[[i]][j]) && as.character(alignment[[1]][j])!='-'){
        c<-c+1
      }
    }
    counts<-append(counts,c)
  }
  Ref<-nchar(turnWoGaps(alignment[[1]]))
  counts<-c(Ref,counts)
  identity_percentage<-as.matrix((counts)/Ref)
  return(identity_percentage)
}