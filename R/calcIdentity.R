
#'Calculate %ID between sequences of alignments
#'Only %ID between Scer and other sequences are calculated
#'Gaps are treated as difference if one sequence has gap and other does not. If both sequence have gap in the same position it will not be count as similarity
#'@param alignment DNAStringSet or AAStringSet object of alignments
#'@param percent boolean for calculating percent ID or number of identical columns
#'@return nx1 matrix with with %IDs

#'@export

calcIdentity<-function(alignment,percent=T){
  #this identity does not exclude gaps !!!!!!!
  #calculate identity matrix with percentage
  counts<-calc_identity(alignment%>%as.vector(),length(alignment))
  Ref<-nchar(turnWoGaps(alignment[[1]]))
  counts<-c(Ref,counts)
  if(percent){
    id<-as.matrix((counts)/Ref)
  }else{
    id <- as.matrix(counts)
  }
  return(id)
}
