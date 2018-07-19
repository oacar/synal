
#'Export subalingment of the region between start and stop Scer orf
#'@param start the start position of proto-gene on the DNAStr object
#'@param stop the stop position of the proto-gene on the DNAStr object (the position of last nucleotide of the ORF)
#'@param DNAStr aligned DNAStringSet object
#'@return list of the subalignment as DNAStringSet object between the positions (start,stop) and identity percentage w.r.t. first sequence (SCER)
#'
#'@export

subalignment<-function(start, stop, DNAStr){
  if((round(start)!=start || round(stop)!=stop)){
    stop('Start and Stop should be integers!')
  }
  else if(start<1 || stop>length(DNAStr[[1]])){
    stop('Start should be positive and stop should be smaller than length of aligned sequence including gaps')
  }
  else if(start>stop-3){
    stop('Stop should be bigger than Start with at least 3 characters')
  }
  #export subalignment with start and stop codon found above
  subalign<-subseq(DNAStr,start,stop)
  #calculate identity matrix with percentage
  dist_mat<-stringDist(subalign, method = 'hamming')
  dd_<-neditAt(subalign[[1]], subalign)
  Ref<-width(subalign)[1]
  identity_percentage<-as.matrix((Ref-dd_)/Ref)
  return(list(subalign,identity_percentage))
}