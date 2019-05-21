
#'Align syntenic blocks
#'This function take align syntenic blocks and then searches for proto-gene sequence on the Scer sequence using regular expressions
#'@param sequenceSet The sequences of syntenic blocks as DNAStringSet object. At leas 3 sequences are needed and last sequence should be the proto-gene
#'@param algorithm this function uses msa package for alignment which can use c("ClustalW", "ClustalOmega", "Muscle") as alignment algorithms. I used 'Muscle' as default. There is also standalone muscle package which is also used in this package
#'@param smorfSeq if this is given, the {sequenceSet} is taken as full syntenic blocks, if not the last sequence of the {sequenceSet} is taken as the orf of interest sequence
#'@param aligned if this is TRUE, then the given sequences will be taken as aligned, so no realignment will be applied.
#'@return list of alignment object as DNAStringSet with start and stop positions of proto-gene ORF on the aligned Scer sequence
#'@import Biostrings
#'@import msa
#'@import muscle
#'@export
#'

alignWoSmorf <- function(sequenceSet, algorithm='Muscle',smorfSeq=NULL, aligned=F) {

  ## requirement : at least 2 sequences
  l<-length(sequenceSet)

  if(l<2){
    stop("Sequence Set should have 2 or more sequences")
  }
  if(is.null(smorfSeq)){
    sub<-sequenceSet[1:l-1] ##take first l-1 sequences
    reg<-regSeq(sequenceSet[[l]])

  }else{
    if(is(smorfSeq,'XStringSet')){
      smorfSeq <- smorfSeq[[1]]
    }
    sub <- sequenceSet
    reg<-regSeq(smorfSeq)

  }
if(aligned==F){
  alignment <- tryCatch({alignment<-msa(sub, algorithm, order='input', type = 'dna')
  alignment<-DNAStringSet(alignment)
  alignment<-chartr('N','-',alignment)
  alignment},
  error=function(cond){
    alignment <- muscle(sub,quiet = T)
    alignment <- DNAStringSet(alignment)
    alignment
  })
}else{
  alignment=sub
}



  regexpr_result<-regexpr(reg, alignment[[1]], perl = TRUE)
  start<-regexpr_result[1]
  lenght_reg<-attr(regexpr_result, 'match.length')
  end<-start+lenght_reg-1
  return(list(alignment,start, end))

}
