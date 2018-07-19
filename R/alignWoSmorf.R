
#'Align syntenic blocks
#'This function take align syntenic blocks and then searches for proto-gene sequence on the Scer sequence using regular expressions
#'@param sequenceSet The sequences of syntenic blocks as DNAStringSet object. At leas 3 sequences are needed and last sequence should be the proto-gene
#'@param algorithm this function uses msa package for alignment which can use c("ClustalW", "ClustalOmega", "Muscle") as alignment algorithms. I used 'Muscle' as default. There is also standalone muscle package which is also used in this package
#'@return list of alignment object as DNAStringSet with start and stop positions of proto-gene ORF on the aligned Scer sequence
#'@import Biostrings
#'@import msa
#'@export
#'

alignWoSmorf <- function(sequenceSet, algorithm='Muscle') {

  ##Input: DNA Sequence set that was read from fasta files, not aligned set object
  ## requirement : at least 3 sequences
  l<-length(sequenceSet)
  if(l<3){
    stop("Sequence Set should have 3 or more sequences")
  }
  sub<-sequenceSet[1:l-1] ##take first 3 sequences
  alignment<-msa(sub, algorithm, order='input', type = 'dna')
  alignment<-DNAStringSet(alignment)
  alignment<-chartr('N','-',alignment)

  reg<-regSeq(sequenceSet[[l]])
  regexpr_result<-regexpr(reg, alignment[[1]], perl = TRUE)
  start<-regexpr_result[1]
  lenght_reg<-attr(regexpr_result, 'match.length')
  end<-start+lenght_reg-1
  return(list(alignment,start, end))

}
