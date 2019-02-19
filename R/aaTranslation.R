#'AminoAcid translation of the DNAStringSet object
#'@param subalign the subalignment of DNA sequences that will be used for translation as DNAStringSet object
#'@param DNAStr  DNAStringSet object that has the multiple sequence alignment (it is used for just names of the sequences)
#'
#'@return AAStringSet object of translated and aligned AminoAcid sequences if there is no error.
#' if there is error, it will show a message and return FALSE. It should be checked in the main script
#'@export

aaTranslation <-function(subalign, DNAStr){
  out<-tryCatch({
    s_wogaps=''
    DNA_wogaps<-DNAStringSet()
    for(i in 1:length(subalign)){
      for (j in 1:length(subalign[[i]])){
        if(as.character(subalign[[i]][j])!='-'){
          s_wogaps=paste(s_wogaps, as.character(subalign[[i]][j]), sep = '')
        }
      }
      DNA_wogaps<-append(DNA_wogaps, DNAStringSet(s_wogaps))
      s_wogaps=''
    }
    names(DNA_wogaps)=names(DNAStr)
    attr(GENETIC_CODE, "alt_init_codons")=character(0)
    aatr<-suppressWarnings(translate(DNA_wogaps,genetic.code = GENETIC_CODE))
    for(i in 1:length(subalign)){
      if(width(aatr)[i]==0){
        return(FALSE)
      }
    }
    aatr<-chartr('*','X',aatr)
    aa_alignment<-muscle(aatr,quiet = TRUE)
    aa_alignment<-AAStringSet(aa_alignment)

    return(aa_alignment)},
    error=function(cond){
      message('No region on one of sequence found!')
      return(FALSE)
    })
  return(out)
}
