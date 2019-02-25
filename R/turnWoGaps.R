
#'Remove gaps (\'-\')
#'@param sequence DNAString object instance
#'@return character string of the sequence without gaps
#'@import stringr
#'@export
turnWoGaps <- function(sequence) {
  s_wogaps<-str_replace_all(sequence,'-','')
  # for (j in 1:length(sequence)){
  #   if(as.character(sequence[j])!='-'){
  #     s_wogaps=paste(s_wogaps, as.character(sequence[j]), sep = '')
  #   }
  # }
  return(s_wogaps)
}
