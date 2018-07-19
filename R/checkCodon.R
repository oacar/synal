
#'Check aminoacid character
#'@param aa_seq Aminoacid sequence as character string. Change it by as.character if using AAString
#'@param checkCharacter the single letter aminoacid code for the searched character. 'M' for start and 'X' for stop
#'@return the first occurence of the character.
#'@export
checkCodon<- function(aa_seq, checkCharacter){

  ##input: aminoacid sequence and specific aminoacid character to be searched for (M for start or X for stop codon)
  ##returns POSITION of the first occurence of the char if the checkCharacter exists on the sequence, if not returns  '-1'
  ## warning: does not check multiple occurences
  pos=-1
  for(i in 1:nchar(as.character(aa_seq))){
    if(as.character(aa_seq[i])==checkCharacter){
      pos=i
      if(pos==i){
        break
      }
    }
  }
  return(pos)
}
