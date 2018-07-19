
#'Correct macse
#'Changes some macse spesific characters for the general purposes of this package
#'since macse reorder the sequences after alignment, I use previously used alignment file to change the order back to scer,spar,smik,skud,sbay
#'changes '!' character to '-' and '*' to 'X'
#'@param macseAlignment macse aligned AAStringSet object
#'@param seqSet previously used DNAStringSet or AAStringSet object of the same sequences
#'@return corrected AAStringSet object 
#'@export

macseAdjustment<-function(macseAlignment, seqSet){
  aaSet<-AAStringSet()
  namesCorrect<-names(seqSet)
  namesWrong<-names(macseAlignment)
  for (i in 1:length(namesCorrect)){
    for(j in 1:length(namesWrong)){
      if(namesCorrect[i]==namesWrong[j]){
        aaSet<-append(aaSet, macseAlignment[j])
      }
    }
  }
  aaSet<-chartr('!','-',aaSet)
  aaSet<-chartr('*','X',aaSet)
  
  
  return(aaSet)
}
