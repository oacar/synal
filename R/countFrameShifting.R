#'Count frameshifts with MACSE
#'@param aaSeq macse aligned AA sequence
#'@return number of '!' characters since MACSE uses them to introduce frameshifting indels to the alignments
#'@import stringr
#'@export

countFrameShifting<-function(aaSeq){
  return (str_count(aaSeq, '!'))
}
