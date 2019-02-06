
#'Check Orf
#'
#'Checks if Amino Acid sequence has M --- X or I --- X ORF by using regular expressions.
#'This can be changed to only check M --- X ORF, or for any other start codons
#'
#'@param aa_seq Aminoacid sequence as character string. Change it by as.character if using AAString
#'@return if there is an ORF then return TRUE
#' If there is no ORF check what is missing an M or X
#' If M is missing return 'M'  --- if X is missing return 'X' --- if both are missing return 'BOTH', if there is not orf but X ---- M sequence, it will return 'Neither'
#'@export

checkORF<- function(aa_seq){
  check=FALSE
  mx=gregexpr("M[^X]*X",aa_seq)[[1]][1]
  ix=-1#gregexpr("I[^X]*X",aa_seq)[[1]][1]
  if(mx!=-1 || ix!=-1){
    check=TRUE
  }


  if(check==TRUE){
    return(check)
  }
  else{
    startCheck<-grepl('M',aa_seq)
    startCheckI<-F#grepl('I',aa_seq)
    stopCheck<-grepl('X',aa_seq)
    if(startCheck==FALSE && startCheckI==FALSE && stopCheck==FALSE){
      return('BOTH')
    }else if(startCheck==FALSE && startCheckI==FALSE && stopCheck){
      return('M')
    }else if(stopCheck==FALSE && (startCheck || startCheckI)){
      return('X')
    }else{
      return('Neither')
    }
  }}
