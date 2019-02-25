#'Search for stop codon to get the overlapping orf
#'@param DNAStr aligned DNAStringSet object
#'@param speciesIndex the index for which sequence in the DNAStringSet object to be searched for ATG and stop codons
#'@param start the start position of the S. Cerevisiae ORF over the alignment object
#'@param orfName the corresponding orf name. Used for a clear warning message
#'@param stop stop position of the S. Cerevisiae ORF over the alignment object. It is used for when the Start codon not found
#'
#'@return The list of boolean value for which will be used to determine if another subalignment and AA translation will be written and new start and stop positions
#'
#'
#'
#'@export

findLongerOverlap <- function(DNAStr,j,start,stop,dirName, aa){
  library(stringr)

  find.orf.start<-findPrevStart(DNAStr,j,start,stop,dirName)
  find.orf.stop<-findNextStop(DNAStr,j,start,stop,dirName,l=length(aa[[j]]))

  if(is.null(find.orf.start$writeFile) | is.null(find.orf.stop$writeFile)){#check if any of them returned NULL
    if(is.null(find.orf.start$writeFile) && is.null(find.orf.stop$writeFile)){
      return(FALSE)
    }else if(is.null(find.orf.stop$writeFile)==T && find.orf.start$writeFile==T){
      find.orf<-find.orf.start
    }else if(is.null(find.orf.start$writeFile)==T && find.orf.stop$writeFile==T){
      find.orf <- find.orf.stop
    }else{
      return(FALSE)
    }
  }else{#if both worked nicely
    if(find.orf.start$writeFile==T && find.orf.stop$writeFile==T){
      if(compareMandX(aa[[j]])=='M'){
        find.orf<-find.orf.start
      }else{
        find.orf<-find.orf.stop
      }
    }else if(find.orf.start$writeFile==T && find.orf.stop$writeFile==F){
      find.orf<-find.orf.start
    }else{
      find.orf <- find.orf.stop
    }
  }
  find.orf
}

#' this compares if there is X ---- M sequence, which one has larger overlap with Scer orf by checking the length between scer start to X and M to end of scer
#' @param aa_seq aminoacid sequence as AAString
#' @return 'M' or 'X'
#'
#' @export
compareMandX <- function(aa_seq){
  aa_wogaps <- AAString(turnWoGaps(aa_seq))
  l <- length(aa_wogaps)
  x <- findX(aa_wogaps)[1]
  x <- as.integer(str_split(x,',')[[1]])[1]
  m <- findM(aa_wogaps)[1]
  m <- as.integer(str_split(m,',')[[1]])[1]
  #i <- findI(aa_wogaps)[1]
  #i <- as.integer(str_split(i,',')[[1]])[1]
  #assume both m(or i) and x found
  #if(m==-1){ m <- i}
  if(x>l-m){
    return('X')
  }else{
    return('M')
  }
}
