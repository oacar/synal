
#'Read syntetic block sequence files
#'@param path folder path that has sequence fasta files in it
#'@param specNames species names to look for
#' every sequence file needs to contain its species name . (i.e. S. Cerevisiae file should contain 'scer')
#' also assumes there is only 1 file for every species
#'@return DNAStringSet object with the sequences in it 
#'order of DNAStringSet will be order of specNames
#'
#'TODO check for multiple files having the search string in them
#'@import stringr
#'@export


readFiles<-function(path,specNames){
  if (dir.exists(path) != TRUE) {
    stop("Path does not exist. Give a correct path")
  }
  mySequences <- DNAStringSet()
  fls <- list.files(path)
  for(i in 1:length(specNames)){
    sName <- specNames[i]
    boolRes <- str_detect(fls,sName)
    if(any(boolRes)){
      fileName <- fls[boolRes]
      dna <- readDNAStringSet(paste0(path,'/',fileName))
      names(dna) <- sName
      mySequences <- append(mySequences,dna)
    }
  }
  
  return(mySequences)}
