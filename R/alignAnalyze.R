#'This will be THE function for alignment-analysis output
#'@param filename alignment file name or directory name
#'@param outputDirectory directory to save alignment output
#'@param orfName identifier for ORF. SGD name if annotated
#'@param annotated boolean shows if the given sequence is annotated or not
#'@param specNames species names given by user provided tree. it is needed for correct output names and analysis
#'@param ... placeholder for sequence file names
#'
#'
#'@export


alignAnalyze <- function(filename, orfName, annotated=T,outputDirectory,specNames,...){
  cl <- makeCluster(6)
  registerDoParallel(cl)

  ygeneSeq <-  findYGeneSeq(orfName)

  #Input check----------
  #checkInput<-file.exists(filename)
  if(file_test("-f",filename)){
    mySequences<-readDNAStringSet(filename)
    path<- outputDirectory#str_remove(filename,paste(orfName,'_alignment.fa',sep=''))

  }else if(file_test("-d",filename)){
    mySequences<-readFiles(filename)
    path<-outputDirectory#filename
  }else{
    stop("Filename is wrong. Make sure it is a directory or alignment file that exists")
  }
  #Input read------
  vec<-speciesVector(names(mySequences),specNames)



  #Alignment--------
  dnalist <- align(mySequences,orfName, path, ygeneSeq)
  DNAStr <- dnalist$dnaAlignmentList[[1]]
  start <- dnalist$dnaAlignmentList[[2]]
  stop <- dnalist$dnaAlignmentList[[3]]
  aa_alignment <- dnalist$aa_alignment

  #overlapping orfs----------------


  types<-specNames#c('scer','Spar','Smik','Skud','Sbay','Sarb')
  types <- types[vec>0]
  findHomolog(DNAStr, aa_alignment, start, stop, ygeneSeq, types, path, orfName)

  dataTable <- analyze(path,orfName,types)

}








