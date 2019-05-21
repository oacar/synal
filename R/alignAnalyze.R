#' This will be THE function for alignment-analysis output
#' @param filename alignment file name or directory name
#' @param outputDirectory directory to save alignment output
#' @param orfName identifier for ORF. SGD name if annotated
#' @param annotated boolean shows if the given sequence is annotated or not
#' @param specNames species names given by user provided tree. it is needed for correct output names and analysis
#' @param ... placeholder for sequence file names
#' @return save data into csv and return dataframe
#' @importFrom utils file_test write.csv
#' @export


alignAnalyze <- function(filename, orfName, annotated = T, outputDirectory = NULL, specNames, ...) {
  # cl <- makeCluster(6)
  # registerDoParallel(cl)

  ygeneSeq <-  findYGeneSeq(orfName)

  # Input check----------
  # checkInput<-file.exists(filename)
  if (file_test("-f", filename)) {
    mySequences <- readDNAStringSet(filename)
  } else if (file_test("-d", filename)) {
    mySequences <- readFiles(filename)
  } else {
    stop("Filename is wrong. Make sure it is a directory or alignment file that exists")
  }
  # Input read------
  vec <- speciesVector(names(mySequences), specNames)



  # Alignment--------
  dnalist <- align(mySequences, orfName, outputDirectory, ygeneSeq)
  DNAStr <- dnalist$dnaAlignmentList[[1]]
  start <- dnalist$dnaAlignmentList[[2]]
  stop <- dnalist$dnaAlignmentList[[3]]
 # aa_alignment <- dnalist$aa_alignment

  # overlapping orfs----------------


  types <- specNames # c('scer','Spar','Smik','Skud','Sbay','Sarb')
  types <- types[vec > 0]
  findHomolog(DNAStr, start, stop, ygeneSeq, types, outputDirectory, orfName)

  dataTable <- analyze(outputDirectory, orfName, types)
  if (is.null(outputDirectory) == F) write.csv(dataTable, file = paste0(outputDirectory, "/", orfName, "_data.csv"))
  dataTable
}
