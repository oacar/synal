##### !/usr/bin/Rscript --vanilla

library(optparse)
suppressPackageStartupMessages(library(Biostrings))
library(synal)


option_list <- list(
  make_option(c("-p", "--specNames"),
    action = "store", default = NA, type = "character",
    help = "file that contains species names."
  ), # explainfile format
  make_option(c("-f", "--filename"),
    action = "store", default = NA, type = "character",
    help = "File name for multiple alignments"
  ),
  make_option(c("-v", "--verbose"),
    action = "store_true", default = TRUE,
    help = "Should the program print extra stuff out? [default %default]"
  ),
  make_option(c("-q", "--quiet"),
    action = "store_false", dest = "verbose",
    help = "Make the program not be verbose."
  ),
  make_option(c("-a", "--annotated"),
    action = "store_true", default = FALSE,
    help = "If the gene of interest is annotated, give annotation"
  ),
  make_option(c("-s", "--sequence"),
    action = "store", default = NA, type = "character",
    help = "File name for orf sequence"
  ),
  make_option(c("-n", "--name"),
    action = "store", default = NA, type = "character",
    help = "If annotated is true this should be yeast gene name, otherwise an identifier"
  ),

  make_option(c("-o", "--output"),
    action = "store", default = NA, type = "character",
    help = "Output directory for files to be written"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.na(opt$specNames)) {
  stop(paste0("File containing species name for ", opt$name, " should be given. Stopping"), call. = F)
} else if (file.access(opt$specNames, mode = 4) == -1) {
  stop(paste0(opt$specNames, " doesn't exists or you don't have read permissions. Stopping"), call. = F)
} else {
  specNames <- opt$specNames
}
if (is.na(opt$filename)) {
  stop(paste0("Alignment file for ", opt$name, " should be given. Stopping"), call. = F)
} else if (file.access(opt$filename, mode = 4) == -1) {
  stop(paste0(opt$filename, " doesn't exists or you don't have read permissions. Stopping"), call. = F)
} else {
  aln <- readDNAStringSet(opt$filename)
}
if (opt$annotated) {
  # check if given name is in yeast genome
  seq <- tryCatch(findYGeneSeq(opt$name), warning = function(e) stop(paste0(opt$name, " cannot be found. Stopping"), call. = F))
} else {
  if (is.na(opt$sequence)) {
    stop(paste0("Sequence file for ", opt$name, " should be given. Stopping"), call. = F)
  } else if (file.access(opt$sequence, mode = 4) == -1) {
    stop(paste0(opt$sequence, " cannot be found, please provide correct file name for ", opt$name, ". Stopping"), call. = F)
  } else {
    seq <- readDNAStringSet(opt$sequence)
    names(seq) <- opt$name
  }
}

if (is.na(opt$output)) {
  stop(paste0("Output folder should be given for ", opt$name, ". Stopping"), call. = F)
} else if (file.access(opt$output, mode = 2) == -1) {
  stop(paste0(opt$output, " doesn't exists or you don't have write permissions. Stopping"), call. = F)
} else {
  outputDirectory <- opt$output
}

specNames <- readLines(specNames)
vec <- speciesVector(names(aln), specNames)
orfName <- opt$name


# Alignment--------
dnalist <- align(aln, orfName, outputDirectory, seq)
DNAStr <- dnalist$dnaAlignmentList[[1]]
start <- dnalist$dnaAlignmentList[[2]]
stop <- dnalist$dnaAlignmentList[[3]]
aa_alignment <- dnalist$aa_alignment

# overlapping orfs----------------


types <- specNames # c('scer','Spar','Smik','Skud','Sbay','Sarb')
types <- types[vec > 0]
homologs <- findHomolog(DNAStr, start, stop, seq, types, outputDirectory, orfName, best = TRUE)

dataTable <- analyze(outputDirectory, orfName, types)
if (is.null(outputDirectory) == F) write.csv(dataTable, file = paste0(outputDirectory, "/", orfName, "_data.csv"))
# dataTable
