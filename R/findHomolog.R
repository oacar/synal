#'This function compares all overlapping ORFs with the Scer ORF and saves best ORF in terms of number of AAs (highest number of identical aminoacids)
#'@param DNAStr aligned DNAStringSet with all sequences (including the ORF to be searched) in it, this is the full syntenic block not subalignment
#'@param aa_alignment Amino acid alignment over the subalignment
#'@param start start codon position of ygeneSeq in DNAStr
#'@param stop last nucleotide position of ygeneSeq in DNAStr
#'@param ygeneSeq sequence of interest
#'@param types this is a vector which contains the species names in the alignment for writing correct output files
#'@param orfName ORF identifier name that will be used for file writing
#'@param path path for files to be written in
#'@return nothing returns at the moment, all files are written in path with 3 files for each species. pairwise nucleotide alignment, pairwise AA alignment and pairwise AA alignment that shows only overlap of two ORFs
#'@export

findHomolog <- function(DNAStr, aa_alignment, start, stop, ygeneSeq, types, path, orfName) {
  for(j in 2:(length(DNAStr))){
    r=400
    r=ifelse(start<r,start-1,r)
    bo <- findBestOverlap(DNAStr, j, r, start, stop, ygeneSeq, types)

    if(is.null(bo)==F){
      writeXStringSet(bo$dna, file=paste(paste(path,orfName, sep="/"),"_subalignment_",types[j],".fa",sep = ""))
      writeXStringSet(bo$aa,file=paste(paste(path,orfName, sep="/"),"_AATranslation_",types[j],".fa",sep = ""))
      writeXStringSet(bo$aaOverlap,file=paste(paste(path,orfName, sep="/"),"_AATranslation_overlap_",types[j],".fa",sep = ""))

    }
  }
}
