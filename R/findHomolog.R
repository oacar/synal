#'This function compares all overlapping ORFs with the Scer ORF and saves best ORF in terms of number of AAs (highest number of identical aminoacids)
#'@param DNAStr aligned DNAStringSet with all sequences (including the ORF to be searched) in it, this is the full syntenic block not subalignment
#'@param aa_alignment Amino acid alignment over the subalignment
#'@param start start codon position of ygeneSeq in DNAStr
#'@param stop last nucleotide position of ygeneSeq in DNAStr
#'@param ygeneSeq sequence of interest
#'@param types this is a vector which contains the species names in the alignment for writing correct output files
#'@param orfName ORF identifier name that will be used for file writing
#'@param outputDirectory path for files to be written in
#'@param best should the function returns all overlapping ORFs or only the best one
#'@return nothing returns at the moment, all files are written in path with 3 files for each species. pairwise nucleotide alignment, pairwise AA alignment and pairwise AA alignment that shows only overlap of two ORFs
#'@export

findHomolog <- function(DNAStr, aa_alignment, start, stop, ygeneSeq, types, outputDirectory=NULL, orfName,best=FALSE) {
  all <- list()
  map_ygene <- map_alignment_sequence(DNAStr[[1]]%>%as.character(),turnWoGaps(DNAStr[[1]]%>%as.character()))
  for(j in 2:(length(DNAStr))){

    bo <- findBestOverlap(DNAStr, j, start, stop, ygeneSeq, types,map_ygene)

    if(is.null(bo)==F){
      if(best){#if user asks only for only best overlap, it will be returned/written to file
        if(is.null(outputDirectory)==F){
          if(dir.exists(outputDirectory)){
            dir.create(paste0(outputDirectory,'/',types[j]))
            bestId <- bo$id
            writeXStringSet(bo$seq[[bestId]]$dna, filepath=paste(paste(outputDirectory,types[j],orfName, sep="/"),"_subalignment_",types[j],"_best.fa",sep = ""))
            writeXStringSet(bo$seq[[bestId]]$aa,filepath=paste(paste(outputDirectory,types[j],orfName, sep="/"),"_AATranslation_",types[j],"_best.fa",sep = ""))
            writeXStringSet(bo$seq[[bestId]]$aaOverlap,filepath=paste(paste(outputDirectory,types[j],orfName, sep="/"),"_AATranslation_overlap_",types[j],"_best.fa",sep = ""))
            writeXStringSet(bo$seq[[bestId]]$dnaOverlap,filepath=paste(paste(outputDirectory,types[j],orfName, sep="/"),"_subalignment_overlap_",types[j],"_best.fa",sep = ""))
            writeXStringSet(bo$seq[[bestId]]$orfAA,filepath=paste(paste(outputDirectory,types[j],orfName, sep="/"),"_orf_aa_",types[j],"_best.fa",sep = ""))

          }else{
            stop(paste0('outputDirectory for writing homologs (in findHomolog) is wrong/does not exist. Please check: ',outputDirectory))
          }}
          all[[types[j]]] <- bo$seq[[bestId]]

      }else{#all overlapping ORFs are returned otherwise
        if(is.null(outputDirectory)==F){
          if(dir.exists(outputDirectory)){
            dir.create(paste0(outputDirectory,'/',types[j]))
            for(itr in 1:length(bo$seq)){
              if(is.null(bo$seq[[itr]])) next
              if(itr==bo$id){#add best to file name for best homolog
                writeXStringSet(bo$seq[[itr]]$dna, filepath=paste(paste(outputDirectory,types[j],orfName, sep="/"),"_subalignment_",types[j],"_best.fa",sep = ""))
                writeXStringSet(bo$seq[[itr]]$aa,filepath=paste(paste(outputDirectory,types[j],orfName, sep="/"),"_AATranslation_",types[j],"_best.fa",sep = ""))
                writeXStringSet(bo$seq[[itr]]$aaOverlap,filepath=paste(paste(outputDirectory,types[j],orfName, sep="/"),"_AATranslation_overlap_",types[j],"_best.fa",sep = ""))
                writeXStringSet(bo$seq[[itr]]$dnaOverlap,filepath=paste(paste(outputDirectory,types[j],orfName, sep="/"),"_subalignment_overlap_",types[j],"_best.fa",sep = ""))
                writeXStringSet(bo$seq[[itr]]$orfAA,filepath=paste(paste(outputDirectory,types[j],orfName, sep="/"),"_orf_aa_",types[j],"_best.fa",sep = ""))
              }else{
                writeXStringSet(bo$seq[[itr]]$dna, filepath=paste(paste(outputDirectory,types[j],orfName, sep="/"),"_subalignment_",types[j],"_",itr,".fa",sep = ""))
                writeXStringSet(bo$seq[[itr]]$aa,filepath=paste(paste(outputDirectory,types[j],orfName, sep="/"),"_AATranslation_",types[j],"_",itr,".fa",sep = ""))
                writeXStringSet(bo$seq[[itr]]$aaOverlap,filepath=paste(paste(outputDirectory,types[j],orfName, sep="/"),"_AATranslation_overlap_",types[j],"_",itr,".fa",sep = ""))
                writeXStringSet(bo$seq[[itr]]$dnaOverlap,filepath=paste(paste(outputDirectory,types[j],orfName, sep="/"),"_subalignment_overlap_",types[j],"_",itr,".fa",sep = ""))
                writeXStringSet(bo$seq[[itr]]$orfAA,filepath=paste(paste(outputDirectory,types[j],orfName, sep="/"),"_orf_aa_",types[j],"_",itr,".fa",sep = ""))
              }

            }

          }else{
            stop(paste0('outputDirectory for writing homologs (in findHomolog) is wrong/does not exist. Please check: ',outputDirectory))
          }}
          all[[types[j]]] <- bo$seq

      }

    }
  }
  all
}
