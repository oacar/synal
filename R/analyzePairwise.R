#' Analyze pairwise alignment
#' this is helper function for analyze script. Creates a df for every overlapping orf
#' @param typeName species name
#' @param outputDirectory folder path for reading input files
#' @param orfName identifier for ORF. SGD name if annotated
#' @param startAA aminoacid start codons calculated by \code{startCodon}
#' @param startDNA dna start codons calculated by \code{startCodon}
#' @param stopAA aa stop codons calculated by \code{stopCodon}
#' @param subalign dna level subalignment
#' @param aa aminoacid translation using subalignment
#' @param df df for reference orf values
#' @param identifier identifier for the overlapping orf name. it will be used for file names as well.
#' @return dataFrame object


analyzePairwise <- function(typeName, outputDirectory, orfName,startAA, startDNA, stopAA, subalign, aa,df,identifier='best') {
  rw <- 1
  i <- 2
  cn <- c( 'Identifier',paste("What is the Start Codon in ",typeName, "?",sep=''),paste("What is the Start Codon(DNA) in ",typeName,"?",sep=''),paste("Does start codon align in ",typeName,"?",sep=''),
           paste("What is the Stop Codon in ",typeName,"?",sep=''),
           paste(typeName," % DNA ID over Smorf Frame",sep=''),	paste("Length of ",typeName," DNA Sequence over Smorf Frame",sep=''),paste("Number of ",typeName," Gaps over Smorf",sep=''),
           paste("Is there an ORF in the ",typeName," Amino Acid Sequence?",sep=''),
           paste("Is length of ORF same in ",typeName,"?",sep=''),
           paste("Length of ",typeName," Amino Acid start to finish with Gaps",sep=''),	paste("Length of ",typeName," Amino Acid Start to Finish without Gaps",sep=''), paste(typeName, " Length Ratio",sep=''),
           paste(typeName," % Amino Acid over ",typeName," Overlap",sep=''),	paste(typeName, " % Amino Acid over Smorf Frame",sep=''),	paste(typeName, " Number of identical Amino Acid over ",typeName," Overlap",sep=''))
  #coltypeNames <- append(coltypeNames,cn)
  
  
  dataTable<-data.frame(matrix(ncol = length(cn), nrow = 1))#length(list.files(p_))))
  colnames(dataTable)<-cn
  dataTable[['Identifier']][rw] <- identifier
  #Read inputs####
  
  

  AAFileName<-paste0(outputDirectory,'/',typeName,'/',orfName, '_AATranslation_',typeName,'_',identifier,'.fa')
  DNAFileName<-paste0(outputDirectory,'/',typeName,'/',orfName, '_subalignment_',typeName,'_',identifier,'.fa')
  AAOverlapFileName<-paste0(outputDirectory,'/',typeName,'/',orfName, '_AATranslation_overlap_',typeName,'_',identifier,'.fa')
  check<-FALSE
  
  if(file.exists(AAFileName)){
    AA<-readAAStringSet(AAFileName)
    overlap <- readAAStringSet(AAOverlapFileName)
    check<-TRUE
  }else{
    warning(paste0('File: ',AAFileName, 'does not exists. Skipping...'))
    return(NULL)
  }
  
  
  if(check){
    dataTable[[paste("Is there an ORF in the ",typeName," Amino Acid Sequence?",sep='')]][rw]<-TRUE
  }else{
    dataTable[[paste("Is there an ORF in the ",typeName," Amino Acid Sequence?",sep='')]][rw]<-FALSE
    AA<-FALSE
  }
  
  dataTable[[paste("What is the Start Codon in ",typeName, "?",sep='')]][rw]<-as.character(startAA[[i]])
  dataTable[[paste("What is the Start Codon(DNA) in ",typeName,"?",sep='')]][rw]<-as.character(startDNA[[i]])
  dataTable[[paste("Does start codon align in ",typeName,"?",sep='')]][rw]<-(as.character(startAA[[i]])==as.character(startAA[[1]]))
  dataTable[[paste("What is the Stop Codon in ",typeName,"?",sep='')]][rw]<-as.character(stopAA[[i]])
  
  dataTable[[paste(typeName," % DNA ID over Smorf Frame",sep='')]][rw]<-calcIdentity(subalign)[i]
  
  dataTable[[paste("Length of ",typeName," DNA Sequence over Smorf Frame",sep='')]][rw]<-nchar(turnWoGaps(subalign[[i]]))
  
  dataTable[[paste("Number of ",typeName," Gaps over Smorf",sep='')]][rw]<-length(subalign[[i]])-dataTable[[paste("Length of ",typeName," DNA Sequence over Smorf Frame",sep='')]][rw]
  
  if(is.logical(AA)!=TRUE){
    dataTable[[paste("Length of ",typeName," Amino Acid start to finish with Gaps",sep='')]][rw]<-nchar(AA[[2]])
    dataTable[[paste("Length of ",typeName," Amino Acid Start to Finish without Gaps",sep='')]][rw]<-nchar(turnWoGaps(AA[[2]]))
    dataTable[[paste(typeName, " Length Ratio",sep='')]][rw] <- nchar(turnWoGaps(AA[[2]]))/df$`Length of Amino Acid Sequence ORF`[rw]
    dataTable[[paste(typeName," % Amino Acid over ",typeName," Overlap",sep='')]][rw]<-calcIdentity(overlap)[2]
    dataTable[[paste(typeName, " Number of identical Amino Acid over ",typeName," Overlap",sep='')]][rw] <- calcIdentity(overlap,percent = F)[2]
  }
  
  dataTable[[paste(typeName, " % Amino Acid over Smorf Frame",sep='')]][rw]<-calcIdentity(aa)[i]
  dataTable[[paste("Is length of ORF same in ",typeName,"?",sep='')]][rw]<-df$`Length of Amino Acid Sequence ORF`[rw]==dataTable[[paste("Length of ",typeName," Amino Acid Start to Finish without Gaps",sep='')]][rw]
  dataTable
}
