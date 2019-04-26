#'Analyze msa
#'This function reads all output files previously created and analyzes the outputs
#'@param orfName orf identifier for the orf of interest
#'@param outputDirectory folder path for reading input files
#'@param types this is a vector to assess which sequences exists
#'@export

analyze <- function(outputDirectory,orfName,types) {
  #Analysis------------
  if(dir.exists(outputDirectory)==F){
    stop(paste0(outputDirectory, ' does not exist. Please check the input'))
  }

  colNames<-c("ORF Name",	"Do all Start Codon Align?" ,"Do all Stop Codon Align","Length of ORF (DNA)","Length of Amino Acid Sequence ORF")

  for(i in 2:length(types)){
    name <- types[i]
    cn <- c( paste("What is the Start Codon in ",name, "?",sep=''),paste("What is the Start Codon(DNA) in ",name,"?",sep=''),paste("Does start codon align in ",name,"?",sep=''),
             paste("What is the Stop Codon in ",name,"?",sep=''),
             paste(name," % DNA ID over Smorf Frame",sep=''),	paste("Length of ",name," DNA Sequence over Smorf Frame",sep=''),paste("Number of ",name," Gaps over Smorf",sep=''),
             paste("Is there an ORF in the ",name," Amino Acid Sequence?",sep=''),
             paste("Is length of ORF same in ",name,"?",sep=''),
             paste("Length of ",name," Amino Acid start to finish with Gaps",sep=''),	paste("Length of ",name," Amino Acid Start to Finish without Gaps",sep=''), paste(name, " Length Ratio",sep=''),
             paste(name," % Amino Acid over ",name," Overlap",sep=''),	paste(name, " % Amino Acid over Smorf Frame",sep=''),	paste(name, " Number of identical Amino Acid over ",name," Overlap",sep=''))
    colNames <- append(colNames,cn)
  }

  dataTable<-data.frame(matrix(ncol = length(colNames), nrow = 1))#length(list.files(p_))))
  colnames(dataTable)<-colNames

  #Read inputs####
  print(orfName)
  subalign<-readDNAStringSet(paste(outputDirectory,'/',orfName, '_subalignment.fa',sep=''))
  aa<-readAAStringSet(paste(outputDirectory,'/',orfName,'_AATranslation.fa',sep = ''))
  rw=1
  dataTable$`ORF Name`[rw]<-orfName

  #If your alignment files do not have proto-gene sequence as last sequence, please append an arbitrary sequence to end, almost all functions assumes the last sequence is the proto-gene, thus makes calculations accordingly
  if(aa[1]!=aa[length(aa)]){#this number should be changed according to (number of aligned species+1) if above is the case
    aa<-append(aa,aa[1])
  }
  if(subalign[1]!=subalign[length(subalign)]){#this number should be changed according to (number of aligned species+1) if above is the case
    subalign<-append(subalign,subalign[1])
  }


  startDNA<-startCodon(subalign)
  startAA<-startCodon(aa)
  stopDNA<-stopCodon(subalign)
  stopAA<-stopCodon(aa)


  dataTable$`Do all Start Codon Align?`[rw]<-checkAlignment(startAA)
  dataTable$`Do all Stop Codon Align`[rw]<-checkAlignment(stopAA)
  dataTable$`Length of ORF (DNA)`[rw]<-nchar(turnWoGaps(subalign[[1]]))
  dataTable$`Length of Amino Acid Sequence ORF`[rw]<-nchar(turnWoGaps(aa[[1]]))

  for(i in 2:length(types)){
    name=types[i]

    AAFileName<-paste0(outputDirectory,'/',name,'/',orfName, '_AATranslation_',name,'_best.fa')
    DNAFileName<-paste0(outputDirectory,'/',name,'/',orfName, '_subalignment_',name,'_best.fa')
    AAOverlapFileName<-paste0(outputDirectory,'/',name,'/',orfName, '_AATranslation_overlap_',name,'_best.fa')
    check<-FALSE

    if(file.exists(AAFileName)){
      AA<-readAAStringSet(AAFileName)
      overlap <- readAAStringSet(AAOverlapFileName)
      check<-TRUE
    }


    if(check){
      dataTable[[paste("Is there an ORF in the ",name," Amino Acid Sequence?",sep='')]][rw]<-TRUE
    }else{
      dataTable[[paste("Is there an ORF in the ",name," Amino Acid Sequence?",sep='')]][rw]<-FALSE
      AA<-FALSE
    }

    dataTable[[paste("What is the Start Codon in ",name, "?",sep='')]][rw]<-as.character(startAA[[i]])
    dataTable[[paste("What is the Start Codon(DNA) in ",name,"?",sep='')]][rw]<-as.character(startDNA[[i]])
    dataTable[[paste("Does start codon align in ",name,"?",sep='')]][rw]<-(as.character(startAA[[i]])==as.character(startAA[[1]]))
    dataTable[[paste("What is the Stop Codon in ",name,"?",sep='')]][rw]<-as.character(stopAA[[i]])

    dataTable[[paste(name," % DNA ID over Smorf Frame",sep='')]][rw]<-calcIdentity(subalign)[i]

    dataTable[[paste("Length of ",name," DNA Sequence over Smorf Frame",sep='')]][rw]<-nchar(turnWoGaps(subalign[[i]]))
    
    dataTable[[paste("Number of ",name," Gaps over Smorf",sep='')]][rw]<-length(subalign[[i]])-dataTable[[paste("Length of ",name," DNA Sequence over Smorf Frame",sep='')]][rw]

    if(is.logical(AA)!=TRUE){
      dataTable[[paste("Length of ",name," Amino Acid start to finish with Gaps",sep='')]][rw]<-nchar(AA[[2]])
      dataTable[[paste("Length of ",name," Amino Acid Start to Finish without Gaps",sep='')]][rw]<-nchar(turnWoGaps(AA[[2]]))
      dataTable[[paste(name, " Length Ratio",sep='')]][rw] <- nchar(turnWoGaps(AA[[2]]))/dataTable$`Length of Amino Acid Sequence ORF`[rw]
      dataTable[[paste(name," % Amino Acid over ",name," Overlap",sep='')]][rw]<-calcIdentity(overlap)[2]
      dataTable[[paste(name, " Number of identical Amino Acid over ",name," Overlap",sep='')]][rw] <- calcIdentity(overlap,percent = F)[2]
    }

    dataTable[[paste(name, " % Amino Acid over Smorf Frame",sep='')]][rw]<-calcIdentity(aa)[i]
    dataTable[[paste("Is length of ORF same in ",name,"?",sep='')]][rw]<-dataTable$`Length of Amino Acid Sequence ORF`[rw]==dataTable[[paste("Length of ",name," Amino Acid Start to Finish without Gaps",sep='')]][rw]
  }

  return(dataTable)
}
