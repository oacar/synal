#'Analyze msa
#'This function reads all output files previously created and analyzes the outputs
#'@param orfName orf identifier for the orf of interest
#'@param outputDirectory folder path for reading input files
#'@param types this is a vector to assess which species sequences exists
#'@param best specifies whether to analyze only best overlapping orf or all.
#'@importFrom dplyr bind_cols
#'@export

analyze <- function(outputDirectory,orfName,types,best=T) {
  #Analysis------------
  if(dir.exists(outputDirectory)==F){
    stop(paste0(outputDirectory, ' does not exist. Please check the input'))
  }
  #print(orfName)

  colNames<-c("ORF Name",	"Do all Start Codon Align?" ,"Do all Stop Codon Align","Length of ORF (DNA)","Length of Amino Acid Sequence ORF")
  df <- data.frame(matrix(ncol = length(colNames), nrow = 1))
  colnames(df) <- colNames
  subalign<-readDNAStringSet(paste(outputDirectory,'/',orfName, '_subalignment.fa',sep=''))
  aa<-readAAStringSet(paste(outputDirectory,'/',orfName,'_AATranslation.fa',sep = ''))
  vec <- speciesVector(names(aa),types)
  types <- types[vec>0]
  rw=1
  df$`ORF Name`[rw]<-orfName

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


  df$`Do all Start Codon Align?`[rw]<-checkAlignment(startAA)
  df$`Do all Stop Codon Align`[rw]<-checkAlignment(stopAA)
  df$`Length of ORF (DNA)`[rw]<-nchar(turnWoGaps(subalign[[1]]))
  df$`Length of Amino Acid Sequence ORF`[rw]<-nchar(turnWoGaps(aa[[1]]))

  if(best){
    for(i in 2:length(types)){
      typeName <- types[i]
      dataTable <- analyzePairwise(typeName, outputDirectory, orfName,startAA, startDNA, stopAA, subalign[c(1,i)], aa[c(1,i)],df)
      df <- bind_cols(df,dataTable)
    }
  }else{
    for(i in 2:length(types)){
      typeName <- types[i]

      ids <- list.files(paste0(outputDirectory,'/',typeName))%>%str_split('_')%>%sapply(tail,1)%>%unique()%>%str_sub(1,str_length(.)-3)
      numOfOrfs <- length(ids)
      for(j in 1:numOfOrfs){
        dataTable <- analyzePairwise(typeName, outputDirectory, orfName,startAA, startDNA, stopAA, subalign[c(1,i)], aa[c(1,i)],df,ids[j])
        if(is.null(dataTable)==F){
          df <- bind_cols(df,dataTable)
        }
      }

    }
  }


  return(df)
}
