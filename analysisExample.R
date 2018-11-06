library(openxlsx)
library(Biostrings)
source('functions.R')
#p_ is the directory containing the alignment files' folders
p_<-'data/analysisInput/alignmentFiles//'
fls<-list.files(p_)
colNames<-c("ORF Name",	"Do all Start Codon Align?" ,"Do all Stop Codon Align","Length of ORF",
            "Length of Amino Acid Sequence ORF",	"What is the Start Codon in Para?","What is the Start Codon(DNA) in Para?","What is the Start Codon in Para with Macse?"		,"Does start codon align in Para?","What is the Stop Codon in Para?","What is the Stop Codon in Para with Macse?",
            "Para % DNA ID over Smorf Frame",	"Length of Para DNA Sequence over Smorf Frame",
            "Number of Para Gaps over Smorf",	"Is there an ORF in the Para Amino Acid Sequence?",
            "Does Para ORF overlap with Scer?",	"Is it a different size para?",
            "Length of Para Amino Acid start to finish with Gaps",	"Length of Para Amino Acid Start to Finish without Gaps",
            "Para % Amino Acid over Para frame",	"Para % Amino Acid over Smorf Frame",	"Para % Amino Acid over smorf with Macse","Para % Amino Acid over Para with Macse",
            "Where are other \"M\" in Para over the Smorf",	"Where are other \"X\" in the Para over the Smorf",	"Para Number of FrameShifts",

            "What is the Start Codon in Mik?","What is the Start Codon(DNA) in Mik?","What is the Start Codon in Mik with Macse?"		,"Does start codon align in Mik?","What is the Stop Codon in Mik?","What is the Stop Codon in Mik with Macse?",
            "Mik % DNA ID over Smorf Frame",	"Length of Mik DNA Sequence over Smorf Frame",
            "Number of Mik Gaps over Smorf","Is there an ORF in the Mik Amino Acid Sequence?",
            "Does Mik ORF overlap with Scer?",	"Is it a different size mik?",
            "Length of Mik Amino Acid start to finish with Gaps",	"Length of Mik Amino Acid Start to Finish without Gaps",
            "Mik % Amino Acid over Mik frame",	"Mik % Amino Acid over Smorf Frame",	"Mik % Amino Acid over smorf with Macse","Mik % Amino Acid over Mik with Macse",
            "Where are other \"M\" in Mik over the Smorf",	"Where are other \"X\" in the Mik over the Smorf","Mik Number of FrameShifts",

            "What is the Start Codon in Bay?","What is the Start Codon(DNA) in Bay?","What is the Start Codon in Bay with Macse?"		,"Does start codon align in Bay?","What is the Stop Codon in Bay?","What is the Stop Codon in Bay with Macse?",
            "Bay % DNA ID over Smorf Frame",	"Length of Bay DNA Sequence over Smorf Frame",
            "Number of Bay Gaps over Smorf","Is there an ORF in the Bay Amino Acid Sequence?",
            "Does Bay ORF overlap with Scer?",	"Is it a different size Bay?",
            "Length of Bay Amino Acid start to finish with Gaps",	"Length of Bay Amino Acid Start to Finish without Gaps",
            "Bay % Amino Acid over Bay frame",	"Bay % Amino Acid over Smorf Frame",	"Bay % Amino Acid over smorf with Macse","Bay % Amino Acid over Bay with Macse",
            "Where are other \"M\" in Bay over the Smorf",	"Where are other \"X\" in the Bay over the Smorf","Bay Number of FrameShifts",
            "What is the Start Codon in Kud?","What is the Start Codon(DNA) in Kud?","What is the Start Codon in Kud with Macse?"		,"Does start codon align in Kud?","What is the Stop Codon in Kud?","What is the Stop Codon in Kud with Macse?",
            "Kud % DNA ID over Smorf Frame",	"Length of Kud DNA Sequence over Smorf Frame",
            "Number of Kud Gaps over Smorf","Is there an ORF in the Kud Amino Acid Sequence?",
            "Does Kud ORF overlap with Scer?",	"Is it a different size Kud?",
            "Length of Kud Amino Acid start to finish with Gaps",	"Length of Kud Amino Acid Start to Finish without Gaps",
            "Kud % Amino Acid over Kud frame",	"Kud % Amino Acid over Smorf Frame",	"Kud % Amino Acid over smorf with Macse","Kud % Amino Acid over Kud with Macse",
            "Where are other \"M\" in Kud over the Smorf",	"Where are other \"X\" in the Kud over the Smorf","Kud Number of FrameShifts")
dataTable<-data.frame(matrix(ncol = length(colNames), nrow = length(list.files(p_))))
colnames(dataTable)<-colNames


for ( i in 1:length(fls)){
  #Read inputs####
  dirName<-fls[i]
  path<-paste(p_,'/',dirName,'/',sep = '')
  print(dirName)
  alignment <- readDNAStringSet(paste(path,dirName,'_alignment.fa',sep = ''))
  subalign<-readDNAStringSet(paste(path,dirName, '_subalignment.fa',sep=''))
  aa<-readAAStringSet(paste(path,dirName,'_AATranslation.fa',sep = ''))
  aa_macse.p<-readAAStringSet(paste(path,dirName,'_subalignment_macse_AA.fasta',sep = ''))
  aa_macse<-macseAdjustment(aa_macse.p,subalign)

  paraFrameAAFileName<-paste(path,dirName, '_AATranslation_para.fa',sep = '')
  paraFrameDNAFileName<-paste(path,dirName, '_subalignment_para.fa',sep = '')

  parasecondFrameAAFileName<-paste(path,dirName, '_AATranslation_para_2ndFrame.fa',sep = '')
  parathirdFrameAAFileName<-paste(path,dirName, '_AATranslation_para_3rdFrame.fa',sep = '')

  mikFrameDNAFileName<-paste(path,dirName, '_subalignment_mik.fa',sep = '')
  mikFrameAAFileName<-paste(path,dirName, '_AATranslation_mik.fa',sep = '')
  miksecondFrameAAFileName<-paste(path,dirName, '_AATranslation_mik_2ndFrame.fa',sep = '')
  mikthirdFrameAAFileName<-paste(path,dirName, '_AATranslation_mik_3rdFrame.fa',sep = '')

  bayFrameDNAFileName<-paste(path,dirName, '_subalignment_bay.fa',sep = '')
  bayFrameAAFileName<-paste(path,dirName, '_AATranslation_bay.fa',sep = '')

  kudFrameDNAFileName<-paste(path,dirName, '_subalignment_kud.fa',sep = '')
  kudFrameAAFileName<-paste(path,dirName, '_AATranslation_kud.fa',sep = '')

  overlappingMik<-FALSE
  overlappingPara<-FALSE
  overlappingBay<-FALSE
  overlappingKud<-FALSE
  #these checks are for checking if there is a para(or mik,kud,bay) frame file exists
  #following fileChecks are only for larger para or mik frames
  #at this point smaller para-mik ORFs are not saved but are analyzed within this script
  paraCheck<-FALSE
  mikCheck<-FALSE
  bayCheck<-FALSE
  kudCheck<-FALSE
  if(file.exists(paraFrameAAFileName)){
    paraFrameAA<-readAAStringSet(paraFrameAAFileName)
    paraCheck<-TRUE
  }
  if(file.exists(mikFrameAAFileName)){
    mikFrameAA<-readAAStringSet(mikFrameAAFileName)
    mikCheck<-TRUE
  }

  if(file.exists(kudFrameAAFileName)){
    kudFrameAA<-readAAStringSet(kudFrameAAFileName)
    kudCheck<-TRUE
  }

  if(file.exists(bayFrameAAFileName)){
    bayFrameAA<-readAAStringSet(bayFrameAAFileName)
    bayCheck<-TRUE
  }

  dataTable$`ORF Name`[i]<-dirName
  #If your alignment files do not have proto-gene sequence as last sequence, please append an arbitrary sequence to end, almost all functions assumes the last sequence is the proto-gene, thus makes calculations accordingly
  if(length(aa)==5){#this number should be changed according to (number of aligned species+1) if above is the case
    aa<-append(aa,aa[1])
    subalign<-append(subalign,subalign[1])
    aa_macse<-append(aa_macse,aa_macse[1])
    aa<-append(aa,aa[1])
  }
  startDNA<-startCodon(subalign)
  startAA<-startCodon(aa)
  stopDNA<-stopCodon(subalign)
  stopAA<-stopCodon(aa)
  startMacseAA<-startCodon(aa_macse)
  stopMacseAA<-stopCodon(aa_macse)
  dataTable$`Do all Start Codon Align?`[i]<-checkAlignment(startAA)
  dataTable$`Do all Stop Codon Align`[i]<-checkAlignment(stopAA)
  dataTable$`Length of ORF`[i]<-nchar(turnWoGaps(subalign[[1]]))
  dataTable$`Length of Amino Acid Sequence ORF`[i]<-nchar(turnWoGaps(aa[[1]]))


  #Para analysis####
  paraFrameShifts<-countFrameShifting(aa_macse.p[2])

  if(paraCheck){
    dataTable$`Is there an ORF in the Para Amino Acid Sequence?`[i]<-TRUE
    overlappingPara<-findSmallFrame(aa,'Para',paraFrameAA[[2]])
  }else{
    if(checkORF(aa[[2]])==TRUE){
      dataTable$`Is there an ORF in the Para Amino Acid Sequence?`[i]<-TRUE
      paraFrameAA<-findSmallFrame(aa,'Para')

    }else{
      dataTable$`Is there an ORF in the Para Amino Acid Sequence?`[i]<-FALSE
      paraFrameAA<-FALSE
    }
  }
  dataTable$`What is the Start Codon in Para?`[i]<-as.character(startAA[[2]])
  dataTable$`What is the Start Codon(DNA) in Para?`[i]<-as.character(startDNA[[2]])
  dataTable$`What is the Start Codon in Para with Macse?`[i]<-as.character(startMacseAA[[2]])
  dataTable$`What is the Stop Codon in Para with Macse?`[i]<-as.character(stopMacseAA[[2]])
  dataTable$`Does start codon align in Para?`[i]<-(as.character(startAA[[2]])==as.character(startAA[[1]]))
  dataTable$`What is the Stop Codon in Para?`[i]<-as.character(stopAA[[2]])
  dataTable$`Para % DNA ID over Smorf Frame`[i]<-calcIdentity(subalign)[2]
  dataTable$`Length of Para DNA Sequence over Smorf Frame`[i]<-nchar(turnWoGaps(subalign[[2]]))
  dataTable$`Number of Para Gaps over Smorf`[i]<-length(subalign[[2]])-dataTable$`Length of Para DNA Sequence over Smorf Frame`[i]

  if(is.logical(paraFrameAA)!=TRUE){
    dataTable$`Length of Para Amino Acid start to finish with Gaps`[i]<-nchar(paraFrameAA[[2]])
    dataTable$`Length of Para Amino Acid Start to Finish without Gaps`[i]<-nchar(turnWoGaps(paraFrameAA[[2]]))
    if(is.logical(overlappingPara)){
      dataTable$`Para % Amino Acid over Para frame`[i]<-calcIdentity(paraFrameAA)[2]
    }else{
      dataTable$`Para % Amino Acid over Para frame`[i]<-calcIdentity(overlappingPara)[2]
    }

  }

  dataTable$`Para % Amino Acid over Smorf Frame`[i]<-calcIdentity(aa)[2]
  dataTable$`Para % Amino Acid over smorf with Macse`[i]<-calcIdentity(aa_macse)[2]
  macseParaAA<-findSmallFrame(aa_macse,'Para')
  dataTable$`Para % Amino Acid over Para with Macse`[i]<-ifelse(is.logical(macseParaAA) , NA, calcIdentity(macseParaAA)[2])
  dataTable$`Where are other "M" in Para over the Smorf`[i]<-findM(aa[[2]])
  dataTable$`Where are other "X" in the Para over the Smorf`[i]<- findX(aa[[2]])
  dataTable$`Is it a different size para?`[i]<-dataTable$`Length of Amino Acid Sequence ORF`[i]!=dataTable$`Length of Para Amino Acid Start to Finish without Gaps`[i]
  dataTable$`Para Number of FrameShifts`[i]<-paraFrameShifts

  #Mik analysis####
  mikFrameShifts<-countFrameShifting(aa_macse.p[3])

  if(mikCheck){
    dataTable$`Is there an ORF in the Mik Amino Acid Sequence?`[i]<-TRUE
    overlappingMik<-findSmallFrame(aa,'Mik',mikFrameAA[[3]])

  }else{
    if(checkORF(aa[[3]])==TRUE){
      dataTable$`Is there an ORF in the Mik Amino Acid Sequence?`[i]<-TRUE
      mikFrameAA<-findSmallFrame(aa,'Mik')
    }else{
      dataTable$`Is there an ORF in the Mik Amino Acid Sequence?`[i]<-FALSE
      mikFrameAA<-FALSE
    }
  }
  dataTable$`What is the Start Codon in Mik?`[i]<-as.character(startAA[[3]])
  dataTable$`What is the Start Codon(DNA) in Mik?`[i]<-as.character(startDNA[[3]])
  dataTable$`What is the Start Codon in Mik with Macse?`[i]<-as.character(startMacseAA[[3]])
  dataTable$`What is the Stop Codon in Mik with Macse?`[i]<-as.character(stopMacseAA[[3]])
  dataTable$`Does start codon align in Mik?`[i]<-(as.character(startAA[[3]])==as.character(startAA[[1]]))
  dataTable$`What is the Stop Codon in Mik?`[i]<-as.character(stopAA[[3]])
  dataTable$`Mik % DNA ID over Smorf Frame`[i]<-calcIdentity(subalign)[3]
  dataTable$`Length of Mik DNA Sequence over Smorf Frame`[i]<-nchar(turnWoGaps(subalign[[3]]))
  dataTable$`Number of Mik Gaps over Smorf`[i]<-length(subalign[[3]])-dataTable$`Length of Mik DNA Sequence over Smorf Frame`[i]

  if(is.logical(mikFrameAA)!=TRUE){
    dataTable$`Length of Mik Amino Acid start to finish with Gaps`[i]<-nchar(mikFrameAA[[2]])
    dataTable$`Length of Mik Amino Acid Start to Finish without Gaps`[i]<-nchar(turnWoGaps(mikFrameAA[[2]]))
    if(is.logical(overlappingMik)){
      dataTable$`Mik % Amino Acid over Mik frame`[i]<-calcIdentity(mikFrameAA)[2]
    }else{
      dataTable$`Mik % Amino Acid over Mik frame`[i]<-calcIdentity(overlappingMik)[2]

    }

  }
  dataTable$`Mik % Amino Acid over Smorf Frame`[i]<-calcIdentity(aa)[3]
  dataTable$`Mik % Amino Acid over smorf with Macse`[i]<-calcIdentity(aa_macse)[3]
  macseMikAA<-findSmallFrame(aa_macse,'Mik')
  dataTable$`Mik % Amino Acid over Mik with Macse`[i]<-ifelse(is.logical(macseMikAA) , NA, calcIdentity(macseMikAA)[2])
  dataTable$`Where are other "M" in Mik over the Smorf`[i]<-findM(aa[[3]])
  dataTable$`Where are other "X" in the Mik over the Smorf`[i]<- findX(aa[[3]])
  dataTable$`Mik Number of FrameShifts`[i]<-mikFrameShifts

  dataTable$`Is it a different size mik?`[i]<-dataTable$`Length of Amino Acid Sequence ORF`[i]==dataTable$`Length of Mik Amino Acid Start to Finish without Gaps`[i]

  #Kud analysis####
  if(length(subalign)==5){
    kudFrameShifts<-countFrameShifting(aa_macse.p[4])
    if(kudCheck){
      dataTable$`Is there an ORF in the Kud Amino Acid Sequence?`[i]<-TRUE
      overlappingKud<-findSmallFrame(aa,'Kud', kudFrameAA[[4]])

    }else{
      if(checkORF(aa[[4]])==TRUE){
        dataTable$`Is there an ORF in the Kud Amino Acid Sequence?`[i]<-TRUE
        kudFrameAA<-findSmallFrame(aa,'Kud')
      }else{
        dataTable$`Is there an ORF in the Kud Amino Acid Sequence?`[i]<-FALSE
        kudFrameAA<-FALSE
      }
    }
    dataTable$`What is the Start Codon in Kud?`[i]<-as.character(startAA[[4]])
    dataTable$`What is the Start Codon(DNA) in Kud?`[i]<-as.character(startDNA[[4]])
    dataTable$`What is the Start Codon in Kud with Macse?`[i]<-as.character(startMacseAA[[4]])
    dataTable$`What is the Stop Codon in Kud with Macse?`[i]<-as.character(stopMacseAA[[4]])
    dataTable$`Does start codon align in Kud?`[i]<-(as.character(startAA[[4]])==as.character(startAA[[1]]))
    dataTable$`What is the Stop Codon in Kud?`[i]<-as.character(stopAA[[4]])
    dataTable$`Kud % DNA ID over Smorf Frame`[i]<-calcIdentity(subalign)[4]
    dataTable$`Length of Kud DNA Sequence over Smorf Frame`[i]<-nchar(turnWoGaps(subalign[[4]]))
    dataTable$`Number of Kud Gaps over Smorf`[i]<-length(subalign[[4]])-dataTable$`Length of Kud DNA Sequence over Smorf Frame`[i]

    if(is.logical(kudFrameAA)!=TRUE){
      dataTable$`Length of Kud Amino Acid start to finish with Gaps`[i]<-nchar(kudFrameAA[[2]])
      dataTable$`Length of Kud Amino Acid Start to Finish without Gaps`[i]<-nchar(turnWoGaps(kudFrameAA[[2]]))
      if(is.logical(overlappingKud)){
        dataTable$`Kud % Amino Acid over Kud frame`[i]<-calcIdentity(kudFrameAA)[2]
      }else{
        dataTable$`Kud % Amino Acid over Kud frame`[i]<-calcIdentity(overlappingKud)[2]

      }

    }
    dataTable$`Kud % Amino Acid over Smorf Frame`[i]<-calcIdentity(aa)[4]
    dataTable$`Kud % Amino Acid over smorf with Macse`[i]<-calcIdentity(aa_macse)[4]
    macseKudAA<-findSmallFrame(aa_macse,'Kud')
    dataTable$`Kud % Amino Acid over Kud with Macse`[i]<-ifelse(is.logical(macseKudAA) , NA, calcIdentity(macseKudAA)[2])
    dataTable$`Where are other "M" in Kud over the Smorf`[i]<-findM(aa[[4]])
    dataTable$`Where are other "X" in the Kud over the Smorf`[i]<- findX(aa[[4]])
    dataTable$`Kud Number of FrameShifts`[i]<-kudFrameShifts

    dataTable$`Is it a different size Kud?`[i]<-dataTable$`Length of Amino Acid Sequence ORF`[i]==dataTable$`Length of Kud Amino Acid Start to Finish without Gaps`[i]
  }

  #Bay analysis####

  if(length(subalign)==6){
    bayFrameShifts<-countFrameShifting(aa_macse.p[5])
    if(bayCheck){
      dataTable$`Is there an ORF in the Bay Amino Acid Sequence?`[i]<-TRUE
      overlappingBay<-findSmallFrame(aa,'Bay',bayFrameAA[[5]])

    }else{
      if(checkORF(aa[[5]])==TRUE){
        dataTable$`Is there an ORF in the Bay Amino Acid Sequence?`[i]<-TRUE
        bayFrameAA<-findSmallFrame(aa,'Bay')
      }else{
        dataTable$`Is there an ORF in the Bay Amino Acid Sequence?`[i]<-FALSE
        bayFrameAA<-FALSE
      }
    }
    dataTable$`What is the Start Codon in Bay?`[i]<-as.character(startAA[[5]])
    dataTable$`What is the Start Codon(DNA) in Bay?`[i]<-as.character(startDNA[[5]])
    dataTable$`What is the Start Codon in Bay with Macse?`[i]<-as.character(startMacseAA[[5]])
    dataTable$`What is the Stop Codon in Bay with Macse?`[i]<-as.character(stopMacseAA[[5]])
    dataTable$`Does start codon align in Bay?`[i]<-(as.character(startAA[[5]])==as.character(startAA[[1]]))
    dataTable$`What is the Stop Codon in Bay?`[i]<-as.character(stopAA[[5]])
    dataTable$`Bay % DNA ID over Smorf Frame`[i]<-calcIdentity(subalign)[5]
    dataTable$`Length of Bay DNA Sequence over Smorf Frame`[i]<-nchar(turnWoGaps(subalign[[5]]))
    dataTable$`Number of Bay Gaps over Smorf`[i]<-length(subalign[[5]])-dataTable$`Length of Bay DNA Sequence over Smorf Frame`[i]

    if(is.logical(bayFrameAA)!=TRUE){
      dataTable$`Length of Bay Amino Acid start to finish with Gaps`[i]<-nchar(bayFrameAA[[2]])
      dataTable$`Length of Bay Amino Acid Start to Finish without Gaps`[i]<-nchar(turnWoGaps(bayFrameAA[[2]]))
      if(is.logical(overlappingBay)){
        dataTable$`Bay % Amino Acid over Bay frame`[i]<-calcIdentity(bayFrameAA)[2]
      }else{
        dataTable$`Bay % Amino Acid over Bay frame`[i]<-calcIdentity(overlappingBay)[2]

      }

    }
    dataTable$`Bay % Amino Acid over Smorf Frame`[i]<-calcIdentity(aa)[5]
    dataTable$`Bay % Amino Acid over smorf with Macse`[i]<-calcIdentity(aa_macse)[5]
    macseBayAA<-findSmallFrame(aa_macse,'Bay')
    dataTable$`Bay % Amino Acid over Bay with Macse`[i]<-ifelse(is.logical(macseBayAA) , NA, calcIdentity(macseBayAA)[2])
    dataTable$`Where are other "M" in Bay over the Smorf`[i]<-findM(aa[[5]])
    dataTable$`Where are other "X" in the Bay over the Smorf`[i]<- findX(aa[[5]])
    dataTable$`Bay Number of FrameShifts`[i]<-bayFrameShifts

    dataTable$`Is it a different size bay?`[i]<-dataTable$`Length of Amino Acid Sequence ORF`[i]==dataTable$`Length of Bay Amino Acid Start to Finish without Gaps`[i]

  }
}

write.xlsx(dataTable, 'data.example.xlsx')
