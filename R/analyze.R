#'This function reads all output files previously created and analyzes the outputs
#'@param orfName orf identifier for the orf of interest
#'@param path folder path for output files
#'@param vec this is a vector to assess which sequences exists in this order scer,spar,smik,skud,sbay
#'@export

analyze <- function(orfName, path, vec) {
  #Analysis------------
  types<-c('scer','para','mik','kud','bay','arb')
  types <- types[vec>0]
  #p_ is the directory containing the alignment files' folders
  #p_<-'data/analysisInput/alignmentFiles//'
  #p_ <- '../macbook_data/anne/syntenic_lethal_network/aaronAlignments/y_genes/noProblem/'
  #fls<-list.files(p_)



  colNames<-c("ORF Name",	"Do all Start Codon Align?" ,"Do all Stop Codon Align","Length of ORF",
              "Length of Amino Acid Sequence ORF",	"What is the Start Codon in Para?","What is the Start Codon(DNA) in Para?","What is the Start Codon in Para with Macse?"		,"Does start codon align in Para?","What is the Stop Codon in Para?","What is the Stop Codon in Para with Macse?",
              "Para % DNA ID over Smorf Frame",	"Length of Para DNA Sequence over Smorf Frame",
              "Number of Para Gaps over Smorf",	"Is there an ORF in the Para Amino Acid Sequence?",
              "Is length of ORF same in para?",
              "Length of Para Amino Acid start to finish with Gaps",	"Length of Para Amino Acid Start to Finish without Gaps", "Para Length Ratio",
              "Para % Amino Acid over Para Overlap",	"Para % Amino Acid over Smorf Frame",	"Para Number of identical Amino Acid over Para Overlap",
              "Para % Amino Acid over smorf with Macse","Para % Amino Acid over Para with Macse",
              "Where are other \"M\" in Para over the Smorf",	"Where are other \"X\" in the Para over the Smorf",	"Para Number of FrameShifts",

              "What is the Start Codon in Mik?","What is the Start Codon(DNA) in Mik?","What is the Start Codon in Mik with Macse?","Does start codon align in Mik?","What is the Stop Codon in Mik?","What is the Stop Codon in Mik with Macse?",
              "Mik % DNA ID over Smorf Frame",	"Length of Mik DNA Sequence over Smorf Frame",
              "Number of Mik Gaps over Smorf","Is there an ORF in the Mik Amino Acid Sequence?",
              "Is length of ORF same in mik?",
              "Length of Mik Amino Acid start to finish with Gaps",	"Length of Mik Amino Acid Start to Finish without Gaps","Mik Length Ratio",
              "Mik % Amino Acid over Mik Overlap",	"Mik % Amino Acid over Smorf Frame" ,"Mik Number of identical Amino Acid over Mik Overlap",
              "Mik % Amino Acid over smorf with Macse","Mik % Amino Acid over Mik with Macse",
              "Where are other \"M\" in Mik over the Smorf",	"Where are other \"X\" in the Mik over the Smorf","Mik Number of FrameShifts",

              "What is the Start Codon in Bay?","What is the Start Codon(DNA) in Bay?","What is the Start Codon in Bay with Macse?"		,"Does start codon align in Bay?","What is the Stop Codon in Bay?","What is the Stop Codon in Bay with Macse?",
              "Bay % DNA ID over Smorf Frame",	"Length of Bay DNA Sequence over Smorf Frame",
              "Number of Bay Gaps over Smorf","Is there an ORF in the Bay Amino Acid Sequence?",
              "Is length of ORF same in Bay?",
              "Length of Bay Amino Acid start to finish with Gaps",	"Length of Bay Amino Acid Start to Finish without Gaps","Bay Length Ratio",
              "Bay % Amino Acid over Bay Overlap",	"Bay % Amino Acid over Smorf Frame","Bay Number of identical Amino Acid over Bay Overlap",
              "Bay % Amino Acid over smorf with Macse","Bay % Amino Acid over Bay with Macse",
              "Where are other \"M\" in Bay over the Smorf",	"Where are other \"X\" in the Bay over the Smorf","Bay Number of FrameShifts",

              "What is the Start Codon in Kud?","What is the Start Codon(DNA) in Kud?","What is the Start Codon in Kud with Macse?"		,"Does start codon align in Kud?","What is the Stop Codon in Kud?","What is the Stop Codon in Kud with Macse?",
              "Kud % DNA ID over Smorf Frame",	"Length of Kud DNA Sequence over Smorf Frame",
              "Number of Kud Gaps over Smorf","Is there an ORF in the Kud Amino Acid Sequence?",
              "Is length of ORF same in Kud?",
              "Length of Kud Amino Acid start to finish with Gaps",	"Length of Kud Amino Acid Start to Finish without Gaps","Kud Length Ratio",
              "Kud % Amino Acid over Kud Overlap",	"Kud % Amino Acid over Smorf Frame","Kud Number of identical Amino Acid over Kud Overlap",
              "Kud % Amino Acid over smorf with Macse","Kud % Amino Acid over Kud with Macse",
              "Where are other \"M\" in Kud over the Smorf",	"Where are other \"X\" in the Kud over the Smorf","Kud Number of FrameShifts")
  dataTable<-data.frame(matrix(ncol = length(colNames), nrow = 1))#length(list.files(p_))))
  colnames(dataTable)<-colNames


  #for ( i in 1:length(fls)){
  #Read inputs####
  #orfName<-fls[i]
  #path<-paste(p_,'/',orfName,'/',sep = '')
  print(orfName)
  #alignment <- readDNAStringSet(paste(path,'/',orfName,'_alignment.fa',sep = ''))
  subalign<-readDNAStringSet(paste(path,'/',orfName, '_subalignment.fa',sep=''))
  aa<-readAAStringSet(paste(path,'/',orfName,'_AATranslation.fa',sep = ''))
  # aa_macse.p<-readAAStringSet(paste(path,'/',orfName,'_subalignment_macse_AA.fasta',sep = ''))
  # aa_macse<-macseAdjustment(aa_macse.p,subalign)

  paraFrameAAFileName<-paste(path,'/',orfName, '_AATranslation_Spar.fa',sep = '')
  paraFrameDNAFileName<-paste(path,'/',orfName, '_subalignment_Spar.fa',sep = '')
  paraFrameAAOverlapFileName<-paste(path,'/',orfName, '_AATranslation_overlap_Spar.fa',sep = '')

  mikFrameDNAFileName<-paste(path,'/',orfName, '_subalignment_Smik.fa',sep = '')
  mikFrameAAFileName<-paste(path,'/',orfName, '_AATranslation_Smik.fa',sep = '')
  mikFrameAAOverlapFileName<-paste(path,'/',orfName, '_AATranslation_overlap_Smik.fa',sep = '')

  bayFrameDNAFileName<-paste(path,'/',orfName, '_subalignment_Sbay.fa',sep = '')
  bayFrameAAFileName<-paste(path,'/',orfName, '_AATranslation_Sbay.fa',sep = '')
  bayFrameAAOverlapFileName<-paste(path,'/',orfName, '_AATranslation_overlap_Sbay.fa',sep = '')

  kudFrameDNAFileName<-paste(path,'/',orfName, '_subalignment_Skud.fa',sep = '')
  kudFrameAAFileName<-paste(path,'/',orfName, '_AATranslation_Skud.fa',sep = '')
  kudFrameAAOverlapFileName<-paste(path,'/',orfName, '_AATranslation_overlap_Skud.fa',sep = '')


  #these checks are for checking if there is a para(or mik,kud,bay) frame file exists
  #following fileChecks are only for larger para or mik frames
  #at this point smaller para-mik ORFs are not saved but are analyzed within this script
  paraCheck<-FALSE
  mikCheck<-FALSE
  bayCheck<-FALSE
  kudCheck<-FALSE
  if(file.exists(paraFrameAAFileName)){
    paraFrameAA<-readAAStringSet(paraFrameAAFileName)
    overlappingPara <- readAAStringSet(paraFrameAAOverlapFileName)
    paraCheck<-TRUE
  }
  if(file.exists(mikFrameAAFileName)){
    mikFrameAA<-readAAStringSet(mikFrameAAFileName)
    overlappingMik <- readAAStringSet(mikFrameAAOverlapFileName)

    mikCheck<-TRUE
  }

  if(file.exists(kudFrameAAFileName)){
    kudFrameAA<-readAAStringSet(kudFrameAAFileName)
    overlappingKud <- readAAStringSet(kudFrameAAOverlapFileName)

    kudCheck<-TRUE
  }

  if(file.exists(bayFrameAAFileName)){
    bayFrameAA<-readAAStringSet(bayFrameAAFileName)
    overlappingBay <- readAAStringSet(bayFrameAAOverlapFileName)
    bayCheck<-TRUE
  }
  i=1
  dataTable$`ORF Name`[i]<-orfName
  #If your alignment files do not have proto-gene sequence as last sequence, please append an arbitrary sequence to end, almost all functions assumes the last sequence is the proto-gene, thus makes calculations accordingly
  if(aa[1]!=aa[length(aa)]){#this number should be changed according to (number of aligned species+1) if above is the case
    aa<-append(aa,aa[1])
  }
  if(subalign[1]!=subalign[length(subalign)]){#this number should be changed according to (number of aligned species+1) if above is the case
    subalign<-append(subalign,subalign[1])
  }
  # if(aa_macse[1]!=aa_macse[length(aa_macse)]){#this number should be changed according to (number of aligned species+1) if above is the case
  #   aa_macse<-append(aa_macse,aa_macse[1])
  # }
  # aa<-append(aa,aa[1])

  startDNA<-startCodon(subalign)
  startAA<-startCodon(aa)
  stopDNA<-stopCodon(subalign)
  stopAA<-stopCodon(aa)
  # startMacseAA<-startCodon(aa_macse)
  # stopMacseAA<-stopCodon(aa_macse)
  dataTable$`Do all Start Codon Align?`[i]<-checkAlignment(startAA)
  dataTable$`Do all Stop Codon Align`[i]<-checkAlignment(stopAA)
  dataTable$`Length of ORF`[i]<-nchar(turnWoGaps(subalign[[1]]))
  dataTable$`Length of Amino Acid Sequence ORF`[i]<-nchar(turnWoGaps(aa[[1]]))


  #Para analysis####
  if(vec[2]==1){
    #paraFrameShifts<-countFrameShifting(aa_macse.p[2])

    if(paraCheck){
      dataTable$`Is there an ORF in the Para Amino Acid Sequence?`[i]<-TRUE
      #overlappingPara<-#findSmallFrame(aa,2,paraFrameAA[[2]])
      # paraFrameAA <-# findSmallFrame(paraFrameAA,2)

    }else{
      dataTable$`Is there an ORF in the Para Amino Acid Sequence?`[i]<-FALSE
      paraFrameAA<-FALSE
    }
    dataTable$`What is the Start Codon in Para?`[i]<-as.character(startAA[[2]])
    dataTable$`What is the Start Codon(DNA) in Para?`[i]<-as.character(startDNA[[2]])
    #dataTable$`What is the Start Codon in Para with Macse?`[i]<-as.character(startMacseAA[[2]])
    #dataTable$`What is the Stop Codon in Para with Macse?`[i]<-as.character(stopMacseAA[[2]])
    dataTable$`Does start codon align in Para?`[i]<-(as.character(startAA[[2]])==as.character(startAA[[1]]))
    dataTable$`What is the Stop Codon in Para?`[i]<-as.character(stopAA[[2]])
    dataTable$`Para % DNA ID over Smorf Frame`[i]<-calcIdentity(subalign)[2]
    dataTable$`Length of Para DNA Sequence over Smorf Frame`[i]<-nchar(turnWoGaps(subalign[[2]]))
    dataTable$`Number of Para Gaps over Smorf`[i]<-length(subalign[[2]])-dataTable$`Length of Para DNA Sequence over Smorf Frame`[i]

    if(is.logical(paraFrameAA)!=TRUE){
      dataTable$`Length of Para Amino Acid start to finish with Gaps`[i]<-nchar(paraFrameAA[[2]])
      dataTable$`Length of Para Amino Acid Start to Finish without Gaps`[i]<-nchar(turnWoGaps(paraFrameAA[[2]]))
      dataTable$`Para Length Ratio`[i] <- nchar(turnWoGaps(paraFrameAA[[2]]))/dataTable$`Length of Amino Acid Sequence ORF`[i]
      dataTable$`Para % Amino Acid over Para Overlap`[i]<-calcIdentity(overlappingPara)[2]
      dataTable$`Para Number of identical Amino Acid over Para Overlap`[i] <- calcIdentity(overlappingPara,percent = F)[2]
    }

    dataTable$`Para % Amino Acid over Smorf Frame`[i]<-calcIdentity(aa)[2]
    # dataTable$`Para % Amino Acid over smorf with Macse`[i]<-calcIdentity(aa_macse)[2]
    # macseParaAA<-findSmallFrame(aa_macse,2)
    # dataTable$`Para % Amino Acid over Para with Macse`[i]<-ifelse(is.logical(macseParaAA) , NA, calcIdentity(macseParaAA)[2])
    dataTable$`Where are other "M" in Para over the Smorf`[i]<-findM(aa[[2]])
    dataTable$`Where are other "X" in the Para over the Smorf`[i]<- findX(aa[[2]])
    dataTable$`Is length of ORF same in para?`[i]<-dataTable$`Length of Amino Acid Sequence ORF`[i]==dataTable$`Length of Para Amino Acid Start to Finish without Gaps`[i]
    #dataTable$`Para Number of FrameShifts`[i]<-paraFrameShifts
  }
  #Mik analysis####

  if(vec[3]==1){
    mikID <- which(types=='mik')
    # mikFrameShifts<-countFrameShifting(aa_macse.p[mikID])

    if(mikCheck){
      dataTable$`Is there an ORF in the Mik Amino Acid Sequence?`[i]<-TRUE
      #overlappingMik<-findSmallFrame(aa,mikID,mikFrameAA[[mikID]])
      # mikFrameAA <- findSmallFrame(mikFrameAA,mikID)

    }else{

      dataTable$`Is there an ORF in the Mik Amino Acid Sequence?`[i]<-FALSE
      mikFrameAA<-FALSE
    }

    dataTable$`What is the Start Codon in Mik?`[i]<-as.character(startAA[[mikID]])
    dataTable$`What is the Start Codon(DNA) in Mik?`[i]<-as.character(startDNA[[mikID]])
    #dataTable$`What is the Start Codon in Mik with Macse?`[i]<-as.character(startMacseAA[[mikID]])
    #dataTable$`What is the Stop Codon in Mik with Macse?`[i]<-as.character(stopMacseAA[[mikID]])
    dataTable$`Does start codon align in Mik?`[i]<-(as.character(startAA[[mikID]])==as.character(startAA[[1]]))
    dataTable$`What is the Stop Codon in Mik?`[i]<-as.character(stopAA[[mikID]])
    dataTable$`Mik % DNA ID over Smorf Frame`[i]<-calcIdentity(subalign)[mikID]
    dataTable$`Length of Mik DNA Sequence over Smorf Frame`[i]<-nchar(turnWoGaps(subalign[[mikID]]))
    dataTable$`Number of Mik Gaps over Smorf`[i]<-length(subalign[[mikID]])-dataTable$`Length of Mik DNA Sequence over Smorf Frame`[i]

    if(is.logical(mikFrameAA)!=TRUE){
      dataTable$`Length of Mik Amino Acid start to finish with Gaps`[i]<-nchar(mikFrameAA[[2]])
      dataTable$`Length of Mik Amino Acid Start to Finish without Gaps`[i]<-nchar(turnWoGaps(mikFrameAA[[2]]))
      dataTable$`Mik Length Ratio`[i] <- nchar(turnWoGaps(mikFrameAA[[2]]))/dataTable$`Length of Amino Acid Sequence ORF`[i]
      dataTable$`Mik % Amino Acid over Mik Overlap`[i]<-calcIdentity(overlappingMik)[2]
      dataTable$`Mik Number of identical Amino Acid over Mik Overlap`[i] <- calcIdentity(overlappingMik,percent = F)[2]

    }
    dataTable$`Mik % Amino Acid over Smorf Frame`[i]<-calcIdentity(aa)[mikID]
    # dataTable$`Mik % Amino Acid over smorf with Macse`[i]<-calcIdentity(aa_macse)[mikID]
    # macseMikAA<-findSmallFrame(aa_macse,mikID)
    # dataTable$`Mik % Amino Acid over Mik with Macse`[i]<-ifelse(is.logical(macseMikAA) , NA, calcIdentity(macseMikAA)[2])
    dataTable$`Where are other "M" in Mik over the Smorf`[i]<-findM(aa[[mikID]])
    dataTable$`Where are other "X" in the Mik over the Smorf`[i]<- findX(aa[[mikID]])
    #dataTable$`Mik Number of FrameShifts`[i]<-mikFrameShifts

    dataTable$`Is length of ORF same in mik?`[i]<-dataTable$`Length of Amino Acid Sequence ORF`[i]==dataTable$`Length of Mik Amino Acid Start to Finish without Gaps`[i]
  }
  #Kud analysis####
  if(vec[4]==1){
    kudID <- which(types=='kud')

    # kudFrameShifts<-countFrameShifting(aa_macse.p[kudID])
    if(kudCheck){
      dataTable$`Is there an ORF in the Kud Amino Acid Sequence?`[i]<-TRUE
      # overlappingKud<-findSmallFrame(aa,kudID, kudFrameAA[[kudID]])
      # kudFrameAA <- findSmallFrame(kudFrameAA,kudID)

    }else{
      dataTable$`Is there an ORF in the Kud Amino Acid Sequence?`[i]<-FALSE
      kudFrameAA<-FALSE
    }
    dataTable$`What is the Start Codon in Kud?`[i]<-as.character(startAA[[kudID]])
    dataTable$`What is the Start Codon(DNA) in Kud?`[i]<-as.character(startDNA[[kudID]])
    # dataTable$`What is the Start Codon in Kud with Macse?`[i]<-as.character(startMacseAA[[kudID]])
    # dataTable$`What is the Stop Codon in Kud with Macse?`[i]<-as.character(stopMacseAA[[kudID]])
    dataTable$`Does start codon align in Kud?`[i]<-(as.character(startAA[[kudID]])==as.character(startAA[[1]]))
    dataTable$`What is the Stop Codon in Kud?`[i]<-as.character(stopAA[[kudID]])
    dataTable$`Kud % DNA ID over Smorf Frame`[i]<-calcIdentity(subalign)[kudID]
    dataTable$`Length of Kud DNA Sequence over Smorf Frame`[i]<-nchar(turnWoGaps(subalign[[kudID]]))
    dataTable$`Number of Kud Gaps over Smorf`[i]<-length(subalign[[kudID]])-dataTable$`Length of Kud DNA Sequence over Smorf Frame`[i]

    if(is.logical(kudFrameAA)!=TRUE){
      dataTable$`Length of Kud Amino Acid start to finish with Gaps`[i]<-nchar(kudFrameAA[[2]])
      dataTable$`Length of Kud Amino Acid Start to Finish without Gaps`[i]<-nchar(turnWoGaps(kudFrameAA[[2]]))
      dataTable$`Kud Length Ratio`[i] <- nchar(turnWoGaps(kudFrameAA[[2]]))/dataTable$`Length of Amino Acid Sequence ORF`[i]
      dataTable$`Kud % Amino Acid over Kud Overlap`[i]<-calcIdentity(overlappingKud)[2]
      dataTable$`Kud Number of identical Amino Acid over Kud Overlap`[i] <- calcIdentity(overlappingKud,percent = F)[2]

    }
    dataTable$`Kud % Amino Acid over Smorf Frame`[i]<-calcIdentity(aa)[kudID]
    # dataTable$`Kud % Amino Acid over smorf with Macse`[i]<-calcIdentity(aa_macse)[kudID]
    # macseKudAA<-findSmallFrame(aa_macse,kudID)
    # dataTable$`Kud % Amino Acid over Kud with Macse`[i]<-ifelse(is.logical(macseKudAA) , NA, calcIdentity(macseKudAA)[2])
    dataTable$`Where are other "M" in Kud over the Smorf`[i]<-findM(aa[[kudID]])
    dataTable$`Where are other "X" in the Kud over the Smorf`[i]<- findX(aa[[kudID]])
    #dataTable$`Kud Number of FrameShifts`[i]<-kudFrameShifts

    dataTable$`Is length of ORF same in Kud?`[i]<-dataTable$`Length of Amino Acid Sequence ORF`[i]==dataTable$`Length of Kud Amino Acid Start to Finish without Gaps`[i]
  }

  #Bay analysis####

  if(vec[5]==1){
    bayID <- which(types=='bay')

    #bayFrameShifts<-countFrameShifting(aa_macse.p[bayID])
    if(bayCheck){
      dataTable$`Is there an ORF in the Bay Amino Acid Sequence?`[i]<-TRUE
      # overlappingBay<-findSmallFrame(aa,bayID,bayFrameAA[[bayID]])
      #bayFrameAA <- findSmallFrame(bayFrameAA,bayID)

    }else{
      dataTable$`Is there an ORF in the Bay Amino Acid Sequence?`[i]<-FALSE
      bayFrameAA<-FALSE

    }
    dataTable$`What is the Start Codon in Bay?`[i]<-as.character(startAA[[bayID]])
    dataTable$`What is the Start Codon(DNA) in Bay?`[i]<-as.character(startDNA[[bayID]])
    # dataTable$`What is the Start Codon in Bay with Macse?`[i]<-as.character(startMacseAA[[bayID]])
    # dataTable$`What is the Stop Codon in Bay with Macse?`[i]<-as.character(stopMacseAA[[bayID]])
    dataTable$`Does start codon align in Bay?`[i]<-(as.character(startAA[[bayID]])==as.character(startAA[[1]]))
    dataTable$`What is the Stop Codon in Bay?`[i]<-as.character(stopAA[[bayID]])
    dataTable$`Bay % DNA ID over Smorf Frame`[i]<-calcIdentity(subalign)[bayID]
    dataTable$`Length of Bay DNA Sequence over Smorf Frame`[i]<-nchar(turnWoGaps(subalign[[bayID]]))
    dataTable$`Number of Bay Gaps over Smorf`[i]<-length(subalign[[bayID]])-dataTable$`Length of Bay DNA Sequence over Smorf Frame`[i]

    if(is.logical(bayFrameAA)!=TRUE){
      dataTable$`Length of Bay Amino Acid start to finish with Gaps`[i]<-nchar(bayFrameAA[[2]])
      dataTable$`Length of Bay Amino Acid Start to Finish without Gaps`[i]<-nchar(turnWoGaps(bayFrameAA[[2]]))
      dataTable$`Bay Length Ratio`[i] <- nchar(turnWoGaps(bayFrameAA[[2]]))/dataTable$`Length of Amino Acid Sequence ORF`[i]
      dataTable$`Bay % Amino Acid over Bay Overlap`[i]<-calcIdentity(overlappingBay)[2]
      dataTable$`Bay Number of identical Amino Acid over Bay Overlap`[i] <- calcIdentity(overlappingBay,percent = F)[2]


    }
    dataTable$`Bay % Amino Acid over Smorf Frame`[i]<-calcIdentity(aa)[bayID]
    # dataTable$`Bay % Amino Acid over smorf with Macse`[i]<-calcIdentity(aa_macse)[bayID]
    # macseBayAA<-findSmallFrame(aa_macse,bayID)
    # dataTable$`Bay % Amino Acid over Bay with Macse`[i]<-ifelse(is.logical(macseBayAA) , NA, calcIdentity(macseBayAA)[2])
    dataTable$`Where are other "M" in Bay over the Smorf`[i]<-findM(aa[[bayID]])
    dataTable$`Where are other "X" in the Bay over the Smorf`[i]<- findX(aa[[bayID]])
    # dataTable$`Bay Number of FrameShifts`[i]<-bayFrameShifts

    dataTable$`Is length of ORF same in Bay?`[i]<-dataTable$`Length of Amino Acid Sequence ORF`[i]==dataTable$`Length of Bay Amino Acid Start to Finish without Gaps`[i]

  }
  return(dataTable)
}
