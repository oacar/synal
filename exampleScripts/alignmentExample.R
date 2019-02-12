library(muscle)
library(msa)
source('~/Main/synal/R/readFiles.R')
#source('functions.R')
library(synal)
options(warn = 1)
#to remove all alternative start codons. Bc I wanted to only ATG to be translated as M, the defaults were different, you can read ?GENETIC_CODE
attr(GENETIC_CODE, "alt_init_codons")=character(0)

p_ <- 'gene/genes_alignment_1/'
gene.rand <- read_csv('rand.nones.csv')
gene.rand <- gene.rand$orf_name

for ( i in 1:length(gene.rand)){
  dirName <- gene.rand[i]
  if(i==35 | i==75 | i==34 | i==19){
    next
  }
  path<-'data/alignmentInput/sequenceFilescopy/YAL016C-A/'#paste('', sep="/")
  print(dirName)
  #read files####
  if(dir.exists(path)==F){
    warning(paste(path,'does not exist'))
    next
  }
  mySequences<-readFiles(path)


  #current script is to align 5 sequences and proto-gene sequence being the 6th. to get the subalignment you need the proto-gene sequence

  #alignment####
  dnaAlignmentList<-synal::alignWoSmorf(mySequences,'ClustalW')

  DNAStr<-dnaAlignmentList[[1]]
  start<-dnaAlignmentList[[2]]
  stop<-dnaAlignmentList[[3]]
  if(start<0 || stop <0){
    warning(paste('No smorf sequence within these gene range for  ',dirName , sep = ''))
    if(TRUE){
      nPath<-paste(noseqrange,dirName, sep="/") #If regex cannot find the smorf sequence #1-Seq might be changed or #2-input files are wrong
      file.rename(path,nPath)
      if(TRUE){next}
    }
  }
  #
  if(start<=50){
    sStart<-1
  }else{
    sStart<-start-50
  }
  if(abs(stop-length(DNAStr[[1]]))<=50){
    sStop<-length(DNAStr[[1]])
  }else{
    sStop<-stop+50
  }
  smallerDNAStr<-subalignment(sStart,sStop, DNAStr)[[1]]


  msaPrettyPrint(msa(smallerDNAStr))
  seqs <- DNAStringSet()
  for(j in 1:length(smallerDNAStr)){
    seqs <- append(seqs,DNAStringSet(turnWoGaps(smallerDNAStr[[j]])))
  }
  names(seqs) <- names(smallerDNAStr)
  seq_arb <- arb[grepl(name,names(arb),perl=T)]
  seqs <- append(seqs,seq_arb)
  names(seqs)[length(seqs)] <- 'Sarb'
  writeXStringSet(seqs,paste('gene/gene_100/',dirName,'.fa',sep=''))

}
