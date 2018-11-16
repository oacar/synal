#source('alignmentFun.R')
library(Biostrings)
library(BiocGenerics)
library(IRanges)
library(S4Vectors)
source('InputHelper.R')
library(readxl)
options(warn = 1)

msg<-file('msg.Rout',open='wt')
out<-file('out.Rout',open='wt')
sink(msg,type='message')
sink(out,type='output')

#Functions to find nearest genes and y-gene seq----
findClosestGenes<- function(beg, stop, chr, path, name,chromosomeVec,startVec,stopVec){
  ##INPUT: beg: start coordinate of smorf, end: end coordinate of orf, chr: chromosome of orf (should be character as shown on SDG ['I','II' ...]))
  ## file path for scer coding orfs
  ##Output: DNAStringSet object containing nearest gene information
  if(file.exists(path)!=TRUE){
    stop('File does not exists. Give correct file name and path')
  }
  if((round(beg)!=beg || round(end)!=end)){
    stop('Beginning and End should be integers!')
  }
  else if(beg<0 || end<0){
    stop('Beginning and end should be positive')
  }
  scer<-readDNAStringSet(path)
  #dummy<-c()
  #coordinates<-c()
  #startVec<-c()
  #stopVec<-c()
  #chromosomeVec<-c()
  #namesVec<-c()
  nearestList<-c()
  minList<-c()
  ##startVec has start coordinates of genes
  ##stopVec has end coordinates of genes
  min<-Inf
  minIndex<-0
  for(i in 1:length(startVec)){
    if(chromosomeVec[i]==chr){
      if(beg>end){
        u_<-beg
        beg<-end
        end<-u_
      }
      rs<-end-1000
      re<-beg+1000
      if((stopVec[i]>=rs && stopVec[i]<=re) || (startVec[i]<=re && rs<=startVec[i]) || (startVec[i]<=rs &&  re<=stopVec[i] )){
        #min<-min(begToStart,endToStart,begToStop,endToStop)
        minIndex<-i
        nearestList<-append(nearestList,minIndex)
        minList<-append(minList,min)
      }
      # begToStart<-abs(beg-startVec[i])
      # endToStart<-abs(end-startVec[i])
      # begToStop<-abs(beg-stopVec[i])
      # endToStop<-abs(end-stopVec[i])
      # if(max(begToStart,endToStart,begToStop,endToStop)<=1000){
      #   min<-min(begToStart,endToStart,begToStop,endToStop)
      #   minIndex<-i
      #   nearestList<-append(nearestList,minIndex)
      #   minList<-append(minList,min)
      #   
      #}
    }
  }
  ng<-scer[nearestList]
  ng_name<-sapply(strsplit(names(ng), ','), '[', 1)
  if((name %in% ng_name)==FALSE){
    ng_name<-c(ng_name, name)
  }
  return(list(ng, ng_name))
  
}
findYGeneSeq<- function(path, geneName, outputPath){
  ##Input:path : name and path for orfGenomic file
  ##      geneName : annotation name that starts with letter Y
  ##      outputPath : path for output file to be saved
  ##output : saves respective sequence file to outputPath
  
  if(file.exists(path)!=TRUE){
    stop('File does not exist. Give a correct path and file name')
  }  
  if(nchar(geneName)<7){
    stop('Gene Name should be equal or bigger than 7')
  }
  index<-0
  
  seqfile<-(readDNAStringSet(path))
  
  for (i in 1:length(seqfile)){
    if(geneName ==  strsplit(strsplit(names(seqfile)[i], ',')[[1]][1], ' ')[[1]][1]){
      index<-i
    }
  }
  if(index==0){
    
    warning(sprintf("%s", paste(geneName, 'cannot be found ', sep=" ")))
    
  }else{  Seq<-seqfile[index]
  
  if(is.null(outputPath)){
    return(Seq) 
  }else{
    fileName<-sprintf("%s", paste(outputPath,geneName, sep="/"))
    fileName<-sprintf("%s", paste(fileName,"sequence", sep="_"))
    fileName<-sprintf("%s", paste(fileName,'.fa', sep=""))
    writeXStringSet(Seq, fileName)
  }
  
  return(Seq) 
  }
}

#Create scer gene DataFrame----
pathGenomicOrf<-"../../alignment/scer/orf_coding.fasta"
genecol.names<-c('orf_name','gene_name','reverse-complement', 'coor1', 'coor2', 'chr','description')
scerGenes<-data.frame(matrix(nrow = 5917, ncol =7 ))
colnames(scerGenes)<-genecol.names


scer<-readDNAStringSet(pathGenomicOrf)
dummy<-c()
coordinates<-c()
startVec<-c()
stopVec<-c()
chromosomeVec<-c()
namesVec<-c()
nearestList<-c()
minList<-c()
#  descVec
for (i in 1:length(scer)){
  if(i==410 || i==1574 || i == 2042){
    namesVec<-c(namesVec, strsplit(strsplit(names(scer)[i],',')[[1]][1],' ')[[1]][1])
    dummy<-c(dummy,strsplit(names(scer)[i],',')[[1]][3])
    coordinates<-c(coordinates,strsplit(dummy[i], ' ')[[1]][5])
    chromosomeVec<-c(chromosomeVec,strsplit(dummy[i], ' ')[[1]][3])
    startVec<-c(startVec, strsplit(coordinates[i], '-')[[1]][1])
    stopVec<-c(stopVec, strsplit(coordinates[i], '-')[[1]][2])
  }
  else{
    namesVec<-c(namesVec, strsplit(strsplit(names(scer)[i],',')[[1]][1],' ')[[1]][1])
    dummy<-c(dummy,strsplit(names(scer)[i],',')[[1]][2])
    coordinates<-c(coordinates,strsplit(dummy[i], ' ')[[1]][5])
    chromosomeVec<-c(chromosomeVec,strsplit(dummy[i], ' ')[[1]][3])
    startVec<-c(startVec, strsplit(coordinates[i], '-')[[1]][1])
    stopVec<-c(stopVec, strsplit(coordinates[i], '-')[[1]][2])
    
  }}
startVec<-as.integer(startVec)
stopVec<-as.integer(stopVec)
for(i in 1:nrow(scerGenes)){
  scerGenes$orf_name[i]<-namesVec[i]
  scerGenes$coor1[i]<-startVec[i]
  scerGenes$coor2[i]<-stopVec[i]
  scerGenes$chr[i]<-chromosomeVec[i]
  scerGenes$gene_name[i]<-strsplit(strsplit(names(scer)[i],',')[[1]][1],' ')[[1]][2]
  t<-strsplit(names(scer)[i],'"')
  scerGenes$description[i]<-t[[1]][length(t[[1]])]
  if(grepl('reverse complement', names(scer)[i])==TRUE){
    scerGenes$`reverse-complement`[i]<-TRUE
  }else{
    scerGenes$`reverse-complement`[i]<-FALSE
    
  }
}


#Define +-1kb file paths----
p<-'../../alignment/scer/orf_genomic_all.fasta'
pathScer<-'../../alignment/scer/orf_genomic_1000_all.fasta'
pathSMik = "../../alignment/smik/orf_genomic_1000.fasta"
pathSPar = "../../alignment/spar/orf_genomic_1000.fasta"
pathSBay ='../../alignment/sbay/orf_genomic_1000.fasta'
pathSKud<-'../../alignment/skud/orf_genomic_1000.fasta'
spar<-readDNAStringSet(pathSPar)
smik<-readDNAStringSet(pathSMik)
sbay<-readDNAStringSet(pathSBay)
skud<-readDNAStringSet(pathSKud)
scer<-readDNAStringSet(pathScer)
#Define protogene list----
protogeneFile<-readxl::read_excel("~/Box Sync/Carvunis Dry Lab/Karlovich_Summer2017/Other Documents/6-6-17_protogenes_2012_annotated_Karlovich1.xlsx")
smorfs<-subset(protogeneFile, grepl('smorf',protogeneFile$orf_name))
smorfs<-smorfs[c(1,8,9,14)]
y_gene<-subset(protogeneFile, grepl('smorf',protogeneFile$orf_name)!=TRUE)
y_gene<-y_gene[c(1,8,9,14)]
dic<-c('I','II','III','IV', 'V', 'VI', 'VII', 'VIII', 'IX','X','XI','XII','XIII','XIV','XV','XVI')
names(dic) <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
y_gene<-adj
colnames(y_gene)<-c('orf_name','gene')
#find nearest genes----
names<-unique(y_gene$orf_name)
clGenes<-c()

for (i in 1:length(names)){
  
  print(i)

  smorfName<-names[i]
  ygeneSeq<-findYGeneSeq(p, smorfName, NULL)
  dummy<-strsplit(names(ygeneSeq),',')[[1]][2]
  start<-as.integer(strsplit(strsplit(dummy, ' ')[[1]][5], '-')[[1]][1])
  end<-as.integer(strsplit(strsplit(dummy, ' ')[[1]][5], '-')[[1]][2])
  chr<-strsplit(dummy, ' ')[[1]][3]
#  chr<-dic[[chrNum]]
  ng<-findClosestGenes(start, end, chr, pathGenomicOrf, smorfName,chromosomeVec, startVec,stopVec)
  clGenes<-c(clGenes,(paste(smorfName, ' -->', ng[[2]])))
}
print('clGenes created')
clGenesDT<-data.frame(matrix(ncol = 2, nrow = length(clGenes)),stringsAsFactors = FALSE)
for (i in 1:length(clGenes)){
  splt<-strsplit(clGenes[i],'  --> ')
  clGenesDT$X1[i]<-splt[[1]][1]
  clGenesDT$X2[i]<-strsplit(splt[[1]][2],' ')[[1]][1]
}
print('clGenesDT created')

#Find the nearest gene that has most homologous region with other species----
count<-0
for(i in 1:length(names)){
  #i<-5
  outputPath<-paste("y_genes-2",names[i], sep="/")
  #newOutputPath<-paste("./","y_genes_2",y_gene$orf_name[i], sep="/")
  #chrNum<- as.integer(smorfs$chromosome[i])
  chrNum<- as.integer(y_gene$chromosome[i])
  smorfName<-names[i]
  #ygeneSeq<-findYGeneSeq(p, smorfName, NULL)
  ngs<-clGenesDT[which(clGenesDT$X1==smorfName),2]
  #print(ngs)
  nums<-c()
  for(i_ in 1:length(ngs)){
    #scerseq<-
    mikcheck<-any(grepl(ngs[i_],names(smik),fixed = TRUE))
    parcheck<-any(grepl(ngs[i_],names(spar),fixed = TRUE))
    baycheck<-any(grepl(ngs[i_],names(sbay),fixed = TRUE))
    kudcheck<-any(grepl(ngs[i_],names(skud),fixed = TRUE))
    nums<-c(nums,sum(mikcheck,parcheck,baycheck,kudcheck))
  }
  if(any(is.na(nums))){
    warning(paste('No nearest gene for ',smorfName,sep = ''))
    if(TRUE){
      next
    }
  }else if(max(nums)==4){
    ng<-ngs[which.max(nums)]
    dir.create(outputPath)
    print(paste(smorfName,' ',ng,sep = ''))
    count<-count+1
  }else{
    warning(paste('Not all species has syntenic region defined for ',smorfName,sep = ''))
    if(TRUE){
      next
    }
  }
  nearestGene<-ng
  #nameNearestGene <-ng[[2]][1]
  #if(scerGenes[which(scerGenes$orf_name==nearestGene),3]!=TRUE && y_gene$strand[i]=='-'){
  if(scerGenes[which(scerGenes$orf_name==nearestGene),3]!=TRUE && y_gene[which(y_gene$orf_name==names[i]),3]=='-'){
      check<-TRUE
  #}else if(scerGenes[which(scerGenes$orf_name==nearestGene),3]==TRUE && y_gene$strand[i]=='+'){
  }else if(scerGenes[which(scerGenes$orf_name==nearestGene),3]==TRUE && y_gene[which(y_gene$orf_name==names[i]),3]=='+'){
    check<-TRUE
  }else{
    check<-FALSE
  }
  checkMik<-findSeq(pathSMik, nearestGene, outputPath,check,smorfName )
  checkPara<-findSeq(pathSPar,nearestGene,outputPath,check,smorfName)
  findSeq(pathSBay,nearestGene,outputPath,check,smorfName)
  findSeq(pathSKud,nearestGene,outputPath,check,smorfName)


  findSeq(pathScer,nearestGene,outputPath,check,smorfName)
  ygeneSeq<-findYGeneSeq(p, smorfName, outputPath)

  
}
