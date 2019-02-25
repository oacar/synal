library(tidyverse)
library(Biostrings)
library(BiocGenerics)
library(IRanges)
library(S4Vectors)
library(openxlsx)
library(gower)
library(cluster)
library(gplots)
library(dendextend)
library(dplyr)
library(tsne)
library(BBmisc)
library(StatMatch)
library(synal)
options(stringsAsFactors = FALSE)
#read input files####
overlapping <- read.csv('../../anne/network_analysis/overlappingorfs.csv')
pgs_nodes<-read.csv('../data/clusteringInput/filtering_csv/protogenes_nonoverlapping.csv')
ppp <- read.csv('~/Main/anne/proto_genes',header = F)
ppp <- ppp$V1
pccFiltered<-read.csv('../data/clusteringInput/filtering_csv/protogenes_pccfiltered.csv')
#dataTable<-read.xlsx('data/clusteringInput/protogenes.xlsx')
dataTable <- read.xlsx('2-7-19_296_youngorfs_inbarflex_omer.xlsx')
dt2 <- read.xlsx('../data.noorfxm.xlsx')
dataTable <- dataTable[dataTable$ORF.Name%in%ppp,]
dt2  <- dt2[dt2$ORF.Name%in%ppp,]
#dt2 <- dt2[(dt2$ORF.Name=='YER175W-A' | dt2$ORF.Name=='YAL037C-A')==F,]
#dataTable <- rbind(dataTable,dt2)
geneDataTable<-read.xlsx('../data/clusteringInput/genes.xlsx')
nonGeneDataTable <- read.xlsx('../data/clusteringInput/nongenes.xlsx')
geneDataTable <- sample_n(geneDataTable,100)
overlappingNongenes <- read.csv('~/Main/anne/overlappingNongenes.csv')
overlappingNongenes <- overlappingNongenes$x
nonGeneDataTable <- nonGeneDataTable[is.na(nonGeneDataTable$ORF.Name)==F & nonGeneDataTable$ORF.Name%in%overlappingNongenes==F & nonGeneDataTable$Length.of.ORF>=83,]
nonGeneDataTable <- sample_n(nonGeneDataTable,100)

colNames<-c("Orf Name","ORF AA Length", #" ORF DNA Length",
            "Is there overlapping ORF in Para", "start codon in Para","stop codon in Para",
            "Para AA% Identitiy over Smorf Frame", "Para AA% Identitiy over Para Frame", "Para length Ratio",
            "Macse start codon in Para",'Macse stop codon in Para','Macse Para Similarity','FrameShifting indels Para',"Para DNA similarity",'Macse Para overlap Identity',
            "Is there overlapping ORF in Mik", "start codon in Mik","stop codon in Mik",
            "Mik AA% Identitiy over Smorf Frame", "Mik AA% Identitiy over Mik Frame", "Mik length Ratio",
            "Macse start codon in Mik",'Macse stop codon in Mik','Macse Mik Similarity','FrameShifting indels Mik' ,"Mik DNA similarity",'Macse Mik overlap Identity',
            "Is there overlapping ORF in Bay", "start codon in Bay","stop codon in Bay",
            "Bay AA% Identitiy over Smorf Frame", "Bay AA% Identitiy over Bay Frame", "Bay length Ratio",
            "Macse start codon in Bay",'Macse stop codon in Bay','Macse Bay Similarity','FrameShifting indels Bay' ,"Bay DNA similarity",'Macse Bay overlap Identity',
            "Is there overlapping ORF in Kud", "start codon in Kud","stop codon in Kud",
            "Kud AA% Identitiy over Smorf Frame", "Kud AA% Identitiy over Kud Frame", "Kud length Ratio",
            "Macse start codon in Kud",'Macse stop codon in Kud','Macse Kud Similarity','FrameShifting indels Kud' ,"Kud DNA similarity",'Macse Kud overlap Identity', 'Is gene?')

clDT<-data.frame(matrix(ncol = length(colNames), nrow = (nrow(dataTable))))#+nrow(geneDataTable)+nrow(dt2)+nrow(nonGeneDataTable))))

colnames(clDT)<-colNames
i_<-1
dataTable <- dt2
dataTable <- nonGeneDataTable
for (i in 1:nrow(dataTable)){
  clDT$`Orf Name`[i_]<-dataTable$ORF.Name[i]
  clDT$`ORF AA Length`[i_]<-dataTable$Length.of.Amino.Acid.Sequence.ORF[i]
  clDT$`Is there overlapping ORF in Para`[i_]<-dataTable$`Is.there.an.ORF.in.the.Para.Amino.Acid.Sequence?`[i]
  clDT$`Para AA% Identitiy over Smorf Frame`[i_]<-dataTable$`Para.%.Amino.Acid.over.Smorf.Frame`[i]
  clDT$`Para AA% Identitiy over Para Frame`[i_]<-dataTable$`Para.%.Amino.Acid.over.Para.frame`[i]
  clDT$`Is there overlapping ORF in Mik`[i_]<-dataTable$`Is.there.an.ORF.in.the.Mik.Amino.Acid.Sequence?`[i]
  clDT$`Mik AA% Identitiy over Smorf Frame`[i_]<-dataTable$`Mik.%.Amino.Acid.over.Smorf.Frame`[i]
  clDT$`Mik AA% Identitiy over Mik Frame`[i_]<-dataTable$`Mik.%.Amino.Acid.over.Mik.frame`[i]
  clDT$`start codon in Para`[i_]<-dataTable$`What.is.the.Start.Codon.in.Para?`[i]
  clDT$`start codon in Mik`[i_]<-dataTable$`What.is.the.Start.Codon.in.Mik?`[i]
  clDT$`Para DNA similarity`[i_]<-dataTable$`Para.%.DNA.ID.over.Smorf.Frame`[i]
  clDT$`Mik DNA similarity`[i_]<-dataTable$`Mik.%.DNA.ID.over.Smorf.Frame`[i]
  clDT$`Para length Ratio`[i_]<-dataTable$`Length.of.Para.Amino.Acid.Start.to.Finish.without.Gaps`[i]/dataTable$Length.of.Amino.Acid.Sequence.ORF[i]
  clDT$`Mik length Ratio`[i_]<-dataTable$`Length.of.Mik.Amino.Acid.Start.to.Finish.without.Gaps`[i]/dataTable$Length.of.Amino.Acid.Sequence.ORF[i]

  clDT$`stop codon in Para`[i_]<-dataTable$`What.is.the.Stop.Codon.in.Para?`[i]
  clDT$`stop codon in Mik`[i_]<-dataTable$`What.is.the.Stop.Codon.in.Mik?`[i]
  clDT$`Macse start codon in Para`[i_]<-dataTable$`What.is.the.Start.Codon.in.Para.with.Macse?`[i]
  clDT$`Macse start codon in Mik`[i_]<-dataTable$`What.is.the.Start.Codon.in.Mik.with.Macse?`[i]
  clDT$`Macse stop codon in Para`[i_]<-dataTable$`What.is.the.Stop.Codon.in.Para.with.Macse?`[i]
  clDT$`Macse stop codon in Mik`[i_]<-dataTable$`What.is.the.Stop.Codon.in.Mik.with.Macse?`[i]
  clDT$`Macse Para Similarity`[i_]<-dataTable$`Para.%.Amino.Acid.over.smorf.with.Macse`[i]
  clDT$`Macse Mik Similarity`[i_]<-dataTable$`Mik.%.Amino.Acid.over.smorf.with.Macse`[i]
  clDT$`FrameShifting indels Para`[i_]<-dataTable$Para.Number.of.FrameShifts[i]
  clDT$`FrameShifting indels Mik`[i_]<-dataTable$Mik.Number.of.FrameShifts[i]
  clDT$`Macse Para overlap Identity`[i_]<-dataTable$`Para.%.Amino.Acid.over.Para.with.Macse`[i]
  clDT$`Macse Mik overlap Identity`[i_]<-dataTable$`Mik.%.Amino.Acid.over.Mik.with.Macse`[i]
  clDT$`Macse Bay overlap Identity`[i_]<-dataTable$`Bay.%.Amino.Acid.over.Bay.with.Macse`[i]
  clDT$`Macse Kud overlap Identity`[i_]<-dataTable$`Kud.%.Amino.Acid.over.Kud.with.Macse`[i]
  #
  clDT$`Is there overlapping ORF in Bay`[i_]<-dataTable$`Is.there.an.ORF.in.the.Bay.Amino.Acid.Sequence?`[i]
  clDT$`Bay AA% Identitiy over Smorf Frame`[i_]<-dataTable$`Bay.%.Amino.Acid.over.Smorf.Frame`[i]
  clDT$`Bay AA% Identitiy over Bay Frame`[i_]<-dataTable$`Bay.%.Amino.Acid.over.Bay.frame`[i]
  clDT$`Is there overlapping ORF in Kud`[i_]<-dataTable$`Is.there.an.ORF.in.the.Kud.Amino.Acid.Sequence?`[i]
  clDT$`Kud AA% Identitiy over Smorf Frame`[i_]<-dataTable$`Kud.%.Amino.Acid.over.Smorf.Frame`[i]
  clDT$`Kud AA% Identitiy over Kud Frame`[i_]<-dataTable$`Kud.%.Amino.Acid.over.Kud.frame`[i]
  clDT$`start codon in Bay`[i_]<-dataTable$`What.is.the.Start.Codon.in.Bay?`[i]
  clDT$`start codon in Kud`[i_]<-dataTable$`What.is.the.Start.Codon.in.Kud?`[i]
  clDT$`Bay DNA similarity`[i_]<-dataTable$`Bay.%.DNA.ID.over.Smorf.Frame`[i]
  clDT$`Kud DNA similarity`[i_]<-dataTable$`Kud.%.DNA.ID.over.Smorf.Frame`[i]
  clDT$`Bay length Ratio`[i_]<-dataTable$`Length.of.Bay.Amino.Acid.Start.to.Finish.without.Gaps`[i]/dataTable$Length.of.Amino.Acid.Sequence.ORF[i]
  clDT$`Kud length Ratio`[i_]<-dataTable$`Length.of.Kud.Amino.Acid.Start.to.Finish.without.Gaps`[i]/dataTable$Length.of.Amino.Acid.Sequence.ORF[i]
  clDT$`stop codon in Bay`[i_]<-dataTable$`What.is.the.Stop.Codon.in.Bay?`[i]
  clDT$`stop codon in Kud`[i_]<-dataTable$`What.is.the.Stop.Codon.in.Kud?`[i]
  clDT$`Macse start codon in Bay`[i_]<-dataTable$`What.is.the.Start.Codon.in.Bay.with.Macse?`[i]
  clDT$`Macse start codon in Kud`[i_]<-dataTable$`What.is.the.Start.Codon.in.Kud.with.Macse?`[i]
  clDT$`Macse stop codon in Bay`[i_]<-dataTable$`What.is.the.Stop.Codon.in.Bay.with.Macse?`[i]
  clDT$`Macse stop codon in Kud`[i_]<-dataTable$`What.is.the.Stop.Codon.in.Kud.with.Macse?`[i]
  clDT$`Is gene?`[i_]<-3
  clDT$`Macse Bay Similarity`[i_]<-dataTable$`Bay.%.Amino.Acid.over.smorf.with.Macse`[i]
  clDT$`Macse Kud Similarity`[i_]<-dataTable$`Kud.%.Amino.Acid.over.smorf.with.Macse`[i]
  clDT$`FrameShifting indels Bay`[i_]<-dataTable$Bay.Number.of.FrameShifts[i]
  clDT$`FrameShifting indels Kud`[i_]<-dataTable$Kud.Number.of.FrameShifts[i]

  i_<-i_+1

}

for (i in 1:nrow(geneDataTable)){
  clDT$`Orf Name`[i_]<-geneDataTable$ORF.Name[i]
  clDT$`ORF AA Length`[i_]<-geneDataTable$Length.of.Amino.Acid.Sequence.ORF[i]
  clDT$`Is there overlapping ORF in Para`[i_]<-geneDataTable$`Is.there.an.ORF.in.the.Para.Amino.Acid.Sequence?`[i]
  clDT$`Para AA% Identitiy over Smorf Frame`[i_]<-geneDataTable$`Para.%.Amino.Acid.over.Smorf.Frame`[i]
  clDT$`Para AA% Identitiy over Para Frame`[i_]<-geneDataTable$`Para.%.Amino.Acid.over.Para.frame`[i]
  clDT$`Is there overlapping ORF in Mik`[i_]<-geneDataTable$`Is.there.an.ORF.in.the.Mik.Amino.Acid.Sequence?`[i]
  clDT$`Mik AA% Identitiy over Smorf Frame`[i_]<-geneDataTable$`Mik.%.Amino.Acid.over.Smorf.Frame`[i]
  clDT$`Mik AA% Identitiy over Mik Frame`[i_]<-geneDataTable$`Mik.%.Amino.Acid.over.Mik.frame`[i]
  clDT$`start codon in Para`[i_]<-geneDataTable$`What.is.the.Start.Codon.in.Para?`[i]
  clDT$`start codon in Mik`[i_]<-geneDataTable$`What.is.the.Start.Codon.in.Mik?`[i]
  clDT$`Para DNA similarity`[i_]<-geneDataTable$`Para.%.DNA.ID.over.Smorf.Frame`[i]
  clDT$`Mik DNA similarity`[i_]<-geneDataTable$`Mik.%.DNA.ID.over.Smorf.Frame`[i]
  clDT$`Para length Ratio`[i_]<-geneDataTable$`Length.of.Para.Amino.Acid.Start.to.Finish.without.Gaps`[i]/geneDataTable$Length.of.Amino.Acid.Sequence.ORF[i]
  clDT$`Mik length Ratio`[i_]<-geneDataTable$`Length.of.Mik.Amino.Acid.Start.to.Finish.without.Gaps`[i]/geneDataTable$Length.of.Amino.Acid.Sequence.ORF[i]
  clDT$`stop codon in Para`[i_]<-geneDataTable$`What.is.the.Stop.Codon.in.Para?`[i]
  clDT$`stop codon in Mik`[i_]<-geneDataTable$`What.is.the.Stop.Codon.in.Mik?`[i]
  clDT$`Macse start codon in Para`[i_]<-geneDataTable$`What.is.the.Start.Codon.in.Para.with.Macse?`[i]
  clDT$`Macse start codon in Mik`[i_]<-geneDataTable$`What.is.the.Start.Codon.in.Mik.with.Macse?`[i]
  clDT$`Macse stop codon in Para`[i_]<-geneDataTable$`What.is.the.Stop.Codon.in.Para.with.Macse?`[i]
  clDT$`Macse stop codon in Mik`[i_]<-geneDataTable$`What.is.the.Stop.Codon.in.Mik.with.Macse?`[i]
  clDT$`Macse Para Similarity`[i_]<-geneDataTable$`Para.%.Amino.Acid.over.smorf.with.Macse`[i]
  clDT$`Macse Mik Similarity`[i_]<-geneDataTable$`Mik.%.Amino.Acid.over.smorf.with.Macse`[i]
  clDT$`FrameShifting indels Para`[i_]<-geneDataTable$Para.Number.of.FrameShifts[i]
  clDT$`FrameShifting indels Mik`[i_]<-geneDataTable$Mik.Number.of.FrameShifts[i]
  clDT$`Macse Para overlap Identity`[i_]<-geneDataTable$`Para.%.Amino.Acid.over.Para.with.Macse`[i]
  clDT$`Macse Mik overlap Identity`[i_]<-geneDataTable$`Mik.%.Amino.Acid.over.Mik.with.Macse`[i]
  clDT$`Macse Bay overlap Identity`[i_]<-geneDataTable$`Bay.%.Amino.Acid.over.Bay.with.Macse`[i]
  clDT$`Macse Kud overlap Identity`[i_]<-geneDataTable$`Kud.%.Amino.Acid.over.Kud.with.Macse`[i]

  clDT$`Is there overlapping ORF in Bay`[i_]<-geneDataTable$`Is.there.an.ORF.in.the.Bay.Amino.Acid.Sequence?`[i]
  clDT$`Bay AA% Identitiy over Smorf Frame`[i_]<-geneDataTable$`Bay.%.Amino.Acid.over.Smorf.Frame`[i]
  clDT$`Bay AA% Identitiy over Bay Frame`[i_]<-geneDataTable$`Bay.%.Amino.Acid.over.Bay.frame`[i]
  clDT$`Is there overlapping ORF in Kud`[i_]<-geneDataTable$`Is.there.an.ORF.in.the.Kud.Amino.Acid.Sequence?`[i]
  clDT$`Kud AA% Identitiy over Smorf Frame`[i_]<-geneDataTable$`Kud.%.Amino.Acid.over.Smorf.Frame`[i]
  clDT$`Kud AA% Identitiy over Kud Frame`[i_]<-geneDataTable$`Kud.%.Amino.Acid.over.Kud.frame`[i]
  clDT$`start codon in Bay`[i_]<-geneDataTable$`What.is.the.Start.Codon.in.Bay?`[i]
  clDT$`start codon in Kud`[i_]<-geneDataTable$`What.is.the.Start.Codon.in.Kud?`[i]
  clDT$`Bay DNA similarity`[i_]<-geneDataTable$`Bay.%.DNA.ID.over.Smorf.Frame`[i]
  clDT$`Kud DNA similarity`[i_]<-geneDataTable$`Kud.%.DNA.ID.over.Smorf.Frame`[i]
  clDT$`Bay length Ratio`[i_]<-geneDataTable$`Length.of.Bay.Amino.Acid.Start.to.Finish.without.Gaps`[i]/geneDataTable$Length.of.Amino.Acid.Sequence.ORF[i]
  clDT$`Kud length Ratio`[i_]<-geneDataTable$`Length.of.Kud.Amino.Acid.Start.to.Finish.without.Gaps`[i]/geneDataTable$Length.of.Amino.Acid.Sequence.ORF[i]
  clDT$`stop codon in Bay`[i_]<-geneDataTable$`What.is.the.Stop.Codon.in.Bay?`[i]
  clDT$`stop codon in Kud`[i_]<-geneDataTable$`What.is.the.Stop.Codon.in.Kud?`[i]
  clDT$`Macse start codon in Bay`[i_]<-geneDataTable$`What.is.the.Start.Codon.in.Bay.with.Macse?`[i]
  clDT$`Macse start codon in Kud`[i_]<-geneDataTable$`What.is.the.Start.Codon.in.Kud.with.Macse?`[i]
  clDT$`Macse stop codon in Bay`[i_]<-geneDataTable$`What.is.the.Stop.Codon.in.Bay.with.Macse?`[i]
  clDT$`Macse stop codon in Kud`[i_]<-geneDataTable$`What.is.the.Stop.Codon.in.Kud.with.Macse?`[i]
  clDT$`Is gene?`[i_]<-2
  clDT$`Macse Bay Similarity`[i_]<-geneDataTable$`Bay.%.Amino.Acid.over.smorf.with.Macse`[i]
  clDT$`Macse Kud Similarity`[i_]<-geneDataTable$`Kud.%.Amino.Acid.over.smorf.with.Macse`[i]
  clDT$`FrameShifting indels Bay`[i_]<-geneDataTable$Bay.Number.of.FrameShifts[i]
  clDT$`FrameShifting indels Kud`[i_]<-geneDataTable$Kud.Number.of.FrameShifts[i]
  i_<-i_+1
}


#create Data matrix####
P = clDT
paraLength<-P$`ORF AA Length`*P$`Para length Ratio`
mikLength<-P$`ORF AA Length`*P$`Mik length Ratio`
bayLength<-P$`ORF AA Length`*P$`Bay length Ratio`
kudLength<-P$`ORF AA Length`*P$`Kud length Ratio`


#normalizations####
#chose one of the below for length ratio normalization
#assign whatever >1 to 1
P$`Bay length Ratio`[P$`Bay length Ratio`>1]<-1
P$`Para length Ratio`[P$`Para length Ratio`>1]<-1
P$`Mik length Ratio`[P$`Mik length Ratio`>1]<-1
P$`Kud length Ratio`[P$`Kud length Ratio`>1]<-1

#take log10 and normalize between [0,1]
# P$`Mik length Ratio`<-normalize(log10(P$`Mik length Ratio`),method = 'range',range = c(0,1))
# P$`Para length Ratio`<-normalize(log10(P$`Para length Ratio`),method = 'range',range = c(0,1))
# P$`Kud length Ratio`<-normalize(log10(P$`Kud length Ratio`),method = 'range',range = c(0,1))
# P$`Bay length Ratio`<-normalize(log10(P$`Bay length Ratio`),method = 'range',range = c(0,1))
# P$`ORF AA Length`<-normalize(log10(P$`ORF AA Length`),method = 'range',range = c(0,1))

#divide number of frameshifts to length of the orf, when there is no orf it will be treated as NA
P$`FrameShifting indels Para`<-normalize((P$`FrameShifting indels Para`/paraLength),method = 'range',range = c(0,1))
P$`FrameShifting indels Mik`<-normalize((P$`FrameShifting indels Mik`/mikLength),method = 'range',range = c(0,1))
P$`FrameShifting indels Bay`<-normalize((P$`FrameShifting indels Bay`/bayLength),method = 'range',range = c(0,1))
P$`FrameShifting indels Kud`<-normalize((P$`FrameShifting indels Kud`/kudLength),method = 'range',range = c(0,1))
P <-P[(P$`Orf Name`=='YHR125W' | P$`Orf Name`=='YGL052' | P$`Orf Name`=='YBR131W')==F,]
P <-P[P$`Orf Name`%in%t==F,]
P[P=="N/A"]<-0 #make "N/A" cells =0
P[is.na(P)]<-0 #make value error cells = 0
#if dna similarity is less than 50%, make everything related with that species 0
P[which(P$`Mik DNA similarity`<0.5), (grep('Mik',colnames(P),ignore.case = TRUE))]<-0
P[which(P$`Para DNA similarity`<0.5), (grep('Para',colnames(P),ignore.case = TRUE))]<-0
P[which(P$`Kud DNA similarity`<0.5), (grep('Kud',colnames(P),ignore.case = TRUE))]<-0
P[which(P$`Bay DNA similarity`<0.5), (grep('Bay',colnames(P),ignore.case = TRUE))]<-0

#change column names####
pcolNames<-c("Orf Name","ORF AA Length", #" ORF DNA Length",
             "Is there overlapping ORF in Para", "start codon in Para","stop codon in Para",
             "Para AA% Identitiy over Smorf Frame", "Para overlap AA% Identitiy", "(Para length Ratio)",
             "Macse start codon in Para",'Macse stop codon in Para','Macse Para Similarity','FrameShifting indels Para',"Para DNA similarity",'Macse Para overlap Identity',
             "Is there overlapping ORF in Mik", "start codon in Mik","stop codon in Mik",
             "Mik AA% Identitiy over Smorf Frame", "Mik overlap AA% Identitiy", "(Mik length Ratio)",
             "Macse start codon in Mik",'Macse stop codon in Mik','Macse Mik Similarity','FrameShifting indels Mik' ,"Mik DNA similarity", 'Macse Mik overlap Identity',
             "Is there overlapping ORF in Bay", "start codon in Bay","stop codon in Bay",
             "Bay AA% Identitiy over Smorf Frame", "Bay overlap AA% Identitiy", "(Bay length Ratio)",
             "Macse start codon in Bay",'Macse stop codon in Bay','Macse Bay Similarity','FrameShifting indels Bay' ,"Bay DNA similarity", 'Macse Bay overlap Identity',
             "Is there overlapping ORF in Kud", "start codon in Kud","stop codon in Kud",
             "Kud AA% Identitiy over Smorf Frame", "Kud overlap AA% Identitiy", "(Kud length Ratio)",
             "Macse start codon in Kud",'Macse stop codon in Kud','Macse Kud Similarity','FrameShifting indels Kud' ,"Kud DNA similarity", 'Macse Kud overlap Identity',
             'Is gene?')
colnames(P)<-pcolNames

col_order <- c(#"ORF AA Length", #" ORF DNA Length",
  "Is there overlapping ORF in Para", "Is there overlapping ORF in Mik","Is there overlapping ORF in Kud", "Is there overlapping ORF in Bay",
  "start codon in Para","start codon in Mik","start codon in Kud","start codon in Bay","stop codon in Para","stop codon in Mik" ,"stop codon in Kud","stop codon in Bay",
  #"Macse start codon in Para","Macse start codon in Mik","Macse start codon in Kud", "Macse start codon in Bay",
  #'Macse stop codon in Para',
  #'Macse stop codon in Mik',
  #'Macse stop codon in Kud','Macse stop codon in Bay',



  #"Para AA% Identitiy over Smorf Frame",
  "Para overlap AA% Identitiy", #'Macse Para Similarity','Macse Para overlap Identity',
  #"Mik AA% Identitiy over Smorf Frame",
  "Mik overlap AA% Identitiy", #'Macse Mik Similarity','Macse Mik overlap Identity',
  #"Kud AA% Identitiy over Smorf Frame",
  "Kud overlap AA% Identitiy",#'Macse Kud Similarity','Macse Kud ovRplot02erlap Identity',

  #"Bay AA% Identitiy over Smorf Frame",
  "Bay overlap AA% Identitiy",#'Macse Bay Similarity','Macse Bay overlap Identity',
  "(Para length Ratio)" ,"(Mik length Ratio)" ,"(Kud length Ratio)", "(Bay length Ratio)" ,

  "Para DNA similarity","Mik DNA similarity","Kud DNA similarity","Bay DNA similarity"
)

#subset the dataframe####
#take subset of data excluding (last column, 'is gene?') , (first 2 columns "Orf Name","ORF AA Length")

P <-unique(P,by=P$`Orf Name`)
P_matrix <- (P[,3:(ncol(P)-1)])

row.names(P_matrix)<-P$`Orf Name`

#remove the protogenes that are not in the pccfiltered file from P_matrix
#mustfiltered<-pgs_nodes$orf_name[pgs_nodes$orf_name%in%pccFiltered$orf_name==FALSE]
#P_matrix<-P_matrix[which(rownames(P_matrix)%in%mustfiltered==FALSE),]
P_matrix<-P_matrix[,col_order]
#colnames(P_matrix)[37:40]<-c('Para length Ratio','Mik length Ratio', 'Kud length Ratio','Bay length Ratio') #to remove'Log' from Log(* length ratio) column names if no log operation applied

#P_matrix<-P_matrix[,-grep('Macse start codon',colnames(P_matrix))]
#P_matrix<-P_matrix[,-grep('Macse stop codon',colnames(P_matrix))]
#P_matrix<-P_matrix[,-grep('overlap Identity',colnames(P_matrix))]


#create color palette####
library(RColorBrewer)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(12)

#clustering####
#calculate distance matrix using gower.dist
gower_dist<-gower.dist(P_matrix)


#cluster the distance matrix with agnes. I use 'complete' or 'ward' linkage methods depending on the clusters and goals
agn1<-agnes(gower_dist, method = 'average')

#create dendrogram
dend1_<-as.dendrogram(agn1)
cols_branches <- c("red", "gold", "violet", 'darkseagreen')
dend1_ <- color_branches(dend1_,k = length(cols_branches), col = cols_branches)


#heatmap####
#change the letter containing columns to numbers to color them on the heatmap. Otherwise all of them becomes NAs
P_matrix1<-P_matrix
{P_matrix1$`start codon in Para`[P_matrix1$`start codon in Para`!='M' & P_matrix1$`start codon in Para`!='I']<-0
P_matrix1$`start codon in Mik`[P_matrix1$`start codon in Mik`!='M' & P_matrix1$`start codon in Mik`!='I']<-0

P_matrix1$`start codon in Para`[P_matrix1$`start codon in Para`=='M']<-1
P_matrix1$`start codon in Mik`[P_matrix1$`start codon in Mik`=='M']<-1

P_matrix1$`start codon in Para`[P_matrix1$`start codon in Para`=='I']<-0.5
P_matrix1$`start codon in Mik`[P_matrix1$`start codon in Mik`=='I']<-0.5

P_matrix1$`stop codon in Para`[P_matrix1$`stop codon in Para`!='X']<-0
P_matrix1$`stop codon in Mik`[P_matrix1$`stop codon in Mik`!='X']<-0

P_matrix1$`stop codon in Para`[P_matrix1$`stop codon in Para`=='X']<-1
P_matrix1$`stop codon in Mik`[P_matrix1$`stop codon in Mik`=='X']<-1

P_matrix1$`start codon in Bay`[P_matrix1$`start codon in Bay`!='M' & P_matrix1$`start codon in Bay`!='I']<-0
P_matrix1$`start codon in Kud`[P_matrix1$`start codon in Kud`!='M' & P_matrix1$`start codon in Kud`!='I']<-0

P_matrix1$`start codon in Bay`[P_matrix1$`start codon in Bay`=='M']<-1
P_matrix1$`start codon in Kud`[P_matrix1$`start codon in Kud`=='M']<-1

P_matrix1$`start codon in Bay`[P_matrix1$`start codon in Bay`=='I']<-0.5
P_matrix1$`start codon in Kud`[P_matrix1$`start codon in Kud`=='I']<-0.5

P_matrix1$`stop codon in Bay`[P_matrix1$`stop codon in Bay`!='X']<-0
P_matrix1$`stop codon in Kud`[P_matrix1$`stop codon in Kud`!='X']<-0

P_matrix1$`stop codon in Bay`[P_matrix1$`stop codon in Bay`=='X']<-1
P_matrix1$`stop codon in Kud`[P_matrix1$`stop codon in Kud`=='X']<-1
}

#to rotate the dendrogram. see ?click_rotate
dend1_<-click_rotate(dend1_)
gene_col<-ifelse(rownames(P_matrix)%in%ppp,'lightblue',ifelse(P$`Is gene`==3,'blue','white'))

pdf('../data/plots/heatmap_withnongene.pdf',width=16,height=12)
out<-heatmap.2(data.matrix(P_matrix1),
               #cellnote = mat_data,  # same data set for cell labels
               #main = "agnes, manhattan", # heat map title
               notecol="black",      # change font color of cell labels to black
               density.info="none",  # turns off density plot inside color legend
               trace="none",         # turns off trace lines inside the heat map
               #breaks = palette.breaks,
               #scale = "row",
               margins =c(12,6),     # widens margins around plot
               col=hmcol,       # use on color palette defined earlier
               #breaks=col_breaks,    # enable color transition at specified limits
               srtCol = 45,
               dendrogram="row",     # only draw a row dendrogram
               Rowv=dend1_,#sort(rotate.dendrogram1(dend1_,groups)),#,type = "average"),
              # RowSideColors=matrix(gene_col),
               Colv=FALSE,
               cexRow = 0.09,
               #rowsep = c(531),
               #colsep = c(3,6,9,12,15,19,23,27,30,33),
               sepcolor = 'red',
               sepwidth = c(0.05,1))
dev.off()
#ggsave('data/plots/heatmap1.pdf',width=16,height=12)

#find conserved####
#here my conserved gene columns were contained by dend1_[[2]][[2]] and dend1_[[2]][[1]]
#these will change according to usage of click_rotate so you either should find another way to get labels of clusters or try and figure out
conserved1<-rownames(P_matrix1[as.numeric(labels(dend1_[[2]][[2]])),])[rownames(P_matrix1[as.numeric(labels(dend1_[[2]][[2]])),])%in%pgs_nodes$orf_name]
conserved2<-rownames(P_matrix1[as.numeric(labels(dend1_[[2]][[1]])),])[rownames(P_matrix1[as.numeric(labels(dend1_[[2]][[1]])),])%in%pgs_nodes$orf_name]

t<-rownames(P_matrix1[as.numeric(labels(dend1_[[2]][[1]])),])[rownames(P_matrix1[as.numeric(labels(dend1_[[2]][[1]])),])%in%geneDataTable$ORF.Name]
#[1] "YLR390W" "YAL059W" "YGR084C" "YBR069C" "YER019W" "YBR194W"
t <- c("YLR390W", "YAL059W" ,"YGR084C" ,"YBR069C" ,"YER019W" ,"YBR194W",'YBR243C')

#only protogenes clustering####
P_matrix2<-P_matrix[which(rownames(P_matrix)%in%pgs_nodes$orf_name),]


gower_dist2<-gower.dist(P_matrix2)
agn2<-agnes(gower_dist2, method = 'average')

dend1_2<-as.dendrogram(agn2)
cols_branches <- c("red", "gold", "violet", 'darkseagreen')
dend1_2 <- color_branches(dend1_2,k = length(cols_branches), col = cols_branches)

#heatmap only protogenes####

{P_matrix2$`start codon in Para`[P_matrix2$`start codon in Para`!='M' & P_matrix2$`start codon in Para`!='I']<-0
P_matrix2$`start codon in Mik`[P_matrix2$`start codon in Mik`!='M' & P_matrix2$`start codon in Mik`!='I']<-0

P_matrix2$`start codon in Para`[P_matrix2$`start codon in Para`=='M']<-1
P_matrix2$`start codon in Mik`[P_matrix2$`start codon in Mik`=='M']<-1

P_matrix2$`start codon in Para`[P_matrix2$`start codon in Para`=='I']<-0.5
P_matrix2$`start codon in Mik`[P_matrix2$`start codon in Mik`=='I']<-0.5

P_matrix2$`stop codon in Para`[P_matrix2$`stop codon in Para`!='X']<-0
P_matrix2$`stop codon in Mik`[P_matrix2$`stop codon in Mik`!='X']<-0

P_matrix2$`stop codon in Para`[P_matrix2$`stop codon in Para`=='X']<-1
P_matrix2$`stop codon in Mik`[P_matrix2$`stop codon in Mik`=='X']<-1

P_matrix2$`start codon in Bay`[P_matrix2$`start codon in Bay`!='M' & P_matrix2$`start codon in Bay`!='I']<-0
P_matrix2$`start codon in Kud`[P_matrix2$`start codon in Kud`!='M' & P_matrix2$`start codon in Kud`!='I']<-0

P_matrix2$`start codon in Bay`[P_matrix2$`start codon in Bay`=='M']<-1
P_matrix2$`start codon in Kud`[P_matrix2$`start codon in Kud`=='M']<-1

P_matrix2$`start codon in Bay`[P_matrix2$`start codon in Bay`=='I']<-0.5
P_matrix2$`start codon in Kud`[P_matrix2$`start codon in Kud`=='I']<-0.5

P_matrix2$`stop codon in Bay`[P_matrix2$`stop codon in Bay`!='X']<-0
P_matrix2$`stop codon in Kud`[P_matrix2$`stop codon in Kud`!='X']<-0

P_matrix2$`stop codon in Bay`[P_matrix2$`stop codon in Bay`=='X']<-1
P_matrix2$`stop codon in Kud`[P_matrix2$`stop codon in Kud`=='X']<-1
}

dend1_2<-click_rotate(dend1_2)


cons_col<-ifelse(rownames(P_matrix2)%in%conserved1,'violet',ifelse(rownames(P_matrix2)%in%conserved2,'darkseagreen','white'))
out<-heatmap.2(data.matrix(P_matrix2[,colnames(P_matrix1)]),
               #cellnote = mat_data,  # same data set for cell labels
               #main = "agnes, manhattan", # heat map title
               notecol="black",      # change font color of cell labels to black
               density.info="none",  # turns off density plot inside color legend
               trace="none",         # turns off trace lines inside the heat map
               #breaks = palette.breaks,
               #scale = "row",
               margins =c(12,6),     # widens margins around plot
               col=hmcol,       # use on color palette defined earlier
               #breaks=col_breaks,    # enable color transition at specified limits
               srtCol = 45,
               dendrogram="row",     # only draw a row dendrogram
               Rowv=dend1_2,#sort(rotate.dendrogram1(dend1_,groups)),#,type = "average"),
               RowSideColors=matrix(cons_col),
               Colv=FALSE,
               cexRow = 0.09,
               #rowsep = c(531),
               #colsep = c(3,6,9,12,15,19,23,27,30,33),
               sepcolor = 'red',
               sepwidth = c(0.05,1))

#find conserved.2####
#again the dend1_2 parts are subject to change with click_rotate().
conserved_combined<-rownames(P_matrix2[as.numeric(labels(dend1_2[[1]][[2]])),])
afterfilters<-setdiff(rownames(P_matrix2),conserved_combined)
#write.csv(afterfilters,file = 'afterfilterProtogenes.csv',quote = FALSE,row.names = FALSE)
