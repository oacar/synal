library(muscle)
library(msa)
#source('functions.R')
library(synal)
options(warn = 1)
#to remove all alternative start codons. Bc I wanted to only ATG to be translated as M, the defaults were different, you can read ?GENETIC_CODE
attr(GENETIC_CODE, "alt_init_codons")=character(0)
#IMPORTANT
#The syntenic blocks should have the proto-gene orf as starting with ATG. This is because all ORFs are reported this way in the ORF file I have downloaded from SGD.
#and proto-gene orf starting with ATG should be added to DNAStringSet object as last sequence to get subalignment and aa translations
# msg <- file("nongenes_alignment_message_1.Rout", open="wt")
# out <- file("nongenes_alignment_output_1.Rout", open="wt")
# sink(msg, type="message")
# sink(out, type="output")
path<-''
newPath<-''

mainFile<-'data/alignmentInput/'
p_<-paste(mainFile,'sequenceFilescopy/',sep = '')
##these are to copy problematic files to different folders
#npp and noorfwithstartandstop are explained in 'Find orf for aligned regions' section
npp<-paste(mainFile,'/','npp',sep = '')
noorfwithstartandstop<-paste(mainFile,'/','noorfwithmandx/',sep = '')

#this script and whole functions are based on having (scer-spar-smik-skud-sbay) sequence files. so everything is assumed in this order. You can use less number of sequences but for example whatever is on the 3rd sequence will be treated as Mik sequence
#In order not to make wrong analysis, when one of the sequences is not available, the folder is moved to 'noseq' folder and the alignment is not performed.
noseq<-paste(mainFile,'/','noseq',sep = '')

#muscle cannot align sequences if one of them has no nucleotides. So in that case an appropriate warning message is shown and the folder is moved to noregion folder
noregion<-paste(mainFile,'/','noregion',sep = '')
#if proto-gene ORF cannot be found on scer sequence, the folder is moved to noseqrange
noseqrange<-paste(mainFile,'/','noseq/noSeqRange',sep = '')
#if all above not the case, then the folder is moved to noproblem folder
noproblem<-paste(mainFile,'/','noProblem/',sep = '')

#all the above is to help tracking the possible problems related with the alignments. so the script is very much open to improvements
iter<-0
dirName<-'YAL016C-B'
dir.create(noorfwithstartandstop)
dir.create(noseq)
dir.create(noregion)
dir.create(npp)
dir.create(noseqrange)
dir.create(noproblem)

for ( dirName in list.files(p_)){
  file.rename(path, newPath)
  path<-paste(p_,dirName, sep="/")
  newPath<-paste(noproblem,dirName, sep="/")
  noseqPath<-paste(noseq,dirName,sep = '/')
  noregionPath<-paste(noregion,dirName,sep = '/')
  nppPath<-paste(npp,dirName,sep = '/')

  #read files####
  mySequences<-readFiles(path)

  print(dirName)
  #current script is to align 5 sequences and proto-gene sequence being the 6th. to get the subalignment you need the proto-gene sequence
  t<-4
  if(length(mySequences)!=t){
    warning(paste(dirName,' ','File is missing one of the ',t-1,' sequences',sep = ' '))
    file.rename(path, noseqPath)
    next
  }
  #alignment####
  dnaAlignmentList<-synal::alignWoSmorf(mySequences)

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
  if(start<=20){
    sStart<-1
  }else{
    sStart<-start-20
  }
  if(abs(stop-length(DNAStr[[1]]))<=20){
    sStop<-length(DNAStr[[1]])
  }else{
    sStop<-stop+20
  }
  smallerDNAStr<-subalignment(sStart,sStop, DNAStr)[[1]]
  smallerDNAStr<-append(smallerDNAStr,mySequences[length(mySequences)])

  sDnaAlignmentList<-synal::alignWoSmorf(smallerDNAStr)

  sDNAStr<-sDnaAlignmentList[[1]]
  sStart<-sDnaAlignmentList[[2]]
  sStop<-sDnaAlignmentList[[3]]
  writeXStringSet(DNAStr, file=paste(paste(path,dirName, sep="/"),"alignment.fa",sep = "_"))
  i_<-length(DNAStr)


  #start <-findSmorfFrame(DNAStr[[i_]],mySequences[i_])[1]
  #stop <-findSmorfFrame(DNAStr[[i_]],mySequences[i_])[2]

  # subalignment####
  subalign<-subalignment(sStart,sStop, sDNAStr)[[1]]
  subalign<-append(subalign, mySequences[length(mySequences)])
  #subalign<-DNAStringSet(msa(subalign, 'Muscle', order = 'input'))
  #subalign<-chartr('N','-',subalign)warni
  subalign<-tryCatch(DNAStringSet(muscle(subalign)), error=function(e) FALSE)
  if(is.logical(subalign)){
    warning(paste("One of the sequences does not have any nucleotide in the aligned region for ",dirName,sep = ''))
    newPath<-paste(noregion,dirName, sep="/")
    if(TRUE){next}
  }
  identity_percentage<-subalignment(1,length(subalign[[1]]), subalign)[[2]][length(mySequences)]

  if(identity_percentage!=1){
    stop(paste('Identity percentage is not 100% for ', dirName,sep = ''))
  }
  writeXStringSet(subalign, file=paste(paste(path,dirName, sep="/"),"subalignment.fa",sep = "_"))

  #Aminoacid translation####
  aa_alignment<-aaTranslation(subalign,DNAStr)
  if(is.logical(aa_alignment)){
    warning(paste("One of the sequences does not have any nucleotide in the aligned region for ",dirName,sep = ''))
    newPath<-paste(noregion,dirName, sep="/")
    if(TRUE){next}
  }
  writeXStringSet(aa_alignment,file=paste(paste(path,dirName, sep="/"),"AATranslation.fa",sep = "_"))

  #2nd and 3rd frame translation####
  #take scer orf and change the translation frame for other species
  #KNOWN ISSUE: if the first 2 or more nucleotides are gaps, then this is not working but since I don't use 2nd and 3rd frame translations at the moment, I didnot worked on it
  subalign2<-subalignment(start,stop+1, DNAStr)[[1]]
  for(itr in 2:(length(subalign)-1)){
    subalign2[[itr]]<-subalign2[[itr]][2:length(subalign2[[1]])]
  }
  subalign2<-subalignment(1,length(subalign2[[1]])-1, subalign2)[[1]]

  subalign3<-subalignment(start,stop+2, DNAStr)[[1]]
  for(itr in 2:(length(subalign)-1)){
    subalign3[[itr]]<-subalign3[[itr]][3:length(subalign3[[1]])]
  }
  subalign3<-subalignment(1,length(subalign3[[1]])-2, subalign3)[[1]]

  aa_alignment2<-aaTranslation(subalign2,DNAStr)
  aa_alignment3<-aaTranslation(subalign3,DNAStr)
  writeXStringSet(aa_alignment2,file=paste(paste(path,dirName, sep="/"),"AATranslation_2ndFrame.fa",sep = "_"))
  writeXStringSet(aa_alignment3,file=paste(paste(path,dirName, sep="/"),"AATranslation_3rdFrame.fa",sep = "_"))



  #Find orf for aligned regions####
  #This for loop checks if ORF exists in one of the following species, if not searchs for start or stop codons, whichever is missing
  #if the search causes an error the (find*) functions returns NULL and the problematic alingment folder is moved to 'npp' folder
  #KNOWN ISSUE: If there is no ORF but there is an X before M there should be 2 way search:
  #1-Search for an ORF that ends with that X
  #2-Search for an ORF that starts with the M and ends after the scer orf ends
  #take whichever is longer.    I didn't have time to do this. So " CheckOrf=='Neither' " shows this case and they are moved to 'noorfwithstartandstop' folder

  types<-c('scer','para','mik','kud','bay')
  for(j in 2:length(DNAStr)){
  PrevStart<-start
  NextStop<-stop
  CheckOrf<-checkORF(aa_alignment[[j]])

  #####If no start is found prevStartPara returns inf, check for that!!!!!!!!!!!!
  if(CheckOrf!=TRUE){
    warning(paste('No ORF found in ',types[j],'! for ', dirName, sep = ''))
    if(CheckOrf=='BOTH'){
      find.orf<-findBoth(DNAStr, j, start, stop, dirName)
      if(is.null(find.orf$writeFile)){
        newPath<-paste(npp,dirName, sep="/")
        next
      }
    }else if(CheckOrf=='M'){

      find.orf<-findPrevStart(DNAStr,j,start,stop,dirName)

      if(is.null(find.orf$writeFile)){
        newPath<-paste(npp,dirName, sep="/")
        next
      }
    } else if (CheckOrf=='X'){
      find.orf<-findNextStop(DNAStr,j,start,stop,dirName)
      if(is.null(find.orf$writeFile)){
        newPath<-paste(npp,dirName, sep="/")
        if(TRUE){next}
      }
    } else if(CheckOrf=='Neither'){
      find.orf.start<-findPrevStart(DNAStr,j,start,stop,dirName)
      find.orf.stop<-findNextStop(DNAStr,j,start,stop,dirName)
      if(find.orf.start$writeFile==T && find.orf.stop$writeFile==T){
        find.orf <- ifelse((find.orf.start$NextStop-find.orf.start$PrevStart)>(find.orf.stop$NextStop-find.orf.stop$PrevStart),find.orf.start,find.orf.stop)
      }else if(find.orf.start$writeFile==T && is.null(find.orf.stop$writeFile)==T){
        find.orf<-find.orf.start
      }else if(is.null(find.orf.start$writeFile)==T && find.orf.stop$writeFile==T){
        find.orf <- find.orf.stop
      }else{
        newPath<-paste(noorfwithstartandstop,dirName,sep = '/')
        if(TRUE){next}
      }
    }


  if (find.orf$writeFile){
    PrevStart<-find.orf$PrevStart
    NextStop<-find.orf$NextStop
    subalign.p<-subalignment(PrevStart, NextStop,DNAStr)
    aatranslation.p<-aaTranslation(subalign.p[[1]], DNAStr)
    writeXStringSet(subalign.p[[1]], file=paste(paste(path,dirName, sep="/"),"_subalignment_",types[j],".fa",sep = ""))
    writeXStringSet(aatranslation.p,file=paste(paste(path,dirName, sep="/"),"_AATranslation_",types[j],".fa",sep = ""))

      }
    }
  }
}
