

findSeq<- function(path, geneName, outputPath, negativeStrand, smorfName){
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
  if(grepl('cer',path)){
    type_<-'Scer'
  } else if(grepl('par',path)){
    type_<-'Spar'
  }else if(grepl('mik',path)){
    type_<-'Smik'
  }else if(grepl('bay',path)){
    type_<-'Sbay'
  }else if(grepl('kud',path)){
    type_<-'Skud'
  }
  seqfile<-(readDNAStringSet(path))

  #if(type_ <- )
  for (i in 1:length(seqfile)){
    if(geneName==strsplit(strsplit(names(seqfile)[i],',')[[1]][1],' ')[[1]][1] || (is.na(strsplit(strsplit(names(seqfile)[i],',')[[1]][1],' ')[[1]][2])==F &geneName==strsplit(strsplit(names(seqfile)[i],',')[[1]][1],' ')[[1]][2])){
      index<-i
    }
  }
  if(index==0){
    if(type_=='Scer'){
      warning(sprintf("%s", paste(geneName, 'cannot be found for',smorfName, "for", type_, sep=" ")))}
    else{
      if(type_=='Spar'){
        #NOTFOUNDSPAR<<-NOTFOUNDSPAR+1
       # return(FALSE)
        #  print(NOTFOUNDSPAR)
      }else {
        #NOTFOUNDSMIK<<-NOTFOUNDSMIK+1
        #return(FALSE)
        # print(NOTFOUNDSPAR)

      }
      warning(sprintf("%s", paste(geneName, 'cannot be found for', smorfName, "for", type_, sep=" ")))
      return(FALSE)
    }
  }else{  Seq<-seqfile[index]

  fileName<-sprintf("%s", paste(outputPath,geneName, sep="/"))
  fileName<-sprintf("%s", paste(fileName,type_, sep="_"))
  fileName<-sprintf("%s", paste(fileName,'.fa', sep=""))
  if(negativeStrand){
    writeXStringSet(reverseComplement(Seq), file=fileName)
  } else{
    writeXStringSet(Seq, file=fileName)

  }
  }
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


findSmorfCoor<-function(smorfSeq, ChromosomeSeq, Name){
  match<-matchPattern(smorfSeq,ChromosomeSeq[[1]])
  if(length(start(match))==0){
    start<-NULL
    end<-NULL
    NUMBEROFNOTFOUND<<-NUMBEROFNOTFOUND+1
    NOTFOUNDSMORF<<-append(NOTFOUNDSMORF, Name)
    warning(sprintf("%s", paste("Smorf Sequence for",Name,"was not found on the Chromosome it might be changed", sep=" ")))
  }else{
    start<-start(match)[1]
    end<-end(match)[1]}
  return(list(start,end))
}
writeSmorfSeq <-function(smorfSeq, Name, path){
  #take SmorfSeq as string, smorfName as String and output Path as string
  #does not return anything, saves smorfSeq as fasta file to the path
  smorf<-DNAStringSet(smorfSeq)
  names(smorf)<-Name
  #path<-sprintf("%s", paste(path,"smorf_", sep="/"))
  path<-sprintf("%s", paste(path,Name, sep="/"))
  path<-sprintf("%s", paste(path,"fa", sep="."))

  writeXStringSet(smorf,path )
}

findNearestGene <- function(beg, end,chr,path, name) {
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
  dummy<-c()
  coordinates<-c()
  startVec<-c()
  stopVec<-c()
  chromosomeVec<-c()
  namesVec<-c()
  nearestList<-c()
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
  ##startVec has start coordinates of genes
  ##stopVec has end coordinates of genes
  min<-Inf
  minIndex<-0
  for(i in 1:length(startVec)){
    if(chromosomeVec[i]==chr && namesVec[i]!=name){
      begToStart<-abs(beg-startVec[i])
      endToStart<-abs(end-startVec[i])
      begToStop<-abs(beg-stopVec[i])
      endToStop<-abs(end-stopVec[i])
      if(begToStart<min){
        min<-begToStart
        minIndex<-i
      }
      if(endToStart<min){
        min<-endToStart
        minIndex<-i
      }
      if(begToStop<min){
        min<-begToStop
        minIndex<-i
      }
      if(endToStop<min){
        min<-endToStop
        minIndex<-i
      }
    }}
  if(min<1000){
    ng<-scer[minIndex]
    ng_name<-strsplit(strsplit(names(ng), ',')[[1]][1],' ')[[1]][1]
    return(list(scer[minIndex], ng_name))
  }
  else{
    return(FALSE)
  }
}
findNearestGenes <- function(beg, end,chr,path, prevNGene=NULL) {
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
  dummy<-c()
  coordinates<-c()
  startVec<-c()
  stopVec<-c()
  chromosomeVec<-c()
  #namesVec<-c()
  nearestList<-c()
  minList<-c()
  for (i in 1:length(scer)){
    if(i==410 || i==1574 || i == 2042){
      # namesVec<-c(namesVec, strsplit(strsplit(names(scer),',')[[1]][1],' ')[[1]][1])
      dummy<-c(dummy,strsplit(names(scer)[i],',')[[1]][3])
      coordinates<-c(coordinates,strsplit(dummy[i], ' ')[[1]][5])
      chromosomeVec<-c(chromosomeVec,strsplit(dummy[i], ' ')[[1]][3])
      startVec<-c(startVec, strsplit(coordinates[i], '-')[[1]][1])
      stopVec<-c(stopVec, strsplit(coordinates[i], '-')[[1]][2])
    }
    else{

      dummy<-c(dummy,strsplit(names(scer)[i],',')[[1]][2])
      coordinates<-c(coordinates,strsplit(dummy[i], ' ')[[1]][5])
      chromosomeVec<-c(chromosomeVec,strsplit(dummy[i], ' ')[[1]][3])
      startVec<-c(startVec, strsplit(coordinates[i], '-')[[1]][1])
      stopVec<-c(stopVec, strsplit(coordinates[i], '-')[[1]][2])

    }}
  startVec<-as.integer(startVec)
  stopVec<-as.integer(stopVec)
  ##startVec has start coordinates of genes
  ##stopVec has end coordinates of genes
  min<-Inf
  minIndex<-0
  for(i in 1:length(startVec)){
    if(chromosomeVec[i]==chr){
      begToStart<-abs(beg-startVec[i])
      endToStart<-abs(end-startVec[i])
      begToStop<-abs(beg-stopVec[i])
      endToStop<-abs(end-stopVec[i])

      if(min(begToStart,endToStart,begToStop,endToStop)<min){
        min<-min(begToStart,endToStart,begToStop,endToStop)
        minIndex<-i
        if(is.null(nearestList) || all(minIndex!=nearestList) ){
          nearestList<-append(nearestList,minIndex)
          minList<-append(minList,min)
          if(length(nearestList)>2){
            nearestList<-nearestList[2:3]
            minList<-minList[2:3]
          }
        }
      }
    #   if(begToStart<min){
    #     min<-begToStart
    #     minIndex<-i
    #     if(is.null(nearestList) || all(minIndex!=nearestList) ){
    #       nearestList<-append(nearestList,minIndex)
    #       minList<-append(minList,min)
    #       if(length(nearestList)>2){
    #         nearestList<-nearestList[2:3]
    #         minList<-minList[2:3]
    #       }
    #     }
    #   }
    #   if(endToStart<min){
    #     min<-endToStart
    #     minIndex<-i
    #     if(is.null(nearestList) || all(minIndex!=nearestList) ){
    #       nearestList<-append(nearestList,minIndex)
    #       minList<-append(minList,min)
    #       if(length(nearestList)>2){
    #         nearestList<-nearestList[2:3]
    #         minList<-minList[2:3]
    #       }
    #     }
    #
    #   }
    #   if(begToStop<min){
    #     min<-begToStop
    #     minIndex<-i
    #     if(is.null(nearestList) || all(minIndex!=nearestList) ){
    #       nearestList<-append(nearestList,minIndex)
    #       minList<-append(minList,min)
    #       if(length(nearestList)>2){
    #         nearestList<-nearestList[2:3]
    #         minList<-minList[2:3]
    #       }
    #     }
    #
    #   }
    #   if(endToStop<min){
    #     min<-endToStop
    #     minIndex<-i
    #     if(is.null(nearestList) || all(minIndex!=nearestList) ){
    #       nearestList<-append(nearestList,minIndex)
    #       minList<-append(minList,min)
    #       if(length(nearestList)>2){
    #         nearestList<-nearestList[2:3]
    #         minList<-minList[2:3]
    #       }
    #     }
    #   }
     }}
    #
  if(all(minList<1000) && length(minList)==2){
    ng<-c(scer[nearestList[1]], scer[nearestList[2]])
    ng_name<-c(strsplit(strsplit(names(ng[1]), ',')[[1]][1],' ')[[1]][1],strsplit(strsplit(names(ng[2]), ',')[[1]][1],' ')[[1]][1])
    return(list(ng, ng_name))
  }else if(length(minList)==2 && minList[2]<1000){
    ng<-scer[nearestList[2]]
    ng_name<-strsplit(strsplit(names(ng), ',')[[1]][1],' ')[[1]][1]
    return(list(ng, ng_name))
  } else if(length(minList)==1 && minList[1]<1000){
    ng<-scer[nearestList[1]]
    ng_name<-strsplit(strsplit(names(ng), ',')[[1]][1],' ')[[1]][1]
    return(list(ng, ng_name))
  }else{
    return(FALSE)
  }
}
findChr<-function(chr, path){
  ##Takes chromosome number and Path of chromosome file
  ## returns DNAStringSet object of the chromosome
  if(dir.exists(path)!=TRUE){
    stop('Path does not exist. Give a correct path')
  }
  chrFName<-''
  for (fname in list.files(path)){
    if (grepl(chr, fname, ignore.case = TRUE)){
      chrFName <-fname
    }
    if(chrFName!=''){
      break
    }
  }
  genomeSeq<-readDNAStringSet(sprintf("%s", paste(path,chrFName, sep="/")))
  return(genomeSeq)
}
