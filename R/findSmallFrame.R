#'Find small ORF frame
#'Finds the ORF region over the other species if there is a smaller ORF, or just returns the overlapping parts of the Scer and the Other species' ORF
#'@param aaAlignment aligned AAStringSet object
#'@param whichSeq string spesifying the other species name. c('Para','Mik', 'Bay', 'Kud') are acceptable
#'@param seqAA.alignment if trying to find overlapping region when there is ....X......M..... situation where X comes before M, this will be used and is the AAString instance for that frame AA
#'@param overlap is false by default. Used only in \code{findBestOverlap} to distinguish with normal usage of function. When this is true, and seqAA.alignment is needed, it will only return the overlapping region from start of sequence to X codon
<<<<<<< HEAD
#'@param longer is false by default. Used only in \code{findBestOverlap} for now. It is to distinguish actual start codon with other start codons in the overlapping region to find correct overlap. It is TRUE if the sequence contains real start codon of proto-gene. False if not
=======
#'@param longer is false by default. Used only in \code{findBestOverlap} for now. It is to distinguish actual start codon with other start codons in the overlapping region to find correct overlap
>>>>>>> c1655aa13be547d7cec4b70c6fa5291cc3ca1da2
#'@return small orf alignment or overlapping region of alignment
#' returns FALSE if alignment causes an error. Should be checked on main script properly
#' @import Biostrings
#' @import muscle
#'@export

<<<<<<< HEAD
findSmallFrame <- function(aaAlignment, fourCodons, seqAA.alignment=NULL,overlap=F,longer=F) {
=======
findSmallFrame <- function(aaAlignment, whichSeq, seqAA.alignment=NULL,overlap=F,longer=F) {
>>>>>>> c1655aa13be547d7cec4b70c6fa5291cc3ca1da2
  out<-tryCatch(
    {
      # id<--1
      # if(is.character(whichSeq)){
      #   if(whichSeq=='Para'){
      #     id<-2
      #   }else if(whichSeq=='Mik'){
      #     id<-3
      #   }else if(whichSeq=='Bay'){
      #     id<-5
      #   }else if(whichSeq=='Kud'){
      #     id<-4
      #   }else{
      #     stop('Specify which sequence to be found as \'Para\' or \'Mik\' or \'Bay\' or \'Kud\' ')
      #   }
      # }else{
      #   id <- whichSeq
      # }
      id<-2

      aaAlignment<-aaAlignment[c(1,id)]
      ##to calculate identity over overlapping region of scer and para\mik orfs,
      ##this function now returns the only the part overlaps with scer orf
      orf<-checkORF(aaAlignment[[id]])
<<<<<<< HEAD
      if(longer){
        #coors<-findORF(aaAlignment[[id]])
        startPos<-str_locate(aaAlignment[[2]],as.character(fourCodons))[1,1]#ifelse(longer,coors[[1]],1)
        xpos <- str_locate(aaAlignment[[2]],'X')[1,1]
        xpos2 <- str_locate(aaAlignment[[1]],'X')[1,1]
        stopPos<-ifelse(is.na(xpos),length(aaAlignment[[1]]),
                        ifelse(xpos<startPos,
                               ifelse(nrow(str_locate_all(aaAlignment[[2]],'X')[[1]])==2,
                                      ifelse(str_locate_all(aaAlignment[[2]],'X')[[1]][2,1]<xpos2,str_locate_all(aaAlignment[[2]],'X')[[1]][2,1],xpos2),length(aaAlignment[[1]])),xpos))
=======
      if(orf==TRUE){
        coors<-findORF(aaAlignment[[id]])
        startPos<-ifelse(longer,coors[[1]],1)
        stopPos<-coors[[2]]
>>>>>>> c1655aa13be547d7cec4b70c6fa5291cc3ca1da2

        if(startPos==-1 || stopPos==-1){
          return(FALSE)
        }
        aaSmall<-subseq(aaAlignment,startPos, stopPos)
        aaSmall<-AAStringSet(muscle(aaSmall,quiet = T))
<<<<<<< HEAD
        # startPos<-ifelse(longer,checkCodon(aaSmall[[id]],'M'),1)
        # # if(startPos==-1){
        # #   startPos<-checkCodon(aaSmall[[id]],'I')
        # # }
        # stopPos<-checkCodon(aaSmall[[id]],'X')
        # aaSmall<-subseq(aaSmall,startPos, stopPos)
        # aaSmall<-AAStringSet(muscle(aaSmall,quiet = T))
        return(aaSmall)
      }else{
        mpos2 <- str_locate(aaAlignment[[1]],'M')[1,1]
        startPos<-ifelse(mpos2>1,mpos2,1)
        xpos <- str_locate(aaAlignment[[2]],'X')[1,1]
        xpos2 <- str_locate(aaAlignment[[1]],'X')[1,1]

        stopPos<-ifelse(is.na(xpos),length(aaAlignment[[1]]),ifelse(xpos>xpos2,xpos2,xpos))
        if(startPos==-1 || stopPos==-1){
          return(FALSE)
        }
        aaSmall<-subseq(aaAlignment,startPos, stopPos)
        aaSmall<-AAStringSet(muscle(aaSmall,quiet = T))
        return(aaSmall)

      }
    },
=======
        startPos<-ifelse(longer,checkCodon(aaSmall[[id]],'M'),1)
        # if(startPos==-1){
        #   startPos<-checkCodon(aaSmall[[id]],'I')
        # }
        stopPos<-checkCodon(aaSmall[[id]],'X')
        aaSmall<-subseq(aaSmall,startPos, stopPos)
        aaSmall<-AAStringSet(muscle(aaSmall,quiet = T))
        return(aaSmall)
      }else{
        if(orf=="BOTH"){
          return(aaAlignment)
        }else if(orf=="M"){
          startPos<-1
          stopPos<-checkCodon(aaAlignment[[id]],'X')

          aaSmall<-subseq(aaAlignment,startPos, stopPos)
          aaSmall<-AAStringSet(muscle(aaSmall,quiet = T))

          stopPos<-checkCodon(aaSmall[[id]],'X')
          aaSmall<-subseq(aaSmall,startPos, stopPos)
          aaSmall<-AAStringSet(muscle(aaSmall,quiet = T))
          return(aaSmall)
        }else if(orf=="X"){
          startPos<-checkCodon(aaAlignment[[id]],'M')
          # if(startPos==-1){
          #   startPos<-checkCodon(aaAlignment[[id]],'I')
          # }
          stopPos<-length(aaAlignment[[1]])

          aaSmall<-subseq(aaAlignment,startPos, stopPos)
          aaSmall<-AAStringSet(muscle(aaSmall,quiet = T))
          startPos<-checkCodon(aaSmall[[id]],'M')
          # if(startPos==-1){
          #   startPos<-checkCodon(aaSmall[[id]],'I')
          # }
          stopPos<-length(aaSmall[[1]])
          aaSmall<-subseq(aaSmall,startPos, stopPos)
          aaSmall<-AAStringSet(muscle(aaSmall,quiet = T))
          return(aaSmall)
        }else if(orf=="Neither"){
          if(is.null(seqAA.alignment)==F){
            if(as.character(seqAA.alignment[1])=='M' | overlap){
              startPos<-1
              stopPos<-checkCodon(aaAlignment[[id]],'X')

              aaSmall<-subseq(aaAlignment,startPos, stopPos)
              aaSmall<-AAStringSet(muscle(aaSmall,quiet = T))

              stopPos<-checkCodon(aaSmall[[id]],'X')
              aaSmall<-subseq(aaSmall,startPos, stopPos)
              aaSmall<-AAStringSet(muscle(aaSmall,quiet = T))
              return(aaSmall)
            }else{
              startPos<-checkCodon(aaAlignment[[id]],'M')
              # if(startPos==-1){
              #   startPos<-checkCodon(aaAlignment[[id]],'I')
              # }
              stopPos<-length(aaAlignment[[1]])

              aaSmall<-subseq(aaAlignment,startPos, stopPos)
              aaSmall<-AAStringSet(muscle(aaSmall,quiet = T))
              startPos<-checkCodon(aaSmall[[id]],'M')
              # if(startPos==-1){
              #   startPos<-checkCodon(aaSmall[[id]],'I')
              # }
              stopPos<-length(aaSmall[[1]])
              aaSmall<-subseq(aaSmall,startPos, stopPos)
              aaSmall<-AAStringSet(muscle(aaSmall,quiet = T))
              return(aaSmall)
            }
          }else{return(FALSE)}
        }else{return(FALSE)}
      }},
>>>>>>> c1655aa13be547d7cec4b70c6fa5291cc3ca1da2
    error=function(cond){
      message('No region on one of sequence found!')
      return(FALSE)
    })
  return(out)
}
#
# findSmallFrame <- function(aaAlignment, whichSeq, seqAA.alignment=NULL) {
#   out<-tryCatch(
#     {id<--1
#     if(whichSeq=='Para'){
#       id<-2
#     }else if(whichSeq=='Mik'){
#       id<-3
#     }else if(whichSeq=='Bay'){
#       id<-5
#     }else if(whichSeq=='Kud'){
#       id<-4
#     }else{
#       stop('Specify which sequence to be found as \'Para\' or \'Mik\' or \'Bay\' or \'Kud\' ')
#     }
#     aaAlignment<-aaAlignment[c(1,id)]
#     id<-2
#     ##to calculate identity over overlapping region of scer and para\mik orfs,
#     ##this function now returns the only the part overlaps with scer orf
#     orf<-checkORF(aaAlignment[[id]])
#     if(orf==TRUE){
#       coors<-findORF(aaAlignment[[id]])
#       startPos<-coors[[1]]
#       stopPos<-coors[[2]]
#
#       if(startPos==-1 || stopPos==-1){
#         return(FALSE)
#       }
#       aaSmall<-subseq(aaAlignment,startPos, stopPos)
#       aaSmall<-AAStringSet(muscle(aaSmall))
#       startPos<-checkCodon(aaSmall[[id]],'M')
#       if(startPos==-1){
#         startPos<-checkCodon(aaSmall[[id]],'I')
#       }
#       stopPos<-checkCodon(aaSmall[[id]],'X')
#       aaSmall<-subseq(aaSmall,startPos, stopPos)
#       aaSmall<-AAStringSet(muscle(aaSmall))
#       return(aaSmall)
#     }else{
#       if(orf=="BOTH"){
#         return(aaAlignment)
#       }else if(orf=="M"){
#         startPos<-1
#         stopPos<-checkCodon(aaAlignment[[id]],'X')
#
#         aaSmall<-subseq(aaAlignment,startPos, stopPos)
#         aaSmall<-AAStringSet(muscle(aaSmall))
#
#         stopPos<-checkCodon(aaSmall[[id]],'X')
#         aaSmall<-subseq(aaSmall,startPos, stopPos)
#         aaSmall<-AAStringSet(muscle(aaSmall))
#         return(aaSmall)
#       }else if(orf=="X"){
#         startPos<-checkCodon(aaAlignment[[id]],'M')
#         if(startPos==-1){
#           startPos<-checkCodon(aaAlignment[[id]],'I')
#         }
#         stopPos<-length(aaAlignment[[1]])
#
#         aaSmall<-subseq(aaAlignment,startPos, stopPos)
#         aaSmall<-AAStringSet(muscle(aaSmall))
#         startPos<-checkCodon(aaSmall[[id]],'M')
#         if(startPos==-1){
#           startPos<-checkCodon(aaSmall[[id]],'I')
#         }
#         stopPos<-length(aaSmall[[1]])
#         aaSmall<-subseq(aaSmall,startPos, stopPos)
#         aaSmall<-AAStringSet(muscle(aaSmall))
#         return(aaSmall)
#       }else if(orf=="Neither"){
#         if(is.null(seqAA.alignment)==F){
#           if(seqAA.alignment[1]=='M' | seqAA.alignment[1]==I){
#             startPos<-1
#             stopPos<-checkCodon(aaAlignment[[id]],'X')
#
#             aaSmall<-subseq(aaAlignment,startPos, stopPos)
#             aaSmall<-AAStringSet(muscle(aaSmall))
#
#             stopPos<-checkCodon(aaSmall[[id]],'X')
#             aaSmall<-subseq(aaSmall,startPos, stopPos)
#             aaSmall<-AAStringSet(muscle(aaSmall))
#             return(aaSmall)
#           }else{
#             startPos<-checkCodon(aaAlignment[[id]],'M')
#             if(startPos==-1){
#               startPos<-checkCodon(aaAlignment[[id]],'I')
#             }
#             stopPos<-length(aaAlignment[[1]])
#
#             aaSmall<-subseq(aaAlignment,startPos, stopPos)
#             aaSmall<-AAStringSet(muscle(aaSmall))
#             startPos<-checkCodon(aaSmall[[id]],'M')
#             if(startPos==-1){
#               startPos<-checkCodon(aaSmall[[id]],'I')
#             }
#             stopPos<-length(aaSmall[[1]])
#             aaSmall<-subseq(aaSmall,startPos, stopPos)
#             aaSmall<-AAStringSet(muscle(aaSmall))
#             return(aaSmall)
#           }
#         }else{return(FALSE)}
#       }else{return(FALSE)}
#     }},
#     error=function(cond){
#       message('No region on one of sequence found!')
#       return(FALSE)
#     })
#   return(out)
# }
