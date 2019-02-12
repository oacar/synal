
#'Find small ORF frame
#'Finds the ORF region over the other species if there is a smaller ORF, or just returns the overlapping parts of the Scer and the Other species' ORF
#'@param aaAlignment aligned AAStringSet object
#'@param whichSeq string spesifying the other species name. c('Para','Mik', 'Bay', 'Kud') are acceptable
#'@return small orf alignment or overlapping region of alignment
#' returns FALSE if alignment causes an error. Should be checked on main script properly
#' @import Biostrings
#' @import muscle
#'@export


#'Find small ORF frame
#'Finds the ORF region over the other species if there is a smaller ORF, or just returns the overlapping parts of the Scer and the Other species' ORF
#'@param aaAlignment aligned AAStringSet object
#'@param whichSeq string spesifying the other species name. c('Para','Mik', 'Bay', 'Kud') are acceptable
#'@return small orf alignment or overlapping region of alignment
#' returns FALSE if alignment causes an error. Should be checked on main script properly
#' @import Biostrings
#' @import muscle
#'@export

findSmallFrame <- function(aaAlignment, whichSeq, seqAA.alignment=NULL) {
  out<-tryCatch(
    {
      id<--1
      if(is.character(whichSeq)){
        if(whichSeq=='Para'){
          id<-2
        }else if(whichSeq=='Mik'){
          id<-3
        }else if(whichSeq=='Bay'){
          id<-5
        }else if(whichSeq=='Kud'){
          id<-4
        }else{
          stop('Specify which sequence to be found as \'Para\' or \'Mik\' or \'Bay\' or \'Kud\' ')
        }
      }else{
        id <- whichSeq
      }
      aaAlignment<-aaAlignment[c(1,id)]
      id<-2
      ##to calculate identity over overlapping region of scer and para\mik orfs,
      ##this function now returns the only the part overlaps with scer orf
      orf<-checkORF(aaAlignment[[id]])
      if(orf==TRUE){
        coors<-findORF(aaAlignment[[id]])
        startPos<-coors[[1]]
        stopPos<-coors[[2]]

        if(startPos==-1 || stopPos==-1){
          return(FALSE)
        }
        aaSmall<-subseq(aaAlignment,startPos, stopPos)
        aaSmall<-AAStringSet(muscle(aaSmall))
        startPos<-checkCodon(aaSmall[[id]],'M')
        # if(startPos==-1){
        #   startPos<-checkCodon(aaSmall[[id]],'I')
        # }
        stopPos<-checkCodon(aaSmall[[id]],'X')
        aaSmall<-subseq(aaSmall,startPos, stopPos)
        aaSmall<-AAStringSet(muscle(aaSmall))
        return(aaSmall)
      }else{
        if(orf=="BOTH"){
          return(aaAlignment)
        }else if(orf=="M"){
          startPos<-1
          stopPos<-checkCodon(aaAlignment[[id]],'X')

          aaSmall<-subseq(aaAlignment,startPos, stopPos)
          aaSmall<-AAStringSet(muscle(aaSmall))

          stopPos<-checkCodon(aaSmall[[id]],'X')
          aaSmall<-subseq(aaSmall,startPos, stopPos)
          aaSmall<-AAStringSet(muscle(aaSmall))
          return(aaSmall)
        }else if(orf=="X"){
          startPos<-checkCodon(aaAlignment[[id]],'M')
          # if(startPos==-1){
          #   startPos<-checkCodon(aaAlignment[[id]],'I')
          # }
          stopPos<-length(aaAlignment[[1]])

          aaSmall<-subseq(aaAlignment,startPos, stopPos)
          aaSmall<-AAStringSet(muscle(aaSmall))
          startPos<-checkCodon(aaSmall[[id]],'M')
          # if(startPos==-1){
          #   startPos<-checkCodon(aaSmall[[id]],'I')
          # }
          stopPos<-length(aaSmall[[1]])
          aaSmall<-subseq(aaSmall,startPos, stopPos)
          aaSmall<-AAStringSet(muscle(aaSmall))
          return(aaSmall)
        }else if(orf=="Neither"){
          if(is.null(seqAA.alignment)==F){
            if(seqAA.alignment[1]=='M'){
              startPos<-1
              stopPos<-checkCodon(aaAlignment[[id]],'X')

              aaSmall<-subseq(aaAlignment,startPos, stopPos)
              aaSmall<-AAStringSet(muscle(aaSmall))

              stopPos<-checkCodon(aaSmall[[id]],'X')
              aaSmall<-subseq(aaSmall,startPos, stopPos)
              aaSmall<-AAStringSet(muscle(aaSmall))
              return(aaSmall)
            }else{
              startPos<-checkCodon(aaAlignment[[id]],'M')
              # if(startPos==-1){
              #   startPos<-checkCodon(aaAlignment[[id]],'I')
              # }
              stopPos<-length(aaAlignment[[1]])

              aaSmall<-subseq(aaAlignment,startPos, stopPos)
              aaSmall<-AAStringSet(muscle(aaSmall))
              startPos<-checkCodon(aaSmall[[id]],'M')
              # if(startPos==-1){
              #   startPos<-checkCodon(aaSmall[[id]],'I')
              # }
              stopPos<-length(aaSmall[[1]])
              aaSmall<-subseq(aaSmall,startPos, stopPos)
              aaSmall<-AAStringSet(muscle(aaSmall))
              return(aaSmall)
            }
          }else{return(FALSE)}
        }else{return(FALSE)}
      }},
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
