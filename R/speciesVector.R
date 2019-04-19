#'This function check a set of strings to create a species vector
#'@param input a list of strings representing the species names in the input file
#'@param specNames species names given by species tree input
#'@return a vector that represents which specNames exists in the alignment file. DOES NOT CONSIDER THE ORDER
#'
#'
#'@export


speciesVector <- function(input,specNames){
  vec <-numeric(length(specNames))# c(1,0,0,0,0,0)
  for(i in 1:length(specNames)){
    if(str_detect(specNames[i],input)%>%any){
      vec[i] <- 1
    }
  }
  vec
}
