#'This function check a set of strings to create a species vector
#'@param input a list of strings
#'@return a vector that represents whether following species exists in this order. Scer,Spar,Smik,Skud,Sbay,Sarb
#'
#'
#'@export


speciesVector <- function(input){
  vec <- c(1,0,0,0,0,0)
  for(i in 1:length(input)){
    if (grepl("par", input[i], ignore.case = TRUE)) {
      #sparname <- fname
      vec[2] <- 1
    } else if (grepl("mik", input[i], ignore.case = TRUE)) {
      vec[3] <- 1

    } else if (grepl("kud", input[i], ignore.case = TRUE)) {
      vec[4] <- 1

    } else if (grepl("bay", input[i], ignore.case = TRUE)) {
      vec[5] <- 1

    }else if (grepl("arb", input[i], ignore.case = TRUE)) {
      vec[6] <- 1

    }
  }
  vec
}
