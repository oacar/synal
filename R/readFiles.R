
#'Read syntetic block sequence files
#'@param path folder path that has sequence fasta files in it
#' every sequence file needs to contain its species name . (i.e. S. Cerevisiae file should contain 'scer', proto-gene sequence should contain 'smorf' or 'sequence')
#' also assumes there is only 1 file for every species
#'@return DNAStringSet object with the sequences in it 
#'order of DNAStringSet : Scer(1) -> Spar(2) -> Smik(3) -> Skud(4) --> Sbay(5) --> Smorf(6) 
#'if one of them does not exist, shows warning 
#'
#'TODO check for multiple files having the search string in them
#'@export



readFiles<-function(path){
  if (dir.exists(path) != TRUE) {
    stop("Path does not exist. Give a correct path")
  }
  scername <- ""
  sparname <- ""
  smikname <- ""
  sbayname <- ""
  skudname <- ""
  smorfname <- ""
  spar <- TRUE
  smik <- TRUE
  sbay <- TRUE
  skud <- TRUE
  for (fname in list.files(path)) {
    if (grepl("Scer", fname, ignore.case = TRUE)) {
      scername <- fname
    } else if (grepl("Spar", fname, ignore.case = TRUE)) {
      sparname <- fname
    } else if (grepl("Smik", fname, ignore.case = TRUE)) {
      smikname <- fname
    } else if (grepl("Sbay", fname, ignore.case = TRUE)) {
      sbayname <- fname
    } else if (grepl("Skud", fname, ignore.case = TRUE)) {
      skudname <- fname
    } else if (grepl("smorf", fname, ignore.case = TRUE) || grepl("sequence", fname, ignore.case = FALSE)){
      
      smorfname <- fname
    }
  }
  if (scername == "") {
    warning("No Scer file is found")
    return(FALSE)
  } else {
    scer_seq <- (readDNAStringSet(sprintf("%s", paste(path, scername, sep = "/"))))
  }
  if (sparname == "") {
    warning("No Spar file is found")
    spar <- FALSE
  } else {
    spar_seq <- (readDNAStringSet(sprintf("%s", paste(path, sparname, sep = "/"))))
  }
  if (smikname == "") {
    warning("No Smik file is found")
    smik <- FALSE
  } else {
    smik_seq <- (readDNAStringSet(sprintf("%s", paste(path, smikname, sep = "/"))))
  }
  if (sbayname == "") {
    warning("No Sbay file is found")
    sbay <- FALSE
  } else {
    sbay_seq <- (readDNAStringSet(sprintf("%s", paste(path, sbayname, sep = "/"))))
  }
  if (skudname == "") {
    warning("No Skud file is found")
    skud <- FALSE
  } else {
    skud_seq <- (readDNAStringSet(sprintf("%s", paste(path, skudname, sep = "/"))))
  }
  if (smorfname == "") {
    stop("No smorf file is found")
  } else {
    smorf <- readDNAStringSet(sprintf("%s", paste(path, smorfname, sep = "/")))
  }
  # read all sequences
  mySequences <- DNAStringSet()
  mySequences <- append(mySequences, (scer_seq))
  if (spar) {
    mySequences <- append(mySequences, (spar_seq))
  }
  if (smik) {
    mySequences <- append(mySequences, (smik_seq))
  }
  if (skud) {
    mySequences <- append(mySequences, (skud_seq))
  }
  if (sbay) {
    mySequences <- append(mySequences, (sbay_seq))
  }
  
  mySequences <- append(mySequences, smorf)
  
  return(mySequences)}