test_that("find homolog doesn't find length<12nt overlap", {
  DNAStr <- Biostrings::DNAStringSet(c("ATGTGGTGA", "ATGGGGTAG"))
  start <- 1
  stop <- 9
  ygeneSeq <- Biostrings::DNAStringSet(c("ATGTGGTGA"))
  types <- c("a", "b")
  names(DNAStr) <- types
  expect_equal(list(), findHomolog(DNAStr, start, stop, ygeneSeq, types))
})

test_that("find homolog works same length, no gap sequences", {
  DNAStr <- Biostrings::DNAStringSet(c("ATGTGGAAATTTGGGTGA",
                                       "ATGGGAAATTTGGGGTAG"))
  start <- 1
  stop <- 18
  ygeneSeq <- Biostrings::DNAStringSet(c("ATGTGGAAATTTGGGTGA"))
  types <- c("a", "b")
  names(DNAStr) <- types
  tt <- unlist(findHomolog(DNAStr, start, stop, ygeneSeq, types))
  expect_equal(DNAStr, tt$b.dna)
  expect_equal(DNAStr, tt$b.dnaOverlap)
})

test_that("find homolog works with gaps same frame", {
  DNAStr <- Biostrings::DNAStringSet(c("ATGTGGTTTAGCCCTGGGTGA",
                                       "ATGGGA---AATTTGGGGTAG"))
  start <- 1
  stop <- 21
  ygeneSeq <- Biostrings::DNAStringSet(c("ATGTGGTTTAGCCCTGGGTGA"))
  types <- c("a", "b")
  names(DNAStr) <- types
  tt <- unlist(findHomolog(DNAStr, start, stop, ygeneSeq, types))
  expect_equal(DNAStr, tt$b.dna)
  expect_equal(DNAStr, tt$b.dnaOverlap)
})

test_that("find homolog works with gaps different frame", {
  DNAStr <- Biostrings::DNAStringSet(c("ATGTGG--TTTAGCCCTGGGTGA",
                                       "ATG-GGAA--T-G----GGGTAG"))
  start <- 1
  stop <- 23
  ygeneSeq <- Biostrings::DNAStringSet(c("ATGTGGTTTAGCCCTGGGTGA"))
  types <- c("a", "b")
  names(DNAStr) <- types
  tt <- unlist(findHomolog(DNAStr, start, stop, ygeneSeq,types))
  expect_equal(tt$b.dnaOverlap, tt$b.dna)
  expect_equal(tt$b.aa, tt$b.aaOverlap )
})

test_that("find homolog works with gaps different frame-different length orfs", {
  DNAStr <- Biostrings::DNAStringSet(c("TATATGGTTTACCCCTGGTGA",
                                       "AATGGGAAT-----GGGGTAG"))
  start <- 4
  stop <- 21
  ygeneSeq <- Biostrings::DNAStringSet(c("ATGGTTTACCCCTGGTGA"))
  types <- c("a", "b")
  names(DNAStr) <- types
  tt <- unlist(findHomolog(DNAStr, start, stop, ygeneSeq, types))
  expect_equal(tt$b.dnaOverlap, tt$b.dna)
  expect_equal(tt$b.aa, tt$b.aaOverlap )
})


test_that("find homolog doesn't work when overlap consists of only gaps", {
  DNAStr <- Biostrings::DNAStringSet(c("ATGAATATGGTTTACCCCTGGTGATAG",
                                       "ATGAAT------------------TAG"))
  start <- 9
  stop <- 24
  ygeneSeq <- Biostrings::DNAStringSet(c("ATGGTTTACCCCTGGTGA"))
  types <- c("a", "b")
  names(DNAStr) <- types
  tt <- unlist(findHomolog(DNAStr, start, stop, ygeneSeq, types))
  expect_equal(NULL, tt$b.dna)
})
