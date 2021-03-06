% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/align.R
\name{align}
\alias{align}
\title{Gets DNAStringSet aligns it}
\usage{
align(mySequences, orfName, outputDirectory = NULL, ygeneSeq = NULL,
  algorithm = "Muscle")
}
\arguments{
\item{mySequences}{DNAStringSet with all sequences in it}

\item{orfName}{ORF identifier name that will be used for file writing}

\item{outputDirectory}{path for files to be written in}

\item{ygeneSeq}{sequence of ORF of interest. If not given, last sequence of {mySequences} considered to be ygeneSeq}

\item{algorithm}{alignment algorithm to be used. Default is 'Muscle' but 'ClustalW' and 'ClustalOmega' are available}
}
\value{
aligned DNAStringset, start and stop positions of orfName on DNAStringSet and AA alignment of the orfName sequence with other species
}
\description{
Gets DNAStringSet aligns it
}
