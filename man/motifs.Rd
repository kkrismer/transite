% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{motifs}
\alias{motifs}
\title{Transite Motif Database}
\format{
A list of lists with the following components:
\tabular{rl}{
  \code{id} \tab motif id\cr
  \code{rbps} \tab gene symbols of RNA-binding proteins
  associated with motif\cr
  \code{matrix} \tab data frame of sequence motif
  (position weight matrix)\cr
  \code{hexamers} \tab all motif-associated hexamers\cr
  \code{heptamers} \tab all motif-associated heptamers\cr
  \code{length} \tab length of motif in nucleotides\cr
  \code{iupac} \tab IUPAC string of sequence motif\cr
  \code{type} \tab type of motif, e.g., RNAcompete\cr
  \code{species} \tab usually human\cr
  \code{src} \tab source of motif, e.g., RNA Zoo
}
}
\usage{
data(motifs)
}
\description{
The Transite motif database contains sequence motifs and associated
\emph{k}-mers of more than 100 different RNA-binding proteins, obtained from
publicly available motif databases.
}
\references{
\url{http://cisbp-rna.ccbr.utoronto.ca/}

\url{http://rbpdb.ccbr.utoronto.ca/}
}
\keyword{datasets}
