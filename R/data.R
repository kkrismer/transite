#' Transite Motif Database
#'
#' The Transite motif database contains sequence motifs and associated
#' \emph{k}-mers of more than 100 different RNA-binding proteins, obtained from
#' publicly available motif databases.
#' @references
#' \url{http://cisbp-rna.ccbr.utoronto.ca/}
#'
#' \url{http://rbpdb.ccbr.utoronto.ca/}
#'
#'
#' @format A list of lists with the following components:
#' \tabular{rl}{
#'   \code{id} \tab motif id\cr
#'   \code{rbps} \tab gene symbols of RNA-binding proteins
#'   associated with motif\cr
#'   \code{matrix} \tab data frame of sequence motif
#'   (position weight matrix)\cr
#'   \code{hexamers} \tab all motif-associated hexamers\cr
#'   \code{heptamers} \tab all motif-associated heptamers\cr
#'   \code{length} \tab length of motif in nucleotides\cr
#'   \code{iupac} \tab IUPAC string of sequence motif\cr
#'   \code{type} \tab type of motif, e.g., RNAcompete\cr
#'   \code{species} \tab usually human\cr
#'   \code{src} \tab source of motif, e.g., RNA Zoo
#' }
"motifs"

#' Toy Gene Expression Data Set
#'
#' This object contains a toy data set based on gene expression measurements
#' and 3'-UTR sequences of 1000 genes. It comprises three data frames with
#' RefSeq identifiers, log fold change values, and 3'-UTR sequences of genes,
#' which are either upregulated or downregulated after some hypothetical
#' treatment, as well as all measured genes. The actual values are not
#' important. This data set merely serves as an example input for various
#' functions.
#'
#' @format A list with the following components:
#' \tabular{rl}{
#'   \code{foreground1.df} \tab data frame that contains down-regulated
#'   genes after treatment\cr
#'   \code{foreground2.df} \tab data frame that contains up-regulated
#'   genes after treatment\cr
#'   \code{background.df} \tab data frame that contains all genes measured
#' }
"ge"


#' Toy Motif Matrix
#'
#' This toy motif matrix is used in code examples for various functions.
#'
#' @format A data frame with four columns (A, C, G, U) and seven
#' rows (position 1 - 7)
"toy.motif.matrix"
