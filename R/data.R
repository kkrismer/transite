#' Transite motifs
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

#' A toy data set with gene expression measurements and 3'-UTR sequences
#' of 1000 genes
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


#' A toy motif matrix.
#'
#' @format A data frame with four columns (A, C, G, U) and seven
#' rows (position 1 - 7)
"toy.motif.matrix"
