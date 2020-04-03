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
#'   \code{foreground1_df} \tab data frame that contains down-regulated
#'   genes after treatment\cr
#'   \code{foreground2_df} \tab data frame that contains up-regulated
#'   genes after treatment\cr
#'   \code{background_df} \tab data frame that contains all genes measured
#' }
"ge"

#' Toy Motif Matrix
#'
#' This toy motif matrix is used in code examples for various functions.
#'
#' @format A data frame with four columns (A, C, G, U) and seven
#' rows (position 1 - 7)
"toy_motif_matrix"

#' Example \emph{k}-mer Enrichment Data
#'
#' This data frame with \emph{k}-mer enrichment data (as produced by
#' \code{\link{run_kmer_tsma}}) is used in a code example for k-mer volcano plot
#' function \code{\link{draw_volcano_plot}}.
#'
#' @format A data frame with the following columns:
#' \tabular{rl}{
#'   \code{kmer} \tab contains all hexamers (AAAAAA to UUUUUU)\cr
#'   \code{foreground_count} \tab absolute \emph{k}-mer frequency in
#'   foreground set\cr
#'   \code{background_count} \tab absolute \emph{k}-mer frequency in
#'   background set\cr
#'   \code{enrichment} \tab enrichment of \emph{k}-mer in foreground
#'   relative to background\cr
#'   \code{p_value} \tab associated p-value of enrichment\cr
#'   \code{adj_p_value} \tab multiple testing corrected p-value\cr
#' }
"kmers_enrichment"
