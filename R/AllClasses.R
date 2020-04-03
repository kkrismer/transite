setClassUnion("dfOrNULL", c("data.frame", "NULL"))

#' An S4 class to represent a RBPMotif
#'
#' @slot id motif id (character vector of length 1)
#' @slot rbps character vector of names of RNA-binding proteins associated
#' with this motif
#' @slot matrix data frame with four columns (A, C, G, U) and 6 - 15 rows
#' (positions),
#' where cell (i, j) contains weight of nucleotide j on position i
#' @slot hexamers character vector of hexamers associated with this motif
#' @slot heptamers character vector of heptamers associated with this motif
#' @slot length length of the motif (i.e., \code{nrow(matrix)})
#' @slot iupac IUPAC code for motif matrix
#' (see \code{\link{generate_iupac_by_matrix}})
#' @slot type type of motif (e.g., \code{'HITS-CLIP'}, \code{'EMSA'},
#' \code{'SELEX'}, etc.)
#' @slot species species where motif was discovered (e.g.,
#' \code{'Homo sapiens'})
#' @slot src source of motif (e.g., \code{'RBPDB v1.3.1'})
#' @return Object of type RBPMotif
#' @examples
#' kmers <- c("AAAAAAA", "CAAAAAA")
#' iupac <- generate_iupac_by_kmers(kmers,
#'   code = init_iupac_lookup_table())
#' hexamers <- generate_kmers_from_iupac(iupac, 6)
#' heptamers <- generate_kmers_from_iupac(iupac, 7)
#' new("RBPMotif", id = "custom_motif", rbps = "RBP1",
#'   matrix = NULL, hexamers = hexamers, heptamers = heptamers, length = 7L,
#'   iupac = iupac, type = "HITS-CLIP", species = "Homo sapiens", src = "user"
#' )
#' @importFrom methods setClass
#' @importFrom methods setClassUnion
#' @importFrom methods new
#' @name RBPMotif-class
#' @rdname RBPMotif-class
#' @exportClass RBPMotif
.RBPMotif <- setClass("RBPMotif", slots = c(id = "character",
                                            rbps = "character",
                                            matrix = "dfOrNULL",
                                            hexamers = "character",
                                            heptamers = "character",
                                            length = "integer",
                                            iupac = "character",
                                            type = "character",
                                            species = "character",
                                            src = "character"))

#' An S4 class to represent a scored spectrum
#'
#' @slot adj_r_squared adjusted \eqn{R^2} of polynomial model
#' @slot degree degree of polynomial (integer between 0 and 5)
#' @slot residuals residuals of the polynomial model
#' @slot slope coefficient of the linear term of the polynomial model
#' (spectrum "direction")
#' @slot f_statistic F statistic from the F test used to determine the
#' degree of the polynomial model
#' @slot f_statistic_p_value p-value associated with the F statistic
#' @slot consistency_score raw local consistency score of the spectrum
#' @slot consistency_score_p_value p-value associated with the local consistency
#' score
#' @slot consistency_score_n number of permutations performed to calculate
#' p-value of local consistency score (permutations performed before early
#' stopping criterion reached)
#' @slot plot spectrum plot
#' @examples
#' new("SpectrumScore",
#'     adj_r_squared = 0,
#'     degree = 0L,
#'     residuals = 0,
#'     slope = 0,
#'     f_statistic = 0,
#'     f_statistic_p_value = 1,
#'     consistency_score = 1,
#'     consistency_score_p_value = 1,
#'     consistency_score_n = 1000L,
#'     plot = NULL
#' )
#' @return Object of type SpectrumScore
#' @importFrom methods setClass
#' @name SpectrumScore-class
#' @rdname SpectrumScore-class
#' @exportClass SpectrumScore
.SpectrumScore <- setClass("SpectrumScore", slots = c(adj_r_squared = "numeric",
                                                      degree = "integer",
                                                      residuals = "numeric",
                                                      slope = "numeric",
                                                      f_statistic = "numeric",
                                                      f_statistic_p_value = "numeric",
                                                      consistency_score = "numeric",
                                                      consistency_score_p_value = "numeric",
                                                      consistency_score_n = "integer",
                                                      plot = "ANY"))
