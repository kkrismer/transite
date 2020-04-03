motif_db <- new.env(parent = emptyenv())

#' transite
#'
#' transite is a computational method that allows comprehensive analysis of
#' the regulatory
#' role of RNA-binding proteins in various cellular processes by leveraging
#' preexisting
#' gene expression data and current knowledge of binding preferences of
# RNA-binding proteins.
#'
#' @docType package
#' @author Konstantin Krismer
#' @useDynLib transite
#' @name transite
#' @importFrom Rcpp sourceCpp
NULL

#' @title Retrieve list of all motifs
#'
#' @description
#' Retrieves all Transite motifs
#'
#' @return A list of objects of class Motif
#' @examples
#' transite_motifs <- getMotifs()
#' @family motif functions
#' @importFrom utils data
#' @export
getMotifs <- function() {
    return(motif_db$motifs)
}

#' @title Set Transite motif database
#'
#' @description
#' Globally sets Transite motif database, use with care.
#'
#' @param value list of Motif objects
#'
#' @return void
#'
#' @examples
#' custom_motif <- createKmerMotif(
#'   "custom_motif", "RBP1",
#'   c("AAAAAAA", "CAAAAAA"), "HITS-CLIP",
#'   "Homo sapiens", "user"
#' )
#' setMotifs(list(custom_motif))
#' @family motif functions
#' @importFrom utils data
#' @export
setMotifs <- function(value) {
    old <- motif_db$motifs
    motif_db$motifs <- value
    invisible(old)
}

.onLoad <- function(libname = find.package("transite"), pkgname = "transite") {
    utils::data("motifs", package = pkgname, envir = motif_db)
    envir <- parent.env(environment())
    utils::data("ge", package = pkgname, envir = envir)
    utils::data("toy_motif_matrix", package = pkgname, envir = envir)
    utils::data("kmers_enrichment", package = pkgname, envir = envir)
}
