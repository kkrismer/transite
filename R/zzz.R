motif_db <- new.env(parent = emptyenv())

#' @title Retrieve list of all motifs
#'
#' @description
#' Retrieves all Transite motifs
#'
#' @return A list of objects of class Motif
#' @examples
#' transite_motifs <- get_motifs()
#' @family motif functions
#' @importFrom utils data
#' @export
get_motifs <- function() {
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
#' custom_motif <- create_kmer_motif(
#'   "custom_motif", "RBP1",
#'   c("AAAAAAA", "CAAAAAA"), "HITS-CLIP",
#'   "Homo sapiens", "user"
#' )
#' set_motifs(list(custom_motif))
#' @family motif functions
#' @importFrom utils data
#' @export
set_motifs <- function(value) {
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
