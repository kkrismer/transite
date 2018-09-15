package.environment <- new.env(parent = emptyenv())

#' transite
#'
#' transite is a computational method that allows comprehensive analysis of the regulatory
#' role of RNA-binding proteins in various cellular processes by leveraging preexisting
#' gene expression data and current knowledge of binding preferences of RNA-binding proteins.
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
#' transite.motifs <- getMotifs()
#' @family motif functions
#' @importFrom utils data
#' @export
getMotifs <- function() {
  return(package.environment$motifs)
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
#' custom.motif <- createKmerMotif(
#'   "custom.motif", "RBP1",
#'   c("AAAAAA", "CAAAAA"), "HITS-CLIP",
#'   "Homo sapiens", "user"
#' )
#' setMotifs(list(custom.motif))
#' @family motif functions
#' @importFrom utils data
#' @export
setMotifs <- function(value) {
  old <- package.environment$motifs
  package.environment$motifs <- value
  invisible(old)
}

.onLoad <- function(libname = find.package("transite"), pkgname = "transite") {
  utils::data("motifs", package = pkgname, envir = package.environment)
}

# .onUnload <- function(libpath) {
#  library.dynam.unload("transite", libpath)
# }
