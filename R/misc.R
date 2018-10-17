
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
#' (see \code{\link{generateIUPACByMatrix}})
#' @slot type type of motif (e.g., \code{'HITS-CLIP'}, \code{'EMSA'},
#' \code{'SELEX'}, etc.)
#' @slot species species where motif was discovered (e.g.,
#' \code{'Homo sapiens'})
#' @slot src source of motif (e.g., \code{'RBPDB v1.3.1'})
#' @return Object of type RBPMotif
#' @examples
#' kmers <- c("AAAAAAA", "CAAAAAA")
#' iupac <- generateIUPACByKmers(kmers,
#'   code = initIUPAClookupTable())
#' hexamers <- generateKmersFromIUPAC(iupac, 6)
#' heptamers <- generateKmersFromIUPAC(iupac, 7)
#' new("RBPMotif", id = "custom.motif", rbps = "RBP1",
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

#' Getter Method motifId
#' @name RBPMotif-class
#' @rdname RBPMotif-class
#' @param object RBPMotif object
#' @importFrom methods setGeneric
#' @exportMethod motifId
setGeneric("motifId", function(object) standardGeneric("motifId"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname RBPMotif-class
#' @aliases motifId
#' @aliases motifId,RBPMotif-method
setMethod("motifId", signature(object = "RBPMotif"),
          function(object) object@id)

#' Getter Method motifRbps
#' @name RBPMotif-class
#' @rdname RBPMotif-class
#' @importFrom methods setGeneric
#' @exportMethod motifRbps
setGeneric("motifRbps", function(object) standardGeneric("motifRbps"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname RBPMotif-class
#' @aliases motifRbps
#' @aliases motifRbps,RBPMotif-method
setMethod("motifRbps", signature(object = "RBPMotif"),
          function(object) object@rbps)

#' Getter Method motifMatrix
#' @name RBPMotif-class
#' @rdname RBPMotif-class
#' @importFrom methods setGeneric
#' @exportMethod motifMatrix
setGeneric("motifMatrix", function(object) standardGeneric("motifMatrix"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname RBPMotif-class
#' @aliases motifMatrix
#' @aliases motifMatrix,RBPMotif-method
setMethod("motifMatrix", signature(object = "RBPMotif"),
          function(object) object@matrix)

#' Getter Method motifHexamers
#' @name RBPMotif-class
#' @rdname RBPMotif-class
#' @importFrom methods setGeneric
#' @exportMethod motifHexamers
setGeneric("motifHexamers", function(object) standardGeneric("motifHexamers"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname RBPMotif-class
#' @aliases motifHexamers
#' @aliases motifHexamers,RBPMotif-method
setMethod("motifHexamers", signature(object = "RBPMotif"),
          function(object) object@hexamers)

#' Getter Method motifHeptamers
#' @name RBPMotif-class
#' @rdname RBPMotif-class
#' @importFrom methods setGeneric
#' @exportMethod motifHeptamers
setGeneric("motifHeptamers",
           function(object) standardGeneric("motifHeptamers"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname RBPMotif-class
#' @aliases motifHeptamers
#' @aliases motifHeptamers,RBPMotif-method
setMethod("motifHeptamers", signature(object = "RBPMotif"),
          function(object) object@heptamers)

#' Getter Method motifLength
#' @name RBPMotif-class
#' @rdname RBPMotif-class
#' @importFrom methods setGeneric
#' @exportMethod motifLength
setGeneric("motifLength",
           function(object) standardGeneric("motifLength"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname RBPMotif-class
#' @aliases motifLength
#' @aliases motifLength,RBPMotif-method
setMethod("motifLength", signature(object = "RBPMotif"),
          function(object) object@length)

#' Getter Method motifIUPAC
#' @name RBPMotif-class
#' @rdname RBPMotif-class
#' @importFrom methods setGeneric
#' @exportMethod motifIUPAC
setGeneric("motifIUPAC", function(object) standardGeneric("motifIUPAC"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname RBPMotif-class
#' @aliases motifIUPAC
#' @aliases motifIUPAC,RBPMotif-method
setMethod("motifIUPAC", signature(object = "RBPMotif"),
          function(object) object@iupac)

#' Getter Method motifType
#' @name RBPMotif-class
#' @rdname RBPMotif-class
#' @importFrom methods setGeneric
#' @exportMethod motifType
setGeneric("motifType", function(object) standardGeneric("motifType"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname RBPMotif-class
#' @aliases motifType
#' @aliases motifType,RBPMotif-method
setMethod("motifType", signature(object = "RBPMotif"),
          function(object) object@type)

#' Getter Method motifSpecies
#' @name RBPMotif-class
#' @rdname RBPMotif-class
#' @importFrom methods setGeneric
#' @exportMethod motifSpecies
setGeneric("motifSpecies", function(object) standardGeneric("motifSpecies"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname RBPMotif-class
#' @aliases motifSpecies
#' @aliases motifSpecies,RBPMotif-method
setMethod("motifSpecies", signature(object = "RBPMotif"),
          function(object) object@species)

#' Getter Method motifSource
#' @name RBPMotif-class
#' @rdname RBPMotif-class
#' @importFrom methods setGeneric
#' @exportMethod motifSource
setGeneric("motifSource", function(object) standardGeneric("motifSource"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname RBPMotif-class
#' @aliases motifSource
#' @aliases motifSource,RBPMotif-method
setMethod("motifSource", signature(object = "RBPMotif"),
          function(object) object@src)

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @importFrom methods is
#' @rdname RBPMotif-class
#' @aliases show,RBPMotif-method
setMethod("show", signature(object = "RBPMotif"), function(object) {
    cat(is(object)[[1]], "\n",
        "  id: ", object@id, "\n",
        "  RBPs:  ", paste(object@rbps, collapse = ", "), "\n",
        "  length (in nt):  ", object@length, "\n",
        "  IUPAC:  ", object@iupac, "\n",
        "  type:  ", object@type, "\n",
        "  species:  ", object@species, "\n",
        "  source:  ", object@src, "\n",
        sep = ""
    )
})

#' @param x RBPMotif object
#' @importFrom methods setMethod
#' @importFrom methods signature
#' @importFrom ggseqlogo ggseqlogo
#' @rdname RBPMotif-class
#' @aliases plot,RBPMotif-method
#' @exportMethod plot
setMethod("plot", signature(x = "RBPMotif"), function(x) {
    ppm <- t(getPPM(x))
    return(ggseqlogo::ggseqlogo(ppm))
})

#' @title Creates Transite motif object from position weight matrix
#'
#' @description
#' Takes a position weight matrix (PWM) and meta info and returns an object of
#' class \code{RBPMotif}.
#'
#' @param id motif id (character vector of length 1)
#' @param rbps character vector of names of RNA-binding proteins associated
#' with this motif
#' @param matrix data frame with four columns (A, C, G, U) and 6 - 15 rows
#' (positions),
#' where cell (i, j) contains weight of nucleotide j on position i
#' @param type type of motif (e.g., \code{'HITS-CLIP'}, \code{'EMSA'},
#' \code{'SELEX'}, etc.)
#' @param species species where motif was discovered (e.g.,
#' \code{'Homo sapiens'})
#' @param src source of motif (e.g., \code{'RBPDB v1.3.1'})
#' @return object of class \code{RBPMotif}
#' @examples
#' custom.motif <- createMatrixMotif(
#'   "custom.motif", "RBP1",
#'   transite:::toy.motif.matrix, "HITS-CLIP",
#'   "Homo sapiens", "user"
#' )
#' @export
createMatrixMotif <- function(id, rbps, matrix, type, species, src) {
    iupac <- generateIUPACByMatrix(matrix, code = initIUPAClookupTable())
    length <- nrow(matrix)
    hexamers <- generateKmersFromIUPAC(iupac, 6)
    heptamers <- generateKmersFromIUPAC(iupac, 7)

    return(.RBPMotif(
        id = id, rbps = rbps, matrix = matrix,
        hexamers = hexamers, heptamers = heptamers,
        length = as.integer(length),
        iupac = iupac, type = type, species = species, src = src
    ))
}

#' @title Creates Transite motif object from character vector of \emph{k}-mers
#'
#' @description
#' Takes a position weight matrix (PWM) and meta info and returns an object of
#' class \code{RBPMotif}.
#'
#' @param id motif id (character vector of length 1)
#' @param rbps character vector of names of RNA-binding proteins associated
#' with this motif
#' @param kmers character vector of \emph{k}-mers that are associated with
#' the motif, set of
#' \emph{k}-mers is
#' valid if
#' (1) all \emph{k}-mers must have the same length, (2) only hexamers or
#' heptamers
#' allowed, (3) allowed characters are A, C, G, U
#' @param type type of motif (e.g., \code{'HITS-CLIP'}, \code{'EMSA'},
#' \code{'SELEX'}, etc.)
#' @param species species where motif was discovered (e.g.,
#' \code{'Homo sapiens'})
#' @param src source of motif (e.g., \code{'RBPDB v1.3.1'})
#' @return object of class \code{RBPMotif}
#' @examples
#' custom.motif <- createKmerMotif(
#'   "custom.motif", "RBP1",
#'   c("AAAAAAA", "CAAAAAA"), "HITS-CLIP",
#'   "Homo sapiens", "user"
#' )
#' @export
createKmerMotif <- function(id, rbps, kmers, type, species, src) {
    kmers <- toupper(kmers)
    kmers <- gsub("T", "U", kmers)
    if (!checkKmers(kmers)) {
        stop("invalid k-mers: (1) all k-mers must have the same length, (2) only hexamers or heptamers allowed, (3) allowed characters are A, C, G, U")
    }

    length <- unique(nchar(kmers))
    if (length == 6) {
        hexamers <- kmers
        heptamers <- ""
    } else if (length == 7) {
        hexamers <- ""
        heptamers <- kmers
    }

    iupac <- generateIUPACByKmers(kmers, code = initIUPAClookupTable())
    hexamers <- generateKmersFromIUPAC(iupac, 6)
    heptamers <- generateKmersFromIUPAC(iupac, 7)

    return(.RBPMotif(
        id = id, rbps = rbps, matrix = NULL,
        hexamers = hexamers, heptamers = heptamers,
        length = as.integer(length), iupac = iupac, type = type,
        species = species, src = src
    ))
}

#' @title Displays motif meta information.
#'
#' @description
#' Generates a data frame with meta information about all Transite motifs.
#'
#' @return A data frame containing meta information for all Transite
#' motifs, with the
#' following columns:
#' \itemize{
#'   \item \code{id}
#'   \item \code{rbps}
#'   \item \code{length}
#'   \item \code{iupac}
#'   \item \code{type}
#'   \item \code{species}
#'   \item \code{src}
#' }
#' @examples
#' motifsMetaInfo()
#' @family motif functions
#' @export
motifsMetaInfo <- function() {
    motifs <- getMotifs()
    id <- vapply(motifs, function(motif) {
        return(motifId(motif))
    }, FUN.VALUE = character(1))
    rbps <- vapply(motifs, function(motif) {
        return(paste0(motifRbps(motif), collapse = ", "))
    }, FUN.VALUE = character(1))
    length <- vapply(motifs, function(motif) {
        return(motifLength(motif))
    }, FUN.VALUE = integer(1))
    iupac <- vapply(motifs, function(motif) {
        return(motifIUPAC(motif))
    }, FUN.VALUE = character(1))
    type <- vapply(motifs, function(motif) {
        return(motifType(motif))
    }, FUN.VALUE = character(1))
    species <- vapply(motifs, function(motif) {
        return(motifSpecies(motif))
    }, FUN.VALUE = character(1))
    src <- vapply(motifs, function(motif) {
        return(motifSource(motif))
    }, FUN.VALUE = character(1))
    return(data.frame(id, rbps, length, iupac, type, species, src))
}

#' @title Retrieve motif objects by id
#'
#' @description
#' Retrieves one or more motif objects identified by motif id.
#'
#' @param id character vector of motif identifiers
#' @return A list of objects of class \code{RBPMotif}
#' @examples
#' getMotifById("M178_0.6")
#'
#' getMotifById(c("M178_0.6", "M188_0.6"))
#' @family motif functions
#' @export
getMotifById <- function(id) {
    motifs <- getMotifs()
    idx <- unlist(lapply(motifs, function(motif) {
        return(sum(tolower(motifId(motif)) %in% tolower(id)) > 0)
    }))
    return(motifs[idx])
}

#' @title Retrieve motif objects by gene symbol
#'
#' @description
#' Retrieves one or more motif objects identified by gene symbol.
#'
#' @param rbp character vector of gene symbols of RNA-binding proteins
#' @return A list of objects of class \code{RBPMotif}
#' @examples
#' getMotifByRBP("ELAVL1")
#'
#' getMotifByRBP(c("ELAVL1", "ELAVL2"))
#' @family motif functions
#' @export
getMotifByRBP <- function(rbp) {
    motifs <- getMotifs()
    idx <- unlist(lapply(motifs, function(motif) {
        return(sum(tolower(motifRbps(motif)) %in% tolower(rbp)) > 0)
    }))
    return(motifs[idx])
}

#' @title Get Position Probability Matrix (PPM) from motif object
#'
#' @description
#' Return the position probability matrix of the specified motif.
#'
#' @param motif object of class \code{RBPMotif}
#' @return The position probability matrix of the specified motif
#' @examples
#' getPPM(getMotifById("M178_0.6")[[1]])
#' @family motif functions
#' @export
getPPM <- function(motif) {
    return(2^motifMatrix(motif) * 0.25)
}

#' @title Generates IUPAC code for motif matrix
#'
#' @description
#' Generates a compact logo of a motif based on IUPAC codes given by a
#' position weight matrix
#'
#' @details
#' IUPAC RNA nucleotide code:
#' \tabular{rl}{
#'   \code{A} \tab Adenine\cr
#'   \code{C} \tab Cytosine\cr
#'   \code{G} \tab Guanine\cr
#'   \code{U} \tab Uracil\cr
#'   \code{R} \tab A or G\cr
#'   \code{Y} \tab C or U\cr
#'   \code{S} \tab G or C\cr
#'   \code{W} \tab A or U\cr
#'   \code{K} \tab G or U\cr
#'   \code{M} \tab A or C\cr
#'   \code{B} \tab C or G or U\cr
#'   \code{D} \tab A or G or U\cr
#'   \code{H} \tab A or C or U\cr
#'   \code{V} \tab A or C or G\cr
#'   \code{N} \tab any base
#' }
#' @param matrix the position probability matrix of an RNA-binding protein
#' @param threshold the threshold probability (nucleotides with lower
#' probabilities are ignored)
#' @param code if IUPAC code table has already been initialized by
#' \code{\link{initIUPAClookupTable}}, it can be specified here
#' @return the IUPAC string of the binding site
#' @examples
#' generateIUPACByMatrix(motifMatrix(getMotifById("M178_0.6")[[1]]))
#' @family motif functions
#' @references \url{http://www.chem.qmul.ac.uk/iubmb/misc/naseq.html}
#' @export
generateIUPACByMatrix <- function(matrix, threshold = 0.215, code = NULL) {
    if (is.null(code)) {
        code <- initIUPAClookupTable()
    }

    codes <- apply(matrix, 1, function(x) {
        value <- ""
        if (x[1] > threshold) {
            value <- "A"
        }
        if (x[2] > threshold) {
            value <- paste0(value, "C")
        }
        if (x[3] > threshold) {
            value <- paste0(value, "G")
        }
        if (x[4] > threshold) {
            value <- paste0(value, "U")
        }
        return(code[[value]])
    })

    return(paste0(codes, collapse = ""))
}

#' @title Generates IUPAC code for a character vector of \emph{k}-mers
#'
#' @description
#' Generates a compact logo of a motif based on IUPAC codes given by a character
#' vector of \emph{k}-mers
#'
#' @details
#' IUPAC RNA nucleotide code:
#' \tabular{rl}{
#'   \code{A} \tab Adenine\cr
#'   \code{C} \tab Cytosine\cr
#'   \code{G} \tab Guanine\cr
#'   \code{U} \tab Uracil\cr
#'   \code{R} \tab A or G\cr
#'   \code{Y} \tab C or U\cr
#'   \code{S} \tab G or C\cr
#'   \code{W} \tab A or U\cr
#'   \code{K} \tab G or U\cr
#'   \code{M} \tab A or C\cr
#'   \code{B} \tab C or G or U\cr
#'   \code{D} \tab A or G or U\cr
#'   \code{H} \tab A or C or U\cr
#'   \code{V} \tab A or C or G\cr
#'   \code{N} \tab any base
#' }
#' @param kmers character vector of \emph{k}-mers
#' @param code if IUPAC code table has already been initialized by
#' \code{\link{initIUPAClookupTable}}, it can be
#' specified here
#' @return the IUPAC string of the binding site
#' @examples
#' generateIUPACByKmers(c("AACCAA", "AACCGG", "CACCGA"))
#' @family motif functions
#' @references \url{http://www.chem.qmul.ac.uk/iubmb/misc/naseq.html}
#' @export
generateIUPACByKmers <- function(kmers, code = NULL) {
    if (is.null(code)) {
        code <- initIUPAClookupTable()
    }

    codes <- lapply(seq_len(nchar(kmers[1])), function(i) {
        return(code[[paste0(sort(unique(unlist(lapply(seq_len(length(kmers)),
                                                      function(j) {
            return(substr(kmers[j], i, i))
        })))), collapse = "")]])
    })

    return(paste0(codes, collapse = ""))
}

#' @title Initializes the IUPAC lookup table
#'
#' @description
#' Initializes a hash table that serves as a IUPAC lookup table for the
#' \code{\link{generateIUPACByMatrix}} function.
#'
#' @details
#' IUPAC RNA nucleotide code:
#' \tabular{rl}{
#'   \code{A} \tab Adenine\cr
#'   \code{C} \tab Cytosine\cr
#'   \code{G} \tab Guanine\cr
#'   \code{U} \tab Uracil\cr
#'   \code{R} \tab A or G\cr
#'   \code{Y} \tab C or U\cr
#'   \code{S} \tab G or C\cr
#'   \code{W} \tab A or U\cr
#'   \code{K} \tab G or U\cr
#'   \code{M} \tab A or C\cr
#'   \code{B} \tab C or G or U\cr
#'   \code{D} \tab A or G or U\cr
#'   \code{H} \tab A or C or U\cr
#'   \code{V} \tab A or C or G\cr
#'   \code{N} \tab any base
#' }
#' @return an environment, the IUPAC lookup hash table
#' @examples
#' generateIUPACByMatrix(motifMatrix(getMotifById("M178_0.6")[[1]]),
#'   code = initIUPAClookupTable())
#' @family motif functions
#' @references \url{http://www.chem.qmul.ac.uk/iubmb/misc/naseq.html}
#' @export
initIUPAClookupTable <- function() {
    table <- new.env(hash = TRUE, parent = emptyenv())
    table[["A"]] <- "A"
    table[["C"]] <- "C"
    table[["G"]] <- "G"
    table[["U"]] <- "U"
    table[["AG"]] <- "R"
    table[["CU"]] <- "Y"
    table[["AC"]] <- "M"
    table[["GU"]] <- "K"
    table[["AU"]] <- "W"
    table[["CG"]] <- "S"
    table[["CGU"]] <- "B"
    table[["AGU"]] <- "D"
    table[["ACU"]] <- "H"
    table[["ACG"]] <- "V"
    table[["ACGU"]] <- "N"

    return(table)
}

#' @title Generates all \emph{k}-mers for IUPAC string
#'
#' @description
#' Generates all possible \emph{k}-mers for a given IUPAC string.
#'
#' @details
#' IUPAC RNA nucleotide code:
#' \tabular{rl}{
#'   \code{A} \tab Adenine\cr
#'   \code{C} \tab Cytosine\cr
#'   \code{G} \tab Guanine\cr
#'   \code{U} \tab Uracil\cr
#'   \code{R} \tab A or G\cr
#'   \code{Y} \tab C or U\cr
#'   \code{S} \tab G or C\cr
#'   \code{W} \tab A or U\cr
#'   \code{K} \tab G or U\cr
#'   \code{M} \tab A or C\cr
#'   \code{B} \tab C or G or U\cr
#'   \code{D} \tab A or G or U\cr
#'   \code{H} \tab A or C or U\cr
#'   \code{V} \tab A or C or G\cr
#'   \code{N} \tab any base
#' }
#' @param iupac IUPAC string
#' @param k length of \emph{k}-mer, \code{6} (hexamers) or \code{7} (heptamers)
#' @return list of \emph{k}-mers
#' @examples
#' generateKmersFromIUPAC(motifIUPAC(getMotifById("M178_0.6")[[1]]), k = 6)
#' @family motif functions
#' @importFrom Biostrings IUPAC_CODE_MAP
#' @references \url{http://www.chem.qmul.ac.uk/iubmb/misc/naseq.html}
#' @export
generateKmersFromIUPAC <- function(iupac, k) {
    if (nchar(iupac) < k) {
        diff <- k - nchar(iupac)
        left.padding.iupac <- paste0(paste(rep("N", diff), collapse = ""), iupac)
        right.padding.iupac <- paste0(iupac, paste(rep("N", diff), collapse = ""))

        kmers.left.padded <- generateKmersFromIUPAC(left.padding.iupac, k)
        kmers.right.padded <- generateKmersFromIUPAC(right.padding.iupac, k)
        return(unique(c(kmers.left.padded, kmers.right.padded)))
    } else {
        iupac <- gsub(pattern = "U", replacement = "T", x = iupac, fixed = TRUE)
        iupac.vector <- strsplit(iupac, "")[[1]]
        expanded <- lapply(iupac.vector, function(letter) {
            strsplit(Biostrings::IUPAC_CODE_MAP[letter], "")[[1]]
        })

        n <- nchar(iupac)
        full.length <- do.call(paste0, expand.grid(expanded,
                                                   stringsAsFactors = FALSE))
        kmers <- unique(as.vector(apply(as.matrix(full.length), 1, function(x) {
            substring(x, seq_len(n - k + 1), k:n)
        })))
        return(gsub(pattern = "T", replacement = "U", x = kmers, fixed = TRUE))
    }
}

lock <- function(lock.file) {
    t <- 3600
    if (file.exists(lock.file)) {
        # check whether file is older than time threshold t, if it is,
        # remove it not
        # possible: could cause resource conflict as well
        # con <- file(lock.file)
        # timestamp <- scan(file = con, what = numeric(0), quiet = TRUE)
        # close(con)
        # if(as.numeric(Sys.time()) - timestamp > t ||
        # as.numeric(Sys.time()) - timestamp
        # < 0) { file.remove(lock.file) return(lock(lock.file)) } else {
        return(getLockObject(lock.file, FALSE))
        # }
    } else {
        if (!dir.exists(dirname(lock.file))) {
            dir.create(dirname(lock.file), showWarnings = FALSE,
                       recursive = TRUE)
        }
        if (file.create(lock.file)) {
            con <- file(lock.file)
            write(as.numeric(Sys.time()), con)
            close(con)
            return(getLockObject(lock.file, TRUE))
        } else {
            return(getLockObject(lock.file, FALSE))
        }
    }
}

getLockObject <- function(lock.file, suc) {
    lock.obj <- list(path = lock.file, success = suc)
    class(lock.obj) <- append(class(lock.obj), "lock.object")
    return(lock.obj)
}

#' @importFrom methods is
unlock <- function(lock.obj) {
    if (!is.null(lock.obj)) {
        stopifnot(methods::is(lock.obj, "lock.object"))
        if (lock.obj$success) {
            if (file.exists(lock.obj$path)) {
                return(file.remove(lock.obj$path))
            } else {
                return(TRUE)
            }
        } else {
            return(FALSE)
        }
    } else {
        return(FALSE)
    }
}
