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
#' custom_motif <- create_matrix_motif(
#'   "custom_motif", "RBP1",
#'   transite:::toy_motif_matrix, "HITS-CLIP",
#'   "Homo sapiens", "user"
#' )
#' @export
create_matrix_motif <- function(id, rbps, matrix, type, species, src) {
    iupac <- generate_iupac_by_matrix(matrix, code = init_iupac_lookup_table())
    length <- nrow(matrix)
    hexamers <- generate_kmers_from_iupac(iupac, 6)
    heptamers <- generate_kmers_from_iupac(iupac, 7)

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
#' custom_motif <- create_kmer_motif(
#'   "custom_motif", "RBP1",
#'   c("AAAAAAA", "CAAAAAA"), "HITS-CLIP",
#'   "Homo sapiens", "user"
#' )
#' @export
create_kmer_motif <- function(id, rbps, kmers, type, species, src) {
    kmers <- toupper(kmers)
    kmers <- gsub("T", "U", kmers)
    if (!check_kmers(kmers)) {
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

    iupac <- generate_iupac_by_kmers(kmers, code = init_iupac_lookup_table())
    hexamers <- generate_kmers_from_iupac(iupac, 6)
    heptamers <- generate_kmers_from_iupac(iupac, 7)

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
#' get_motifs_meta_info()
#' @family motif functions
#' @export
get_motifs_meta_info <- function() {
    motifs <- get_motifs()
    id <- vapply(motifs, function(motif) {
        return(get_id(motif))
    }, FUN.VALUE = character(1))
    rbps <- vapply(motifs, function(motif) {
        return(paste0(get_rbps(motif), collapse = ", "))
    }, FUN.VALUE = character(1))
    length <- vapply(motifs, function(motif) {
        return(get_width(motif))
    }, FUN.VALUE = integer(1))
    iupac <- vapply(motifs, function(motif) {
        return(get_iupac(motif))
    }, FUN.VALUE = character(1))
    type <- vapply(motifs, function(motif) {
        return(get_type(motif))
    }, FUN.VALUE = character(1))
    species <- vapply(motifs, function(motif) {
        return(get_species(motif))
    }, FUN.VALUE = character(1))
    src <- vapply(motifs, function(motif) {
        return(get_source(motif))
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
#' get_motif_by_id("M178_0.6")
#'
#' get_motif_by_id(c("M178_0.6", "M188_0.6"))
#' @family motif functions
#' @export
get_motif_by_id <- function(id) {
    motifs <- get_motifs()
    idx <- unlist(lapply(motifs, function(motif) {
        return(sum(tolower(get_id(motif)) %in% tolower(id)) > 0)
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
#' get_motif_by_rbp("ELAVL1")
#'
#' get_motif_by_rbp(c("ELAVL1", "ELAVL2"))
#' @family motif functions
#' @export
get_motif_by_rbp <- function(rbp) {
    motifs <- get_motifs()
    idx <- unlist(lapply(motifs, function(motif) {
        return(sum(tolower(get_rbps(motif)) %in% tolower(rbp)) > 0)
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
#' get_ppm(get_motif_by_id("M178_0.6")[[1]])
#' @family motif functions
#' @export
get_ppm <- function(motif) {
    return(2^get_motif_matrix(motif) * 0.25)
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
#' \code{\link{init_iupac_lookup_table}}, it can be specified here
#' @return the IUPAC string of the binding site
#' @examples
#' generate_iupac_by_matrix(get_motif_matrix(get_motif_by_id("M178_0.6")[[1]]))
#' @family motif functions
#' @references \url{http://www.chem.qmul.ac.uk/iubmb/misc/naseq.html}
#' @export
generate_iupac_by_matrix <- function(matrix, threshold = 0.215, code = NULL) {
    if (is.null(code)) {
        code <- init_iupac_lookup_table()
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
#' \code{\link{init_iupac_lookup_table}}, it can be
#' specified here
#' @return the IUPAC string of the binding site
#' @examples
#' generate_iupac_by_kmers(c("AACCAA", "AACCGG", "CACCGA"))
#' @family motif functions
#' @references \url{http://www.chem.qmul.ac.uk/iubmb/misc/naseq.html}
#' @export
generate_iupac_by_kmers <- function(kmers, code = NULL) {
    if (is.null(code)) {
        code <- init_iupac_lookup_table()
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
#' \code{\link{generate_iupac_by_matrix}} function.
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
#' generate_iupac_by_matrix(get_motif_matrix(get_motif_by_id("M178_0.6")[[1]]),
#'   code = init_iupac_lookup_table())
#' @family motif functions
#' @references \url{http://www.chem.qmul.ac.uk/iubmb/misc/naseq.html}
#' @export
init_iupac_lookup_table <- function() {
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
#' generate_kmers_from_iupac(get_iupac(get_motif_by_id("M178_0.6")[[1]]), k = 6)
#' @family motif functions
#' @importFrom Biostrings IUPAC_CODE_MAP
#' @references \url{http://www.chem.qmul.ac.uk/iubmb/misc/naseq.html}
#' @export
generate_kmers_from_iupac <- function(iupac, k) {
    if (nchar(iupac) < k) {
        diff <- k - nchar(iupac)
        left_padding_iupac <- paste0(paste(rep("N", diff), collapse = ""), iupac)
        right_padding_iupac <- paste0(iupac, paste(rep("N", diff), collapse = ""))

        kmers_left_padded <- generate_kmers_from_iupac(left_padding_iupac, k)
        kmers_right_padded <- generate_kmers_from_iupac(right_padding_iupac, k)
        return(unique(c(kmers_left_padded, kmers_right_padded)))
    } else {
        iupac <- gsub(pattern = "U", replacement = "T", x = iupac, fixed = TRUE)
        iupac_vector <- strsplit(iupac, "")[[1]]
        expanded <- lapply(iupac_vector, function(letter) {
            strsplit(Biostrings::IUPAC_CODE_MAP[letter], "")[[1]]
        })

        n <- nchar(iupac)
        full_length <- do.call(paste0, expand.grid(expanded,
                                                   stringsAsFactors = FALSE))
        kmers <- unique(as.vector(apply(as.matrix(full_length), 1, function(x) {
            substring(x, seq_len(n - k + 1), k:n)
        })))
        return(gsub(pattern = "T", replacement = "U", x = kmers, fixed = TRUE))
    }
}

lock <- function(lock_file) {
    if (file.exists(lock_file)) {
        # t <- 3600
        # check whether file is older than time threshold t, if it is,
        # remove it not
        # possible: could cause resource conflict as well
        # con <- file(lock_file)
        # timestamp <- scan(file = con, what = numeric(0), quiet = TRUE)
        # close(con)
        # if(as.numeric(Sys.time()) - timestamp > t ||
        # as.numeric(Sys.time()) - timestamp
        # < 0) { file.remove(lock_file) return(lock(lock_file)) } else {
        return(get_lock_object(lock_file, FALSE))
        # }
    } else {
        if (!dir.exists(dirname(lock_file))) {
            dir.create(dirname(lock_file), showWarnings = FALSE,
                       recursive = TRUE)
        }
        if (file.create(lock_file)) {
            con <- file(lock_file)
            write(as.numeric(Sys.time()), con)
            close(con)
            return(get_lock_object(lock_file, TRUE))
        } else {
            return(get_lock_object(lock_file, FALSE))
        }
    }
}

get_lock_object <- function(lock_file, suc) {
    lock_obj <- list(path = lock_file, success = suc)
    class(lock_obj) <- append(class(lock_obj), "lock_object")
    return(lock_obj)
}

#' @importFrom methods is
unlock <- function(lock_obj) {
    if (!is.null(lock_obj)) {
        stopifnot(methods::is(lock_obj, "lock_object"))
        if (lock_obj$success) {
            if (file.exists(lock_obj$path)) {
                return(file.remove(lock_obj$path))
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
