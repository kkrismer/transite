#' @title Scores transcripts with position weight matrices
#'
#' @description
#' This function is used to count the binding sites in a set of sequences for
#' all or a
#' subset of RNA-binding protein sequence
#' motifs and returns the result in a data frame, which is subsequently used by
#' \code{\link{calculate_motif_enrichment}} to
#' obtain binding site enrichment scores.
#'
#' @param motifs a list of motifs that is used to score the specified sequences.
#' If \code{is.null(motifs)} then all Transite motifs are used.
#' @param n_cores the number of cores that are used
#' @param cache either logical or path to a directory where scores are cached.
#' The scores of each
#' motif are stored in a
#' separate file that contains a hash table with RefSeq identifiers and
#' sequence type
#' qualifiers as keys and the number of putative binding sites as values.
#' If \code{cache} is \code{FALSE}, scores will not be cached.
#' @inheritParams score_transcripts_single_motif
#'
#' @return A list with three entries:
#'
#' (1) df: a data frame with the following columns:
#' \tabular{rl}{
#'   \code{motif_id} \tab the motif identifier that is used in the original
#'   motif library\cr
#'   \code{motif_rbps} \tab the gene symbol of the RNA-binding protein(s)\cr
#'   \code{absolute_hits} \tab the absolute frequency of putative binding
#'   sites per motif in all
#'   transcripts \cr
#'   \code{relative_hits} \tab  the relative, i.e., absolute divided by total,
#'   frequency of
#'   binding sites per motif in all transcripts \cr
#'   \code{total_sites} \tab the total number of potential binding sites \cr
#'   \code{one_hit}, \code{two_hits}, ... \tab number of transcripts with one,
#'   two,
#'   three, ... putative binding sites
#' }
#' (2) total_sites: a numeric vector with the total number of potential
#' binding sites
#' per transcript
#'
#' (3) absolute_hits: a numeric vector with the absolute (not relative)
#' number of putative
#' binding sites per transcript
#' @examples
#' foreground_set <- c(
#'   "CAACAGCCUUAAUU", "CAGUCAAGACUCC", "CUUUGGGGAAU",
#'   "UCAUUUUAUUAAA", "AAUUGGUGUCUGGAUACUUCCCUGUACAU",
#'   "AUCAAAUUA", "AGAU", "GACACUUAAAGAUCCU",
#'   "UAGCAUUAACUUAAUG", "AUGGA", "GAAGAGUGCUCA",
#'   "AUAGAC", "AGUUC", "CCAGUAA"
#' )
#' # names are used as keys in the hash table (cached version only)
#' # ideally sequence identifiers (e.g., RefSeq ids) and region labels
#' # (e.g., 3UTR for 3'-UTR)
#' names(foreground_set) <- c(
#'   "NM_1_DUMMY|3UTR", "NM_2_DUMMY|3UTR", "NM_3_DUMMY|3UTR",
#'   "NM_4_DUMMY|3UTR", "NM_5_DUMMY|3UTR", "NM_6_DUMMY|3UTR",
#'   "NM_7_DUMMY|3UTR", "NM_8_DUMMY|3UTR", "NM_9_DUMMY|3UTR",
#'   "NM_10_DUMMY|3UTR", "NM_11_DUMMY|3UTR", "NM_12_DUMMY|3UTR",
#'   "NM_13_DUMMY|3UTR", "NM_14_DUMMY|3UTR"
#' )
#'
#' # specific motifs, uncached
#' motifs <- get_motif_by_rbp("ELAVL1")
#' scores <- score_transcripts(foreground_set, motifs = motifs, cache = FALSE)
#' \dontrun{
#' # all Transite motifs, cached (writes scores to disk)
#' scores <- score_transcripts(foreground_set)
#'
#' # all Transite motifs, uncached
#' scores <- score_transcripts(foreground_set, cache = FALSE)
#'
#' foreground_df <- transite:::ge$foreground1_df
#' foreground_set <- foreground_df$seq
#' names(foreground_set) <- paste0(foreground_df$refseq, "|",
#'    foreground_df$seq_type)
#' scores <- score_transcripts(foreground_set)
#' }
#' @family matrix functions
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel parLapply
#' @export
score_transcripts <- function(sequences, motifs = NULL, max_hits = 5,
                             threshold_method = "p_value",
                             threshold_value = 0.25^6,
                             n_cores = 1, cache = paste0(tempdir(), "/sc/")) {
    # if threshold_method == "p_value": default threshold_value == 0.25^6
    # (0.25^6: lowest p-value that can be achieved by hexamer motifs, the
    #shortest supported motifs)
    # if threshold_method == "relative": default threshold_value == 0.9
    # (0.9: 90% of the maximum PWM score)
    if (max_hits < 1) {
        stop("max_hits must not be smaller than 1")
    }

    if (is.null(motifs)) {
        motifs <- get_motifs()
    }

    if (is.logical(cache)) {
        if (cache) {
            cache_path <- paste0(tempdir(), "/sc/")
        } else {
            cache_path <- NULL
        }
    } else {
        cache_path <- as.character(cache)
    }

    if (!is.null(cache_path)) {
        if (is.null(names(sequences))) {
            stop("names(sequences) cannot be null if cache is used")
        }
        names(sequences) <- paste0(names(sequences), "|",
                                   threshold_method, "|", threshold_value)
        dir.create(file.path(cache_path), showWarnings = FALSE,
                   recursive = TRUE)
    }

    if (n_cores == 1) {
        motif_scores <- lapply(motifs, function(motif) {
            return(score_transcripts_single_motif(
                motif, sequences, max_hits, threshold_method,
                threshold_value, cache_path
            ))
        })
    } else {
        cluster <- parallel::makeCluster(mc <- getOption("cl.cores", n_cores))
        parallel::clusterExport(
            cl = cluster,
            varlist = c(
                "score_transcripts_single_motif", "sequences", "max_hits",
                "threshold_method", "threshold_value", "cache_path",
                "motifId", "motifRbps", "motifMatrix", "motifLength",
                "read_motif_cache", "write_motif_cache", "lock", "unlock",
                "get_lock_object", "score_sequences",
                "cached_score_sequences_helper", "score_sequences_helper"
            ),
            envir = environment()
        )
        motif_scores <- parallel::parLapply(cl = cluster, motifs,
                                            function(motif) {
            return(score_transcripts_single_motif(
                motif, sequences, max_hits, threshold_method,
                threshold_value, cache_path
            ))
        })
    }

    motif_scores_proto_df <- lapply(motif_scores, function(motif_score) {
        return(motif_score$df)
    })

    total_sites <- lapply(motif_scores, function(motif_score) {
        return(motif_score$total_sites)
    })

    absolute_hits <- lapply(motif_scores, function(motif_score) {
        return(motif_score$absolute_hits)
    })

    motif_scores_df <- as.data.frame(do.call("rbind", motif_scores_proto_df),
                                     stringsAsFactors = FALSE
    )
    motif_scores_df$motif_id <- as.character(motif_scores_df$motif_id)
    motif_scores_df$motif_rbps <- as.character(motif_scores_df$motif_rbps)
    motif_scores_df$absolute_hits <- as.numeric(motif_scores_df$absolute_hits)
    motif_scores_df$relative_hits <- as.numeric(motif_scores_df$relative_hits)
    motif_scores_df$total_sites <- as.numeric(motif_scores_df$total_sites)
    motif_scores_df$one_hit <- as.numeric(motif_scores_df$one_hit)
    if (max_hits > 1) {
        motif_scores_df$two_hits <- as.numeric(motif_scores_df$two_hits)
    }
    if (max_hits > 2) {
        motif_scores_df$three_hits <- as.numeric(motif_scores_df$three_hits)
    }
    if (max_hits > 3) {
        motif_scores_df$four_hits <- as.numeric(motif_scores_df$four_hits)
    }
    if (max_hits > 4) {
        motif_scores_df$five_hits <- as.numeric(motif_scores_df$five_hits)
    }
    if (max_hits > 5) {
        motif_scores_df$six_hits <- as.numeric(motif_scores_df$six_hits)
    }
    if (max_hits > 6) {
        motif_scores_df$seven_hits <- as.numeric(motif_scores_df$seven_hits)
    }
    if (max_hits > 7) {
        motif_scores_df$eight_hit <- as.numeric(motif_scores_df$eight_hit)
    }
    if (max_hits > 8) {
        motif_scores_df$nine_hits <- as.numeric(motif_scores_df$nine_hits)
    }
    motif_scores_df$more_hits <- as.numeric(motif_scores_df$more_hits)
    return(list(df = motif_scores_df, total_sites = total_sites,
                absolute_hits = absolute_hits))
}

#' @title Scores transadsadscripts with position weight matrices
#'
#' @description
#' This function is used to count the putative binding sites (i.e., motifs)
#' in a set of
#' sequences for the specified RNA-binding protein sequence
#' motifs and returns the result in a data frame, which is aggregated by
#' \code{\link{score_transcripts}} and
#' subsequently used by \code{\link{calculate_motif_enrichment}} to
#' obtain binding site enrichment scores.
#'
#' @param motif a Transite motif that is used to score the specified sequences
#' @param sequences character vector of named sequences
#' (only containing upper case characters A, C, G, T), where the names are
#' RefSeq identifiers
#' and sequence
#' type qualifiers (\code{"3UTR"}, \code{"5UTR"}, \code{"mRNA"}), e.g.
#' \code{"NM_010356|3UTR"}
#' @param max_hits maximum number of putative binding sites per mRNA
#' that are counted
#' @param cache_path the path to a directory where scores are cached.
#' The scores of each
#' motif are stored in a
#' separate file that contains a hash table with RefSeq identifiers
#' and sequence type
#' qualifiers as keys and the number of binding sites as values.
#' If is.null(cache_path), scores will not be cached.
#' @param threshold_method either \code{"p_value"} (default) or
#' \code{"relative"}.
#' If \code{threshold_method} equals \code{"p_value"}, the default
#' \code{threshold_value}
#'  is \code{0.25^6}, which is
#' lowest p-value that can be achieved by hexamer motifs, the shortest
#' supported motifs.
#' If \code{threshold_method} equals \code{"relative"}, the default
#' \code{threshold_value}
#' is \code{0.9}, which is 90\% of the maximum PWM score.
#' @param threshold_value semantics of the \code{threshold_value} depend on
#' \code{threshold_method} (default is 0.25^6)
#' @return A list with the following items:
#' \tabular{rl}{
#'   \code{motif_id} \tab the motif identifier of the specified motif\cr
#'   \code{motif_rbps} \tab the gene symbol of the RNA-binding protein(s)\cr
#'   \code{absolute_hits} \tab the absolute frequency of binding sites per
#'   motif in all
#'   transcripts \cr
#'   \code{relative_hits} \tab  the relative, i.e., absolute divided by
#'   total, frequency of
#'   binding sites per motif in all transcripts \cr
#'   \code{total_sites} \tab the total number of potential binding sites \cr
#'   \code{one_hit}, \code{two_hits}, ... \tab number of transcripts with
#'   one, two, three,
#'   ... binding sites
#' }
#'
#' @family matrix functions
#' @importFrom TFMPvalue TFMpv2sc
#' @importFrom Biostrings maxScore
score_transcripts_single_motif <- function(motif, sequences, max_hits = 5,
                                        threshold_method = "p_value",
                                        threshold_value = 0.25^6,
                                        cache_path = paste0(tempdir(), "/sc/")) {
    # if threshold_method == "p_value": default threshold_value == 0.25^6
    # (0.25^6: lowest p-value that can be achieved by hexamer motifs,
    # the shortest supported motifs)
    # if threshold_method == "relative": default threshold_value == 0.9
    # (0.9: 90% of the maximum PWM score)

    pwm <- t(motifMatrix(motif))
    rownames(pwm)[4] <- "T"
    if (threshold_method == "p_value") {
        threshold_score <- TFMPvalue::TFMpv2sc(
            pwm, threshold_value,
            c(A = 0.25, C = 0.25, G = 0.25, T = 0.25), "PWM"
        )
    } else if (threshold_method == "relative") {
        threshold_score <- Biostrings::maxScore(pwm) * threshold_value
    } else {
        stop("invalid value threshold method")
    }

    if (threshold_score < 0) {
        warning(paste0(
            motifId(motif), ": threshold score below zero (",
            threshold_score, "), set to 0.1"
        ))
        threshold_score <- 0.1
    }

    total_sites <- unlist(lapply(sequences, function(sequence) {
        if (nchar(sequence) < motifLength(motif)) {
            return(0)
        } else {
            return(nchar(sequence) - motifLength(motif) + 1)
        }
    }))
    sum_total_sites <- sum(total_sites)

    if (!is.null(cache_path)) {
        motif_id_file <- gsub("[^[:alnum:]]", "_", motifId(motif))
        motif_cache_file <- paste0(cache_path, motif_id_file, ".rds")
        if (file.exists(motif_cache_file)) {
            motif_cache <- read_motif_cache(cache_path, motif_id_file)

            cached <- vapply(names(sequences), function(seq_id) {
                return(exists(seq_id, envir = motif_cache, inherits = FALSE))
            }, logical(1))

            cached_ids <- names(sequences)[cached]
            if (length(cached_ids) > 0) {
                cached_absolute_hits <- unlist(lapply(cached_ids,
                                                      function(seq_id) {
                    return(get(seq_id, envir = motif_cache, inherits = FALSE))
                }))
            } else {
                cached_absolute_hits <- 0
            }

            uncached_ids <- names(sequences)[!cached]
            if (length(uncached_ids) > 0) {
                uncached_sequences <- sequences[!cached]
                uncached_absolute_hits <- cached_score_sequences_helper(
                    uncached_sequences, uncached_ids,
                    as.matrix(motifMatrix(motif)),
                    threshold_score, motif_cache,
                    cache_path, motif_id_file
                )
            } else {
                uncached_absolute_hits <- 0
            }

            absolute_hits <- c(cached_absolute_hits, uncached_absolute_hits)
        } else {
            motif_cache <- new.env(hash = TRUE, parent = emptyenv(),
                                   size = length(sequences))
            absolute_hits <- cached_score_sequences_helper(
                sequences, names(sequences),
                as.matrix(motifMatrix(motif)),
                threshold_score, motif_cache,
                cache_path, motif_id_file
            )
        }
    } else {
        absolute_hits <- score_sequences_helper(sequences,
                                              as.matrix(motifMatrix(motif)),
                                              threshold_score)
    }

    if (max_hits > 9) {
        max_hits_col <- 9
    } else {
        max_hits_col <- max_hits
    }
    max_hits_col <- max_hits_col + 1
    multiple_hits <- unlist(lapply(seq_len(max_hits_col), function(hit_count) {
        if (hit_count == max_hits_col) {
            sum(vapply(absolute_hits, function(x) {
                x >= hit_count
            }, logical(1)))
        } else {
            sum(vapply(absolute_hits, function(x) {
                x == hit_count
            }, logical(1)))
        }
    }))

    absolute_hits[absolute_hits > max_hits] <- max_hits
    sum_absolute_hits <- sum(absolute_hits)

    relative_hits <- sum_absolute_hits / sum_total_sites

    if (max_hits_col == 2) {
        return(list(
            df = list(
                motif_id = motifId(motif),
                motif_rbps = paste0(motifRbps(motif), collapse = ", "),
                absolute_hits = sum_absolute_hits,
                relative_hits = relative_hits,
                total_sites = sum_total_sites,
                one_hit = multiple_hits[1], more_hits = multiple_hits[2]
            ),
            total_sites = total_sites,
            absolute_hits = absolute_hits
        ))
    } else if (max_hits_col == 3) {
        return(list(
            df = list(
                motif_id = motifId(motif),
                motif_rbps = paste0(motifRbps(motif), collapse = ", "),
                absolute_hits = sum_absolute_hits,
                relative_hits = relative_hits,
                total_sites = sum_total_sites,
                one_hit = multiple_hits[1], two_hits = multiple_hits[2],
                more_hits = multiple_hits[3]
            ),
            total_sites = total_sites,
            absolute_hits = absolute_hits
        ))
    } else if (max_hits_col == 4) {
        return(list(
            df = list(
                motif_id = motifId(motif),
                motif_rbps = paste0(motifRbps(motif), collapse = ", "),
                absolute_hits = sum_absolute_hits,
                relative_hits = relative_hits,
                total_sites = sum_total_sites,
                one_hit = multiple_hits[1], two_hits = multiple_hits[2],
                three_hits = multiple_hits[3], more_hits = multiple_hits[4]
            ),
            total_sites = total_sites,
            absolute_hits = absolute_hits
        ))
    } else if (max_hits_col == 5) {
        return(list(
            df = list(
                motif_id = motifId(motif),
                motif_rbps = paste0(motifRbps(motif), collapse = ", "),
                absolute_hits = sum_absolute_hits,
                relative_hits = relative_hits,
                total_sites = sum_total_sites,
                one_hit = multiple_hits[1], two_hits = multiple_hits[2],
                three_hits = multiple_hits[3], four_hits = multiple_hits[4],
                more_hits = multiple_hits[5]
            ),
            total_sites = total_sites,
            absolute_hits = absolute_hits
        ))
    } else if (max_hits_col == 6) {
        return(list(
            df = list(
                motif_id = motifId(motif),
                motif_rbps = paste0(motifRbps(motif), collapse = ", "),
                absolute_hits = sum_absolute_hits,
                relative_hits = relative_hits,
                total_sites = sum_total_sites,
                one_hit = multiple_hits[1], two_hits = multiple_hits[2],
                three_hits = multiple_hits[3], four_hits = multiple_hits[4],
                five_hits = multiple_hits[5], more_hits = multiple_hits[6]
            ),
            total_sites = total_sites,
            absolute_hits = absolute_hits
        ))
    } else if (max_hits_col == 7) {
        return(list(
            df = list(
                motif_id = motifId(motif),
                motif_rbps = paste0(motifRbps(motif), collapse = ", "),
                absolute_hits = sum_absolute_hits,
                relative_hits = relative_hits,
                total_sites = sum_total_sites,
                one_hit = multiple_hits[1], two_hits = multiple_hits[2],
                three_hits = multiple_hits[3], four_hits = multiple_hits[4],
                five_hits = multiple_hits[5], six_hits = multiple_hits[6],
                more_hits = multiple_hits[7]
            ),
            total_sites = total_sites,
            absolute_hits = absolute_hits
        ))
    } else if (max_hits_col == 8) {
        return(list(
            df = list(
                motif_id = motifId(motif),
                motif_rbps = paste0(motifRbps(motif), collapse = ", "),
                absolute_hits = sum_absolute_hits,
                relative_hits = relative_hits,
                total_sites = sum_total_sites,
                one_hit = multiple_hits[1], two_hits = multiple_hits[2],
                three_hits = multiple_hits[3], four_hits = multiple_hits[4],
                five_hits = multiple_hits[5], six_hits = multiple_hits[6],
                seven_hits = multiple_hits[7], more_hits = multiple_hits[8]
            ),
            total_sites = total_sites,
            absolute_hits = absolute_hits
        ))
    } else if (max_hits_col == 9) {
        return(list(
            df = list(
                motif_id = motifId(motif),
                motif_rbps = paste0(motifRbps(motif), collapse = ", "),
                absolute_hits = sum_absolute_hits,
                relative_hits = relative_hits,
                total_sites = sum_total_sites,
                one_hit = multiple_hits[1], two_hits = multiple_hits[2],
                three_hits = multiple_hits[3], four_hits = multiple_hits[4],
                five_hits = multiple_hits[5], six_hits = multiple_hits[6],
                seven_hits = multiple_hits[7], eight_hit = multiple_hits[8],
                more_hits = multiple_hits[9]
            ),
            total_sites = total_sites,
            absolute_hits = absolute_hits
        ))
    } else if (max_hits_col == 10) {
        return(list(
            df = list(
                motif_id = motifId(motif),
                motif_rbps = paste0(motifRbps(motif), collapse = ", "),
                absolute_hits = sum_absolute_hits,
                relative_hits = relative_hits,
                total_sites = sum_total_sites,
                one_hit = multiple_hits[1], two_hits = multiple_hits[2],
                three_hits = multiple_hits[3], four_hits = multiple_hits[4],
                five_hits = multiple_hits[5], six_hits = multiple_hits[6],
                seven_hits = multiple_hits[7], eight_hit = multiple_hits[8],
                nine_hits = multiple_hits[9], more_hits = multiple_hits[10]
            ),
            total_sites = total_sites,
            absolute_hits = absolute_hits
        ))
    }
}

score_sequences_helper <- function(sequences, motif_matrix, threshold_score) {
    seq_char_vectors <- lapply(sequences, function(seq) {
        unlist(strsplit(seq, ""))
    })

    scores <- score_sequences(seq_char_vectors, motif_matrix)

    absolute_hits <- unlist(lapply(scores, function(scores_per_seq) {
        sum(scores_per_seq >= threshold_score)
    }))

    return(absolute_hits)
}

cached_score_sequences_helper <- function(sequences, seq_ids, motif_matrix,
                                       threshold_score,
                                       motif_cache, cache_path, motif_id_file) {
    absolute_hits <- score_sequences_helper(sequences, motif_matrix,
                                          threshold_score)

    for (i in seq_len(length(absolute_hits))) {
        assign(seq_ids[i], absolute_hits[i], envir = motif_cache)
    }

    write_motif_cache(motif_cache, cache_path, motif_id_file)
    return(absolute_hits)
}

#' @title Binding Site Enrichment Value Calculation
#'
#' @description
#' This function is used to calculate binding site enrichment / depletion scores
#' between predefined foreground and background sequence sets. Significance
#' levels of
#' enrichment values are obtained by Monte Carlo tests.
#'
#' @param foreground_scores_df result of \code{\link{score_transcripts}} on
#' foreground sequence
#' set (foreground sequence sets must be a subset of the background
#' sequence set)
#' @param background_scores_df result of \code{\link{score_transcripts}}
#' on background sequence set
#' @param background_total_sites number of potential binding sites per sequence
#' (returned by \code{\link{score_transcripts}})
#' @param background_absolute_hits number of putative binding sites per sequence
#' (returned by \code{\link{score_transcripts}})
#' @param n_transcripts_foreground number of sequences in the foreground set
#' @param max_fg_permutations maximum number of foreground permutations
#' performed in
#' Monte Carlo test for enrichment score
#' @param min_fg_permutations minimum number of foreground permutations
#' performed in
#' Monte Carlo test for enrichment score
#' @param e integer-valued stop criterion for enrichment score Monte Carlo test: aborting
#' permutation process after
#' observing \code{e} random enrichment values with more extreme values than
#' the actual
#' enrichment value
#' @param p_adjust_method adjustment of p-values from Monte Carlo tests to
#' avoid alpha error
#'  accumulation, see \code{\link[stats]{p.adjust}}
#' @return A data frame with the following columns:
#' \tabular{rl}{
#'   \code{motif_id} \tab the motif identifier that is used in the original
#'   motif library\cr
#'   \code{motif_rbps} \tab the gene symbol of the RNA-binding protein(s)\cr
#'   \code{enrichment} \tab binding site enrichment between foreground
#'   and background sequences \cr
#'   \code{p_value} \tab unadjusted p-value from Monte Carlo test \cr
#'   \code{p_value_n} \tab number of Monte Carlo test permutations \cr
#'   \code{adj_p_value} \tab adjusted p-value from Monte Carlo test
#'   (usually FDR)
#' }
#' @examples
#' foreground_seqs <- c("CAGUCAAGACUCC", "AAUUGGUGUCUGGAUACUUCCCUGUACAU",
#'   "AGAU", "CCAGUAA")
#' background_seqs <- c(foreground_seqs, "CAACAGCCUUAAUU", "CUUUGGGGAAU",
#'                      "UCAUUUUAUUAAA", "AUCAAAUUA", "GACACUUAAAGAUCCU",
#'                      "UAGCAUUAACUUAAUG", "AUGGA", "GAAGAGUGCUCA",
#'                      "AUAGAC", "AGUUC")
#' foreground_scores <- score_transcripts(foreground_seqs, cache = FALSE)
#' background_scores <- score_transcripts(background_seqs, cache = FALSE)
#' enrichments_df <- calculate_motif_enrichment(foreground_scores$df,
#'   background_scores$df,
#'   background_scores$total_sites, background_scores$absolute_hits,
#'   length(foreground_seqs),
#'   max_fg_permutations = 1000
#' )
#' @family matrix functions
#' @importFrom dplyr filter
#' @importFrom stats p.adjust
#' @export
calculate_motif_enrichment <- function(foreground_scores_df,
                                     background_scores_df,
                                     background_total_sites,
                                     background_absolute_hits,
                                     n_transcripts_foreground,
                                     max_fg_permutations = 1000000,
                                     min_fg_permutations = 1000,
                                     e = 5, p_adjust_method = "BH") {
    # avoid CRAN note
    motif_id <- NULL

    if (!setequal(foreground_scores_df$motif_id,
                  background_scores_df$motif_id)) {
        stop("foreground and background scores used different motifs")
    }

    i <- 1
    enrichments <- lapply(foreground_scores_df$motif_id, function(id) {
        f <- dplyr::filter(foreground_scores_df, motif_id == id)
        b <- dplyr::filter(background_scores_df, motif_id == id)

        if (b$absolute_hits == 0) {
            enrichment <- 1
            p_value <- 1
            p_value_n <- 0
        } else {
            enrichment <- (f$absolute_hits / f$total_sites) / (b$absolute_hits / b$total_sites)
            mc_result <- calculate_transcript_mc(
                background_absolute_hits[[i]], background_total_sites[[i]],
                f$absolute_hits / f$total_sites, n_transcripts_foreground,
                max_fg_permutations, min_fg_permutations, e
            )
            p_value <- mc_result$p_value
            p_value_n <- mc_result$n
        }

        i <<- i + 1
        return(list(
            motif_id = id, motif_rbps = f$motif_rbps,
            enrichment = enrichment, p_value = p_value, p_value_n = p_value_n
        ))
    })

    enrichment_df <- as.data.frame(do.call("rbind", enrichments),
                                   stringsAsFactors = FALSE)
    enrichment_df$motif_id <- as.character(enrichment_df$motif_id)
    enrichment_df$motif_rbps <- as.character(enrichment_df$motif_rbps)
    enrichment_df$enrichment <- as.numeric(enrichment_df$enrichment)
    enrichment_df$p_value <- as.numeric(enrichment_df$p_value)
    enrichment_df$p_value_n <- as.numeric(enrichment_df$p_value_n)
    enrichment_df$adj_p_value <- stats::p.adjust(enrichment_df$p_value,
                                                 method = p_adjust_method)

    return(enrichment_df)
}

read_motif_cache <- function(cache_path, motif_id) {
    motif_cache_file <- paste0(cache_path, motif_id, ".rds")
    succeeded <- FALSE
    while (!succeeded) {
        tryCatch({
            motif_cache_lock <- lock(paste0(cache_path, motif_id, ".lock"))
            if (motif_cache_lock$success) {
                motif_cache <- readRDS(motif_cache_file)
                if (!unlock(motif_cache_lock)) {
                    warning(paste0("cannot unlock cache file lock of ",
                                   motif_id, " motif"))
                }
                succeeded <- TRUE
            }
        }, finally = {
            unlock(motif_cache_lock)
        })
        Sys.sleep(0.05)
    }
    return(motif_cache)
}

write_motif_cache <- function(motif_cache, cache_path, motif_id) {
    motif_cache_file <- paste0(cache_path, motif_id, ".rds")
    succeeded <- FALSE
    while (!succeeded) {
        tryCatch({
            motif_cache_lock <- lock(paste0(cache_path, motif_id, ".lock"))
            if (motif_cache_lock$success) {
                saveRDS(motif_cache, file = motif_cache_file)
                if (!unlock(motif_cache_lock)) {
                    warning(paste0("cannot unlock cache file lock of ",
                                   motif_id, " motif"))
                }
                succeeded <- TRUE
            }
        }, finally = {
            unlock(motif_cache_lock)
        })
        Sys.sleep(0.05)
    }
}
