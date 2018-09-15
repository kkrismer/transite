#' @title Scores transcripts with position weight matrices
#'
#' @description
#' This function is used to count the binding sites in a set of sequences for all or a
#' subset of RNA-binding protein sequence
#' motifs and returns the result in a data frame, which is subsequently used by
#' \code{\link{calculateMotifEnrichment}} to
#' obtain binding site enrichment scores.
#'
#' @param motifs a list of motifs that is used to score the specified sequences.
#' If \code{is.null(motifs)} then all Transite motifs are used.
#' @param n.cores the number of cores that are used
#' @param cache either logical or path to a directory where scores are cached. The scores of each
#' motif are stored in a
#' separate file that contains a hash table with RefSeq identifiers and sequence type
#' qualifiers as keys and the number of putative binding sites as values.
#' If \code{cache} is \code{FALSE}, scores will not be cached.
#' @inheritParams scoreTranscriptsSingleMotif
#'
#' @return A list with three entries:
#'
#' (1) df: a data frame with the following columns:
#' \tabular{rl}{
#'   \code{motif.id} \tab the motif identifier that is used in the original motif library\cr
#'   \code{motif.rbps} \tab the gene symbol of the RNA-binding protein(s)\cr
#'   \code{absolute.hits} \tab the absolute frequency of putative binding sites per motif in all
#'   transcripts \cr
#'   \code{relative.hits} \tab  the relative, i.e., absolute divided by total, frequency of
#'   binding sites per motif in all transcripts \cr
#'   \code{total.sites} \tab the total number of potential binding sites \cr
#'   \code{one.hit}, \code{two.hits}, ... \tab number of transcripts with one, two,
#'   three, ... putative binding sites
#' }
#' (2) total.sites: a numeric vector with the total number of potential binding sites
#' per transcript
#'
#' (3) absolute.hits: a numeric vector with the absolute (not relative) number of putative
#' binding sites per transcript
#' @examples
#' foreground.set <- c("CAACAGCCTTAATT", "CAGTCAAGACTCC", "CTTTGGGGAAT",
#'                                    "TCATTTTATTAAA", "AATTGGTGTCTGGATACTTCCCTGTACAT",
#'                                    "ATCAAATTA", "AGAT", "GACACTTAAAGATCCT",
#'                                    "TAGCATTAACTTAATG", "ATGGA", "GAAGAGTGCTCA",
#'                                    "ATAGAC", "AGTTC", "CCAGTAA")
#' # names are used as keys in the hash table (cached version only)
#' # ideally sequence identifiers (e.g., RefSeq ids) and region labels (e.g., 3UTR for 3'-UTR)
#' names(foreground.set) <- c("NM_1_DUMMY|3UTR", "NM_2_DUMMY|3UTR", "NM_3_DUMMY|3UTR",
#'   "NM_4_DUMMY|3UTR", "NM_5_DUMMY|3UTR", "NM_6_DUMMY|3UTR", "NM_7_DUMMY|3UTR",
#'   "NM_8_DUMMY|3UTR", "NM_9_DUMMY|3UTR", "NM_10_DUMMY|3UTR", "NM_11_DUMMY|3UTR",
#'   "NM_12_DUMMY|3UTR", "NM_13_DUMMY|3UTR", "NM_14_DUMMY|3UTR")
#'
#' # cached (writes scores to disk)
#' scores <- scoreTranscripts(foreground.set)
#'
#' # uncached
#' scores <- scoreTranscripts(foreground.set, cache = FALSE)
#'
#' # specific motifs
#' motifs <- transite::getMotifByRBP("ELAVL1")
#' scores <- scoreTranscripts(foreground.set, motifs = motifs)
#'
#' \dontrun{
#' foreground.df <- ge$foreground1
#' foreground.set <- foreground.df$seq
#' names(foreground.set) <- paste0(foreground.df$refseq, "|", foreground.df$seq.type)
#' scores <- scoreTranscripts(foreground.set)
#' }
#' @family matrix functions
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel parLapply
#' @export
scoreTranscripts <- function(sequences, motifs = NULL, max.hits = 5,
                             threshold.method = "p.value", threshold.value = 0.25^6,
                             n.cores = 1, cache = paste0(getwd(), "/sc/")) {
  # if threshold.method == "p.value": default threshold.value == 0.25^6
  # (0.25^6: lowest p-value that can be achieved by hexamer motifs, the shortest supported motifs)
  # if threshold.method == "relative": default threshold.value == 0.9
  # (0.9: 90% of the maximum PWM score)
  if(max.hits < 1) {
    stop("max.hits must not be smaller than 1")
  }

  if(is.null(motifs)) {
    motifs <- getMotifs()
  }

  if(is.logical(cache)) {
    if(cache) {
      cache.path <- paste0(getwd(), "/sc/")
    } else {
      cache.path <- NULL
    }
  } else {
    cache.path <- as.character(cache)
  }

  if(!is.null(cache.path)) {
    if(is.null(names(sequences))) {
      stop("names(sequences) cannot be null if cache is used")
    }
    names(sequences) <- paste0(names(sequences), "|", threshold.method, "|", threshold.value)
    dir.create(file.path(cache.path), showWarnings = FALSE, recursive = TRUE)
  }

  if(n.cores == 1) {
    motif.scores <- lapply(motifs, function(motif) {
      return(scoreTranscriptsSingleMotif(motif, sequences, max.hits, threshold.method,
                                         threshold.value, cache.path))
    })
  } else {
    cluster <- parallel::makeCluster(mc <- getOption("cl.cores", n.cores))
    parallel::clusterExport(cl = cluster,
                            varlist = c("scoreTranscriptsSingleMotif", "sequences", "max.hits",
                                        "threshold.method", "threshold.value", "cache.path"),
                            envir = environment())
    motif.scores <- parallel::parLapply(cl = cluster, motifs, function(motif) {
      return(scoreTranscriptsSingleMotif(motif, sequences, max.hits, threshold.method,
                                         threshold.value, cache.path))
    })
  }

  motif.scores.proto.df <- lapply(motif.scores, function(motif.score) {
    return(motif.score$df)
  })

  total.sites <- lapply(motif.scores, function(motif.score) {
    return(motif.score$total.sites)
  })

  absolute.hits <- lapply(motif.scores, function(motif.score) {
    return(motif.score$absolute.hits)
  })

  motif.scores.df <- as.data.frame(do.call("rbind", motif.scores.proto.df),
                                   stringsAsFactors = FALSE)
  motif.scores.df$motif.id <- as.character(motif.scores.df$motif.id)
  motif.scores.df$motif.rbps <- as.character(motif.scores.df$motif.rbps)
  motif.scores.df$absolute.hits <- as.numeric(motif.scores.df$absolute.hits)
  motif.scores.df$relative.hits <- as.numeric(motif.scores.df$relative.hits)
  motif.scores.df$total.sites <- as.numeric(motif.scores.df$total.sites)
  motif.scores.df$one.hit <- as.numeric(motif.scores.df$one.hit)
  if(max.hits > 1) {
    motif.scores.df$two.hits <- as.numeric(motif.scores.df$two.hits)
  }
  if(max.hits > 2) {
    motif.scores.df$three.hits <- as.numeric(motif.scores.df$three.hits)
  }
  if(max.hits > 3) {
    motif.scores.df$four.hits <- as.numeric(motif.scores.df$four.hits)
  }
  if(max.hits > 4) {
    motif.scores.df$five.hits <- as.numeric(motif.scores.df$five.hits)
  }
  if(max.hits > 5) {
    motif.scores.df$six.hits <- as.numeric(motif.scores.df$six.hits)
  }
  if(max.hits > 6) {
    motif.scores.df$seven.hits <- as.numeric(motif.scores.df$seven.hits)
  }
  if(max.hits > 7) {
    motif.scores.df$eight.hit <- as.numeric(motif.scores.df$eight.hit)
  }
  if(max.hits > 8) {
    motif.scores.df$nine.hits <- as.numeric(motif.scores.df$nine.hits)
  }
  motif.scores.df$more.hits <- as.numeric(motif.scores.df$more.hits)
  return(list(df = motif.scores.df, total.sites = total.sites, absolute.hits = absolute.hits))
}

#' @title Scores transadsadscripts with position weight matrices
#'
#' @description
#' This function is used to count the putative binding sites (i.e., motifs) in a set of
#' sequences for the specified RNA-binding protein sequence
#' motifs and returns the result in a data frame, which is aggregated by
#' \code{\link{scoreTranscripts}} and
#' subsequently used by \code{\link{calculateMotifEnrichment}} to
#' obtain binding site enrichment scores.
#'
#' @param motif a Transite motif that is used to score the specified sequences
#' @param sequences character vector of named sequences
#' (only containing upper case characters A, C, G, T), where the names are RefSeq identifiers
#' and sequence
#' type qualifiers (\code{"3UTR"}, \code{"5UTR"}, \code{"mRNA"}), e.g. \code{"NM_010356|3UTR"}
#' @param max.hits maximum number of putative binding sites per mRNA that are counted
#' @param cache.path the path to a directory where scores are cached. The scores of each
#' motif are stored in a
#' separate file that contains a hash table with RefSeq identifiers and sequence type
#' qualifiers as keys and the number of binding sites as values.
#' If is.null(cache.path), scores will not be cached.
#' @param threshold.method either \code{"p.value"} (default) or \code{"relative"}.
#' If \code{threshold.method} equals \code{"p.value"}, the default \code{threshold.value}
#'  is \code{0.25^6}, which is
#' lowest p-value that can be achieved by hexamer motifs, the shortest supported motifs.
#' If \code{threshold.method} equals \code{"relative"}, the default \code{threshold.value}
#' is \code{0.9}, which is 90\% of the maximum PWM score.
#' @param threshold.value semantics of the \code{threshold.value} depend on
#' \code{threshold.method} (default is 0.25^6)
#' @return A list with the following items:
#' \tabular{rl}{
#'   \code{motif.id} \tab the motif identifier of the specified motif\cr
#'   \code{motif.rbps} \tab the gene symbol of the RNA-binding protein(s)\cr
#'   \code{absolute.hits} \tab the absolute frequency of binding sites per motif in all
#'   transcripts \cr
#'   \code{relative.hits} \tab  the relative, i.e., absolute divided by total, frequency of
#'   binding sites per motif in all transcripts \cr
#'   \code{total.sites} \tab the total number of potential binding sites \cr
#'   \code{one.hit}, \code{two.hits}, ... \tab number of transcripts with one, two, three,
#'   ... binding sites
#' }
#'
#' @family matrix functions
#' @importFrom TFMPvalue TFMpv2sc
#' @importFrom Biostrings maxScore
scoreTranscriptsSingleMotif <- function(motif, sequences, max.hits = 5,
                                        threshold.method = "p.value", threshold.value = 0.25^6,
                                        cache.path = paste0(getwd(), "/sc/")) {
  # if threshold.method == "p.value": default threshold.value == 0.25^6
  # (0.25^6: lowest p-value that can be achieved by hexamer motifs, the shortest supported motifs)
  # if threshold.method == "relative": default threshold.value == 0.9
  # (0.9: 90% of the maximum PWM score)

  pwm <- t(motif$matrix)
  rownames(pwm)[4] <- "T"
  if(threshold.method == "p.value") {
    threshold.score <- TFMPvalue::TFMpv2sc(pwm, threshold.value,
                                           c(A = 0.25, C = 0.25, G = 0.25, T = 0.25), "PWM")
  } else if(threshold.method == "relative") {
    threshold.score <- Biostrings::maxScore(pwm) * threshold.value
  } else {
    stop("invalid value threshold method")
  }

  if(threshold.score < 0) {
    warning(paste0(motif$motif.id, ": threshold score below zero (",
                   threshold.score, "), set to 0.1"))
    threshold.score <- 0.1
  }

  total.sites <- unlist(lapply(sequences, function(sequence) {
    if(nchar(sequence) < motif$length) {
      return(0)
    } else {
      return(nchar(sequence) - motif$length + 1)
    }
  }))
  sum.total.sites <- sum(total.sites)

  if(!is.null(cache.path)) {
    motif.id.file <- gsub("[^[:alnum:]]", "_", motif$id)
    motif.cache.file <- paste0(cache.path, motif.id.file, ".rds")
    if(file.exists(motif.cache.file)) {
      motif.cache <- readMotifCache(cache.path, motif.id.file)

      cached <- vapply(names(sequences), function(seq.id) {
        return(exists(seq.id, envir = motif.cache, inherits = FALSE))
      })

      cached.ids <- names(sequences)[cached]
      if(length(cached.ids) > 0) {
        cached.absolute.hits <- unlist(lapply(cached.ids, function(seq.id) {
          return(get(seq.id, envir = motif.cache, inherits = FALSE))
        }))
      } else {
        cached.absolute.hits <- 0
      }

      uncached.ids <- names(sequences)[!cached]
      if(length(uncached.ids) > 0) {
        uncached.sequences <- sequences[!cached]
        uncached.absolute.hits <- cachedScoreSequencesHelper(uncached.sequences, uncached.ids,
                                                             as.matrix(motif$matrix),
                                                             threshold.score, motif.cache,
                                                             cache.path, motif.id.file)
      } else {
        uncached.absolute.hits <- 0
      }

      absolute.hits <- c(cached.absolute.hits, uncached.absolute.hits)
    } else {
      motif.cache = new.env(hash = TRUE, parent = emptyenv(), size = length(sequences))
      absolute.hits <- cachedScoreSequencesHelper(sequences, names(sequences),
                                                  as.matrix(motif$matrix),
                                                  threshold.score, motif.cache,
                                                  cache.path, motif.id.file)
    }
  } else {
    absolute.hits <- scoreSequencesHelper(sequences, as.matrix(motif$matrix), threshold.score)
  }

  if(max.hits > 9) {
    MAX.HITS.COL <- 9
  } else {
    MAX.HITS.COL <- max.hits
  }
  MAX.HITS.COL <- MAX.HITS.COL + 1
  multiple.hits <- unlist(lapply(seq_len(MAX.HITS.COL), function(hit.count) {
    if(hit.count == MAX.HITS.COL) {
      sum(vapply(absolute.hits, function(x) {
        x >= hit.count
      }))
    } else {
      sum(vapply(absolute.hits, function(x) {
        x == hit.count
      }))
    }
  }))

  absolute.hits[absolute.hits > max.hits] <- max.hits
  sum.absolute.hits <- sum(absolute.hits)

  relative.hits <- sum.absolute.hits / sum.total.sites

  if(MAX.HITS.COL == 2) {
    return(list(df = list(motif.id = motif$id,
                          motif.rbps = paste0(motif$rbps, collapse = ", "),
                          absolute.hits = sum.absolute.hits,
                          relative.hits = relative.hits,
                          total.sites = sum.total.sites,
                          one.hit = multiple.hits[1], more.hits = multiple.hits[2]),
                total.sites = total.sites,
                absolute.hits = absolute.hits))
  } else if(MAX.HITS.COL == 3) {
    return(list(df = list(motif.id = motif$id,
                          motif.rbps = paste0(motif$rbps, collapse = ", "),
                          absolute.hits = sum.absolute.hits,
                          relative.hits = relative.hits,
                          total.sites = sum.total.sites,
                          one.hit = multiple.hits[1], two.hits = multiple.hits[2],
                          more.hits = multiple.hits[3]),
                total.sites = total.sites,
                absolute.hits = absolute.hits))
  } else if(MAX.HITS.COL == 4) {
    return(list(df = list(motif.id = motif$id,
                          motif.rbps = paste0(motif$rbps, collapse = ", "),
                          absolute.hits = sum.absolute.hits,
                          relative.hits = relative.hits,
                          total.sites = sum.total.sites,
                          one.hit = multiple.hits[1], two.hits = multiple.hits[2],
                          three.hits = multiple.hits[3], more.hits = multiple.hits[4]),
                total.sites = total.sites,
                absolute.hits = absolute.hits))
  } else if(MAX.HITS.COL == 5) {
    return(list(df = list(motif.id = motif$id,
                          motif.rbps = paste0(motif$rbps, collapse = ", "),
                          absolute.hits = sum.absolute.hits,
                          relative.hits = relative.hits,
                          total.sites = sum.total.sites,
                          one.hit = multiple.hits[1], two.hits = multiple.hits[2],
                          three.hits = multiple.hits[3], four.hits = multiple.hits[4],
                          more.hits = multiple.hits[5]),
                total.sites = total.sites,
                absolute.hits = absolute.hits))
  } else if(MAX.HITS.COL == 6) {
    return(list(df = list(motif.id = motif$id,
                          motif.rbps = paste0(motif$rbps, collapse = ", "),
                          absolute.hits = sum.absolute.hits,
                          relative.hits = relative.hits,
                          total.sites = sum.total.sites,
                          one.hit = multiple.hits[1], two.hits = multiple.hits[2],
                          three.hits = multiple.hits[3], four.hits = multiple.hits[4],
                          five.hits = multiple.hits[5], more.hits = multiple.hits[6]),
                total.sites = total.sites,
                absolute.hits = absolute.hits))
  } else if(MAX.HITS.COL == 7) {
    return(list(df = list(motif.id = motif$id,
                          motif.rbps = paste0(motif$rbps, collapse = ", "),
                          absolute.hits = sum.absolute.hits,
                          relative.hits = relative.hits,
                          total.sites = sum.total.sites,
                          one.hit = multiple.hits[1], two.hits = multiple.hits[2],
                          three.hits = multiple.hits[3], four.hits = multiple.hits[4],
                          five.hits = multiple.hits[5], six.hits = multiple.hits[6],
                          more.hits = multiple.hits[7]),
                total.sites = total.sites,
                absolute.hits = absolute.hits))
  } else if(MAX.HITS.COL == 8) {
    return(list(df = list(motif.id = motif$id,
                          motif.rbps = paste0(motif$rbps, collapse = ", "),
                          absolute.hits = sum.absolute.hits,
                          relative.hits = relative.hits,
                          total.sites = sum.total.sites,
                          one.hit = multiple.hits[1], two.hits = multiple.hits[2],
                          three.hits = multiple.hits[3], four.hits = multiple.hits[4],
                          five.hits = multiple.hits[5], six.hits = multiple.hits[6],
                          seven.hits = multiple.hits[7], more.hits = multiple.hits[8]),
                total.sites = total.sites,
                absolute.hits = absolute.hits))
  } else if(MAX.HITS.COL == 9) {
    return(list(df = list(motif.id = motif$id,
                          motif.rbps = paste0(motif$rbps, collapse = ", "),
                          absolute.hits = sum.absolute.hits,
                          relative.hits = relative.hits,
                          total.sites = sum.total.sites,
                          one.hit = multiple.hits[1], two.hits = multiple.hits[2],
                          three.hits = multiple.hits[3], four.hits = multiple.hits[4],
                          five.hits = multiple.hits[5], six.hits = multiple.hits[6],
                          seven.hits = multiple.hits[7], eight.hit = multiple.hits[8],
                          more.hits = multiple.hits[9]),
                total.sites = total.sites,
                absolute.hits = absolute.hits))
  } else if(MAX.HITS.COL == 10) {
    return(list(df = list(motif.id = motif$id,
                          motif.rbps = paste0(motif$rbps, collapse = ", "),
                          absolute.hits = sum.absolute.hits,
                          relative.hits = relative.hits,
                          total.sites = sum.total.sites,
                          one.hit = multiple.hits[1], two.hits = multiple.hits[2],
                          three.hits = multiple.hits[3], four.hits = multiple.hits[4],
                          five.hits = multiple.hits[5], six.hits = multiple.hits[6],
                          seven.hits = multiple.hits[7], eight.hit = multiple.hits[8],
                          nine.hits = multiple.hits[9], more.hits = multiple.hits[10]),
                total.sites = total.sites,
                absolute.hits = absolute.hits))
  }
}

scoreSequencesHelper <- function(sequences, motif.matrix, threshold.score) {
  seq.char.vectors <- lapply(sequences, function(seq) {
    unlist(strsplit(seq, ""))
  })

  scores <- scoreSequences(seq.char.vectors, motif.matrix)

  absolute.hits <- unlist(lapply(scores, function(scores.per.seq) {
    sum(scores.per.seq >= threshold.score)
  }))

  return(absolute.hits)
}

cachedScoreSequencesHelper <- function(sequences, seq.ids, motif.matrix, threshold.score,
                                       motif.cache, cache.path, motif.id.file) {
  absolute.hits <- scoreSequencesHelper(sequences, motif.matrix, threshold.score)

  for(i in seq_len(length(absolute.hits))) {
    assign(seq.ids[i], absolute.hits[i], envir = motif.cache)
  }

  writeMotifCache(motif.cache, cache.path, motif.id.file)
  return(absolute.hits)
}

#' @title Binding Site Enrichment Value Calculation
#'
#' @description
#' This function is used to calculate binding site enrichment / depletion scores
#' between predefined foreground and background sequence sets. Significance levels of
#' enrichment values are obtained by Monte Carlo tests.
#'
#' @param foreground.scores.df result of \code{\link{scoreTranscripts}} on foreground sequence
#' set (foreground sequence sets must be a subset of the background sequence set)
#' @param background.scores.df result of \code{\link{scoreTranscripts}} on background sequence set
#' @param background.total.sites number of potential binding sites per sequence
#' (returned by \code{\link{scoreTranscripts}})
#' @param background.absolute.hits number of putative binding sites per sequence
#' (returned by \code{\link{scoreTranscripts}})
#' @param n.transcripts.foreground number of sequences in the foreground set
#' @param max.fg.permutations maximum number of foreground permutations performed in
#' Monte Carlo test for enrichment score
#' @param min.fg.permutations minimum number of foreground permutations performed in
#' Monte Carlo test for enrichment score
#' @param e stop criterion for enrichment score Monte Carlo test: aborting permutation process after
#' observing \code{e} random enrichment values with more extreme values than the actual
#' enrichment value
#' @param p.adjust.method adjustment of p-values from Monte Carlo tests to avoid alpha error
#'  accumulation, see \code{\link[stats]{p.adjust}}
#' @return A data frame with the following columns:
#' \tabular{rl}{
#'   \code{motif.id} \tab the motif identifier that is used in the original motif library\cr
#'   \code{motif.rbps} \tab the gene symbol of the RNA-binding protein(s)\cr
#'   \code{enrichment} \tab binding site enrichment between foreground and background sequences \cr
#'   \code{p.value} \tab unadjusted p-value from Monte Carlo test \cr
#'   \code{p.value.n} \tab number of Monte Carlo test permutations \cr
#'   \code{adj.p.value} \tab adjusted p-value from Monte Carlo test (usually FDR)
#' }
#' @examples
#' foreground.seqs <- c("CAGTCAAGACTCC", "AATTGGTGTCTGGATACTTCCCTGTACAT", "AGAT", "CCAGTAA")
#' background.seqs <- c("CAACAGCCTTAATT", "CAGTCAAGACTCC", "CTTTGGGGAAT",
#'                                               "TCATTTTATTAAA", "AATTGGTGTCTGGATACTTCCCTGTACAT",
#'                                               "ATCAAATTA", "AGAT", "GACACTTAAAGATCCT",
#'                                               "TAGCATTAACTTAATG", "ATGGA", "GAAGAGTGCTCA",
#'                                               "ATAGAC", "AGTTC", "CCAGTAA")
#' foreground.scores <- scoreTranscripts(foreground.seqs, cache = FALSE)
#' background.scores <- scoreTranscripts(background.seqs, cache = FALSE)
#' enrichments.df <- calculateMotifEnrichment(foreground.scores$df, background.scores$df,
#'   background.scores$total.sites, background.scores$absolute.hits, length(foreground.seqs),
#'   max.fg.permutations = 1000)
#' @family matrix functions
#' @importFrom dplyr filter
#' @importFrom stats p.adjust
#' @export
calculateMotifEnrichment <- function(foreground.scores.df, background.scores.df,
                                background.total.sites, background.absolute.hits,
                                n.transcripts.foreground,
                                max.fg.permutations = 1000000, min.fg.permutations = 1000,
                                e = 5, p.adjust.method = "BH") {
  # avoid CRAN note
  motif.id <- NULL

  if(!setequal(foreground.scores.df$motif.id, background.scores.df$motif.id)) {
    stop("foreground and background scores used different motifs")
  }

  i <- 1
  enrichments <- lapply(foreground.scores.df$motif.id, function(id) {
    f <- dplyr::filter(foreground.scores.df, motif.id == id)
    b <- dplyr::filter(background.scores.df, motif.id == id)

    if(b$absolute.hits == 0) {
      enrichment <- 1
      p.value <- 1
      p.value.n <- 0
    } else {
      enrichment <- (f$absolute.hits / f$total.sites) / (b$absolute.hits / b$total.sites)
      mc.result <- calculateTranscriptMC(background.absolute.hits[[i]], background.total.sites[[i]],
                                         f$absolute.hits / f$total.sites, n.transcripts.foreground,
                                         max.fg.permutations, min.fg.permutations, e)
      p.value <- mc.result$p.value
      p.value.n <- mc.result$n
#       cont <- matrix(c(f$absolute.hits, f$total.sites - f$absolute.hits,
#                        b$absolute.hits, b$total.sites - b$absolute.hits), nrow = 2)
#       p.value <- fisher.test(cont)$p.value
    }

    i <<- i + 1
    return(list(motif.id = id, motif.rbps = f$motif.rbps,
                enrichment = enrichment, p.value = p.value, p.value.n = p.value.n))
  })

  enrichment.df <- as.data.frame(do.call("rbind", enrichments), stringsAsFactors = FALSE)
  enrichment.df$motif.id <- as.character(enrichment.df$motif.id)
  enrichment.df$motif.rbps <- as.character(enrichment.df$motif.rbps)
  enrichment.df$enrichment <- as.numeric(enrichment.df$enrichment)
  enrichment.df$p.value <- as.numeric(enrichment.df$p.value)
  enrichment.df$p.value.n <- as.numeric(enrichment.df$p.value.n)
  enrichment.df$adj.p.value <- stats::p.adjust(enrichment.df$p.value, method = p.adjust.method)

  return(enrichment.df)
}

readMotifCache <- function(cache.path, motif.id) {
  motif.cache.file <- paste0(cache.path, motif.id, ".rds")
  succeeded <- FALSE
  while(!succeeded) {
    tryCatch({
      motif.cache.lock <- lock(paste0(cache.path, motif.id, ".lock"))
      if(motif.cache.lock$success) {
        motif.cache <- readRDS(motif.cache.file)
        if(!unlock(motif.cache.lock)) {
          warning(paste0("cannot unlock cache file lock of ", motif.id, " motif"))
        }
        succeeded <- TRUE
      }
    }, finally = {unlock(motif.cache.lock)})
    Sys.sleep(0.1)
  }
  return(motif.cache)
}

writeMotifCache <- function(motif.cache, cache.path, motif.id) {
  motif.cache.file <- paste0(cache.path, motif.id, ".rds")
  succeeded <- FALSE
  while(!succeeded) {
    tryCatch({
      motif.cache.lock <- lock(paste0(cache.path, motif.id, ".lock"))
      if(motif.cache.lock$success) {
        saveRDS(motif.cache, file = motif.cache.file)
        if(!unlock(motif.cache.lock)) {
          warning(paste0("cannot unlock cache file lock of ", motif.id, " motif"))
        }
        succeeded <- TRUE
      }
    }, finally = {unlock(motif.cache.lock)})
    Sys.sleep(0.1)
  }
}
