#' @title Correction for Homopolymeric Stretches
#'
#' @description
#' Counts all non-overlapping instances of \emph{k}-mers in a given set
#' of sequences.
#'
#' @param sequences character vector of DNA or RNA sequences
#' @param k length of \emph{k}-mer, either \code{6} for hexamers or
#' \code{7} for heptamers
#' @param kmers column sums of return value of
#' \code{Biostrings::oligonucleotideFrequency(sequences)}
#' @param is_rna if \code{sequences} are RNA sequences, this
#' flag needs to be set
#'
#' @return Returns a named numeric vector, where the elements are
#' \emph{k}-mer counts and the names are \emph{k}-mers.
#'
#' @importFrom methods as
#' @importFrom Biostrings maskMotif
#' @importFrom BiocGenerics width
#' @importFrom GenomicRanges gaps
#' @importFrom Biostrings BString
#' @family \emph{k}-mer functions
count_homopolymer_corrected_kmers <-
    function(sequences, k, kmers, is_rna = FALSE) {
        seq <- Biostrings::BString(toString(sequences))

        aMask <- paste(rep("A", k), collapse = "")
        maska <- Biostrings::maskMotif(seq, aMask)
        mwidthsa <-
            BiocGenerics::width(methods::as(
                GenomicRanges::gaps(maska),
                "Views"
            ))
        viewcnta <-
            length(methods::as(GenomicRanges::gaps(maska), "Views"))
        kmers[aMask] <-
            kmers[aMask] - (sum(mwidthsa) - viewcnta * (k - 1)) + sum(floor(mwidthsa / k))

        cMask <- paste(rep("C", k), collapse = "")
        maskc <- Biostrings::maskMotif(seq, cMask)
        mwidthsc <-
            BiocGenerics::width(methods::as(GenomicRanges::gaps(maskc), "Views"))
        viewcntc <-
            length(methods::as(GenomicRanges::gaps(maskc), "Views"))
        kmers[cMask] <-
            kmers[cMask] - (sum(mwidthsc) - viewcntc * (k - 1)) + sum(floor(mwidthsc / k))

        gMask <- paste(rep("G", k), collapse = "")
        maskg <- Biostrings::maskMotif(seq, gMask)
        mwidthsg <-
            BiocGenerics::width(methods::as(GenomicRanges::gaps(maskg), "Views"))
        viewcntg <-
            length(methods::as(GenomicRanges::gaps(maskg), "Views"))
        kmers[gMask] <-
            kmers[gMask] - (sum(mwidthsg) - viewcntg * (k - 1)) + sum(floor(mwidthsg / k))

        if (is_rna) {
            uMask <- paste(rep("U", k), collapse = "")
            masku <- Biostrings::maskMotif(seq, uMask)
            mwidthsu <-
                BiocGenerics::width(methods::as(GenomicRanges::gaps(masku), "Views"))
            viewcntu <-
                length(methods::as(GenomicRanges::gaps(masku), "Views"))
            kmers[uMask] <-
                kmers[uMask] - (sum(mwidthsu) - viewcntu * (k - 1)) + sum(floor(mwidthsu / k))
        } else {
            tMask <- paste(rep("T", k), collapse = "")
            maskt <- Biostrings::maskMotif(seq, tMask)
            mwidthst <-
                BiocGenerics::width(methods::as(GenomicRanges::gaps(maskt), "Views"))
            viewcntt <-
                length(methods::as(GenomicRanges::gaps(maskt), "Views"))
            kmers[tMask] <-
                kmers[tMask] - (sum(mwidthst) - viewcntt * (k - 1)) + sum(floor(mwidthst / k))
        }

        return(kmers)
    }

#' @title \emph{k}-mer Enrichment between Foreground and Background Sets
#'
#' @description
#' Calls \code{\link{compute_kmer_enrichment}} to compute \emph{k}-mer
#' enrichment values
#' for multiple foregrounds. Calculates enrichment for foreground sets in
#' parallel.
#'
#' @param foreground_sets list of foreground sets; a foreground set is a
#' character vector of
#' DNA or RNA sequences (not both) and a strict subset of the
#' \code{background_set}
#' @param background_set character vector of DNA or RNA sequences that
#' constitute the
#' background set
#' @param n_cores number of computing cores to use
#' @inheritParams generate_kmers
#' @inheritParams compute_kmer_enrichment
#'
#' @return A list with two entries:
#'
#' (1) dfs: a list of data frames with results from
#' \code{\link{compute_kmer_enrichment}} for each of the foreground sets
#' (2) kmers: a character vector of all k-mers
#'
#' @examples
#' # define simple sequence sets for foreground and background
#' foreground_set1 <- c(
#'   "CAACAGCCUUAAUU", "CAGUCAAGACUCC", "CUUUGGGGAAU",
#'   "UCAUUUUAUUAAA", "AAUUGGUGUCUGGAUACUUCCCUGUACAU",
#'   "AUCAAAUUA", "AGAU", "GACACUUAAAGAUCCU",
#'   "UAGCAUUAACUUAAUG", "AUGGA", "GAAGAGUGCUCA",
#'   "AUAGAC", "AGUUC", "CCAGUAA"
#' )
#' foreground_set2 <- c("UUAUUUA", "AUCCUUUACA", "UUUUUUU", "UUUCAUCAUU")
#' foreground_sets <- list(foreground_set1, foreground_set2)
#' background_set <- c(foreground_set1, foreground_set2,
#'                     "CCACACAC", "CUCAUUGGAG", "ACUUUGGGACA", "CAGGUCAGCA")
#'
#' # single-threaded
#' kmer_enrichment_values_st <- calculate_kmer_enrichment(foreground_sets,
#'   background_set, 6, n_cores = 1)
#' \dontrun{
#' # multi-threaded
#' kmer_enrichment_values_mt <- calculate_kmer_enrichment(foreground_sets,
#'   background_set, 6)}
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel parLapply
#' @family \emph{k}-mer functions
#' @export
calculate_kmer_enrichment <-
    function(foreground_sets,
             background_set,
             k,
             permutation = FALSE,
             chisq_p_value_threshold = 0.05,
             p_adjust_method = "BH",
             n_cores = 4) {
        if (length(foreground_sets) > 0) {
            background_kmers <- generate_kmers(background_set, k)

            n_cores <- min(n_cores, length(foreground_sets))
            if (n_cores == 1) {
                enrichment_dfs <- lapply(foreground_sets,
                                         function(foreground_set) {
                    foreground_kmers <- generate_kmers(foreground_set, k)
                    return(
                        compute_kmer_enrichment(
                            foreground_kmers,
                            background_kmers,
                            permutation = permutation,
                            chisq_p_value_threshold = chisq_p_value_threshold,
                            p_adjust_method = p_adjust_method
                        )
                    )
                })
            } else {
                cluster <-
                    parallel::makeCluster(mc <-
                                              getOption("cl.cores", n_cores))
                parallel::clusterExport(
                    cl = cluster,
                    varlist = c(
                        "k",
                        "background_kmers",
                        "permutation",
                        "chisq_p_value_threshold",
                        "p_adjust_method",
                        "compute_kmer_enrichment",
                        "generate_kmers",
                        "count_homopolymer_corrected_kmers"
                    ),
                    envir = environment()
                )
                enrichment_dfs <-
                    parallel::parLapply(
                        cl = cluster, foreground_sets,
                        function(foreground_set) {
                            foreground_kmers <- generate_kmers(foreground_set, k)
                            return(
                                compute_kmer_enrichment(
                                    foreground_kmers,
                                    background_kmers,
                                    permutation = permutation,
                                    chisq_p_value_threshold = chisq_p_value_threshold,
                                    p_adjust_method = p_adjust_method
                                )
                            )
                        }
                    )
            }
            return(list(
                dfs = enrichment_dfs,
                kmers = gsub("T", "U", as.character(names(
                    background_kmers
                )))
            ))
        } else {
            return(NULL)
        }
    }

#' @title \emph{k}-mer Counts for Sequence Set
#'
#' @description
#' Counts occurrences of \emph{k}-mers of length \code{k} in the given set of
#' sequences. Corrects for homopolymeric stretches.
#'
#' @inheritParams count_homopolymer_corrected_kmers
#'
#' @return Returns a named numeric vector, where the elements are
#' \emph{k}-mer counts and the
#' names are DNA \emph{k}-mers.
#'
#' @examples
#' # count hexamers in set of RNA sequences
#' rna_sequences <- c(
#'   "CAACAGCCUUAAUU", "CAGUCAAGACUCC", "CUUUGGGGAAU",
#'   "UCAUUUUAUUAAA", "AAUUGGUGUCUGGAUACUUCCCUGUACAU",
#'   "AUCAAAUUA", "AGAU", "GACACUUAAAGAUCCU",
#'   "UAGCAUUAACUUAAUG", "AUGGA", "GAAGAGUGCUCA",
#'   "AUAGAC", "AGUUC", "CCAGUAA",
#'   "UUAUUUA", "AUCCUUUACA", "UUUUUUU", "UUUCAUCAUU",
#'   "CCACACAC", "CUCAUUGGAG", "ACUUUGGGACA", "CAGGUCAGCA"
#' )
#' hexamer_counts <- generate_kmers(rna_sequences, 6)
#'
#'
#' # count heptamers in set of DNA sequences
#' dna_sequences <- c(
#'   "CAACAGCCTTAATT", "CAGTCAAGACTCC", "CTTTGGGGAAT",
#'   "TCATTTTATTAAA", "AATTGGTGTCTGGATACTTCCCTGTACAT",
#'   "ATCAAATTA", "AGAT", "GACACTTAAAGATCCT",
#'   "TAGCATTAACTTAATG", "ATGGA", "GAAGAGTGCTCA",
#'   "ATAGAC", "AGTTC", "CCAGTAA",
#'   "TTATTTA", "ATCCTTTACA", "TTTTTTT", "TTTCATCATT",
#'   "CCACACAC", "CTCATTGGAG", "ACTTTGGGACA", "CAGGTCAGCA"
#' )
#' hexamer_counts <- generate_kmers(dna_sequences, 7)
#' @section Warning:
#' \code{generate_kmers} always returns DNA \emph{k}-mers, even if
#' \code{sequences} contains RNA sequences.
#' RNA sequences are internally converted to DNA sequences. It is not
#' allowed to mix DNA and
#' RNA sequences.
#'
#' @importFrom Biostrings oligonucleotideFrequency
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings RNAStringSet
#' @family \emph{k}-mer functions
#' @export
generate_kmers <- function(sequences, k) {
    kmers <- tryCatch({
        kmer_counts <-
            Biostrings::oligonucleotideFrequency(Biostrings::DNAStringSet(sequences), k)
        count_homopolymer_corrected_kmers(sequences, k, colSums(kmer_counts))
    }, error = function(e) {
        kmer_counts <-
            Biostrings::oligonucleotideFrequency(Biostrings::RNAStringSet(sequences), k)
        rna_kmers <-
            count_homopolymer_corrected_kmers(sequences, k, colSums(kmer_counts), is_rna = TRUE)
        names(rna_kmers) <- gsub("U", "T", names(rna_kmers))
        return(rna_kmers)
    })

    return(kmers)
}

#' @title \emph{k}-mer Enrichment between Foreground and Background Sets
#'
#' @description
#' Compares foreground sequence set to background sequence set and computes
#' enrichment values
#' for each possible \emph{k}-mer.
#'
#' @param foreground_kmers \emph{k}-mer counts of the foreground set
#' (generated by \code{\link{generate_kmers}})
#' @param background_kmers \emph{k}-mer counts of the background set
#' (generated by \code{\link{generate_kmers}})
#' @param permutation if \code{TRUE}, only the enrichment value is returned
#' (efficiency mode
#' used for permutation testing)
#' @param chisq_p_value_threshold threshold below which Fisher's exact test
#' is used instead of Pearson's chi-squared test
#' @param p_adjust_method see \code{\link[stats]{p.adjust}}
#'
#' @return enrichment of \emph{k}-mers in specified foreground sequences.
#' A data frame with the following columns is returned:
#' \tabular{rl}{
#'   \code{foreground_count} \tab foreground counts for each \emph{k}-mer\cr
#'   \code{background_count} \tab background counts for each \emph{k}-mer\cr
#'   \code{enrichment} \tab \emph{k}-mer enrichment\cr
#'   \code{p_value} \tab p-value of \emph{k}-mer enrichment (either from
#'   Fisher's exact test or Pearson's chi-squared test)\cr
#'   \code{adj_p_value} \tab multiple testing corrected p-value\cr
#' }
#'
#' @details
#' Usually uses Pearson's chi-squared test, but recalculates p-values
#' with Fisher's exact test
#' for Pearson's chi-squared test p-values \code{<= chisq_p_value_threshold}.
#' The reason this is done is
#' computational efficiency. Fisher's exact tests are computationally
#' demanding and are only
#' performed in situations, where exact p-values are preferred, e.g.,
#' if expected < 5 or
#' significant p-values.
#'
#' @examples
#' # define simple sequence sets for foreground and background
#' foreground_set <- c(
#'   "CAACAGCCUUAAUU", "CAGUCAAGACUCC", "CUUUGGGGAAU",
#'   "UCAUUUUAUUAAA", "AAUUGGUGUCUGGAUACUUCCCUGUACAU",
#'   "AUCAAAUUA", "AGAU", "GACACUUAAAGAUCCU",
#'   "UAGCAUUAACUUAAUG", "AUGGA", "GAAGAGUGCUCA",
#'   "AUAGAC", "AGUUC", "CCAGUAA"
#' )
#' background_set <- c(
#'   "CAACAGCCUUAAUU", "CAGUCAAGACUCC", "CUUUGGGGAAU",
#'   "UCAUUUUAUUAAA", "AAUUGGUGUCUGGAUACUUCCCUGUACAU",
#'   "AUCAAAUUA", "AGAU", "GACACUUAAAGAUCCU",
#'   "UAGCAUUAACUUAAUG", "AUGGA", "GAAGAGUGCUCA",
#'   "AUAGAC", "AGUUC", "CCAGUAA",
#'   "UUAUUUA", "AUCCUUUACA", "UUUUUUU", "UUUCAUCAUU",
#'   "CCACACAC", "CUCAUUGGAG", "ACUUUGGGACA", "CAGGUCAGCA"
#' )
#' foreground_kmers <- generate_kmers(foreground_set, 6)
#' background_kmers <- generate_kmers(background_set, 6)
#'
#'
#' kmer_enrichment_values <- compute_kmer_enrichment(foreground_kmers,
#'   background_kmers)
#' @importFrom stats chisq.test
#' @importFrom stats p.adjust
#' @importFrom stats fisher.test
#' @family \emph{k}-mer functions
#' @export
compute_kmer_enrichment <-
    function(foreground_kmers,
             background_kmers,
             permutation = FALSE,
             chisq_p_value_threshold = 0.05,
             p_adjust_method = "BH") {
        background_kmers_sum <- sum(background_kmers)
        other_background_kmers <-
            background_kmers_sum - background_kmers

        foreground_kmers_sum <- sum(foreground_kmers)
        other_foreground_kmers <-
            foreground_kmers_sum - foreground_kmers

        # calculate enrichment value
        if (foreground_kmers_sum == 0 ||
            background_kmers_sum == 0) {
            enrichment <- as.numeric(NA)
        } else {
            enrichment <-
                (foreground_kmers / foreground_kmers_sum) / (background_kmers / background_kmers_sum)
        }

        if (permutation) {
            return(enrichment)
        } else {
            # Pearson's Chi-squared test
            chisq_p_values <-
                vapply(seq_len(length(foreground_kmers)), function(i) {
                    cont <- matrix(
                        c(
                            foreground_kmers[i],
                            other_foreground_kmers[i],
                            background_kmers[i],
                            other_background_kmers[i]
                        ),
                        nrow = 2
                    )
                    res <- suppressWarnings(stats::chisq.test(cont))
                    if (sum(res$expected < 5) > 0) {
                        # calculate p-value with Fisher's exact test
                        res$p.value <- 0
                    }
                    return(res$p.value)
                }, numeric(1))

            # calculate Fisher's exact test for Chi-squared test
            # p-values <= chisq_p_value_threshold
            idx <-
                which(
                    stats::p.adjust(chisq_p_values, method = p_adjust_method) <= chisq_p_value_threshold
                )

            if (length(idx) > 0) {
                # Fisher's exact test
                fisher_p_values <- vapply(idx, function(i) {
                    cont <- matrix(
                        c(
                            foreground_kmers[i],
                            other_foreground_kmers[i],
                            background_kmers[i],
                            other_background_kmers[i]
                        ),
                        nrow = 2
                    )
                    return(stats::fisher.test(cont)$p.value)
                }, numeric(1))

                # replace Chi-squared test p-values with Fisher's
                # exact test p-values
                chisq_p_values[idx] <- fisher_p_values
            }

            # adjust for multiple testing
            adj_p_value <-
                stats::p.adjust(chisq_p_values, method = p_adjust_method)
            df <- data.frame(
                foreground_count = foreground_kmers,
                background_count = background_kmers,
                enrichment = enrichment,
                p_value = chisq_p_values,
                adj_p_value = adj_p_value,
                stringsAsFactors = FALSE
            )
            rownames(df) <- NULL
            return(df)
        }
    }

#' @title Significance of Observed Mean
#'
#' @description
#' \code{estimate_significance_core} returns an estimate of the significance
#' of the observed
#' mean, given a vector of means based on random permutations of the data.
#'
#' @param random_means numeric vector of means based on random permutations
#' of the data
#' (empirical null distribution)
#' @param actual_mean observed mean
#' @param alternative side of the test, one of the following:
#' \code{"two_sided"},
#' \code{"less"}, \code{"greater"}
#' @param conf_level confidence level for the returned confidence interval
#'
#' @return A list with the following components:
#' \tabular{rl}{
#'   \code{p_value_estimate} \tab the estimated p-value of the observed mean\cr
#'   \code{conf_int} \tab the confidence interval around that estimate
#' }
#'
#' @examples
#' test_sd <- 1.0
#' test_null_distribution <- rnorm(n = 10000, mean = 1.0, sd = test_sd)
#'
#' estimate_significance_core(test_null_distribution, test_sd * 2, "greater")
#' @importFrom stats binom.test
#' @family \emph{k}-mer functions
#' @export
estimate_significance_core <- function(random_means,
                                       actual_mean,
                                       alternative = c("two_sided",
                                                       "less", "greater"),
                                       conf_level = 0.95) {
    alternative <- match.arg(alternative,
                             choices = c("two_sided", "less", "greater"))
    if (alternative == "greater") {
        # upper tail probability
        k <- sum(random_means >= actual_mean)
    } else if (alternative == "less") {
        # lower tail probability
        k <- sum(random_means <= actual_mean)
    } else {
        # two-tailed probability
        k <- sum(abs(log2(random_means)) >= abs(log2(actual_mean)))
    }
    p_value <- (k + 1) / (length(random_means) + 1)

    conf_int <- stats::binom.test(k,
                                  length(random_means),
                                  p = 0.5,
                                  conf.level = conf_level
    )$conf.int[seq_len(2)]
    return(list(p_value_estimate = p_value, conf_int = conf_int))
}

#' @title Geometric Mean
#'
#' @description
#' Calculates the geometric mean of the specified values.
#'
#' @param x numeric vector of values for which the geometric mean
#' will be computed
#' @param na_rm logical. Should missing values (including NaN) be removed?
#'
#' @return Geometric mean of \code{x} or \code{1} if length of \code{x} is 0
#'
#' @examples
#' geometric_mean(c(0.123, 0.441, 0.83))
#' @export
geometric_mean <- function(x, na_rm = TRUE) {
    if (length(x) > 0) {
        return(exp(sum(log(x), na.rm = na_rm) / length(x)))
    } else {
        return(1)
    }
}
#' @title Generate Random Permutations of the Enrichment Data
#'
#' @description
#' Calculates \emph{k}-mer enrichment values for randomly sampled (without
#' replacement) foreground sets.
#'
#' @param n_transcripts_foreground number of transcripts in the original
#' foreground set
#' @param n_permutations number of permutations to perform
#' @inheritParams calculate_kmer_enrichment
#'
#' @return The result of \code{\link{calculate_kmer_enrichment}} for the
#' random foreground sets.
#'
#' @family \emph{k}-mer functions
#' @export
generate_permuted_enrichments <-
    function(n_transcripts_foreground,
             background_set,
             k,
             n_permutations = 1000,
             n_cores = 4) {
        background_seq_n <- length(background_set)
        random_foreground_sets <-
            lapply(seq_len(n_permutations), function(i) {
                return(background_set[sample.int(background_seq_n,
                                                 n_transcripts_foreground)])
            })
        return(
            calculate_kmer_enrichment(
                random_foreground_sets,
                background_set,
                k,
                permutation = TRUE,
                n_cores = n_cores
            )
        )
    }

#' @title Permutation Test Based Significance of Observed Mean
#'
#' @description
#' \code{estimate_significance} returns an estimate of the significance
#' of the observed
#' mean, given a set of random permutations of the data.
#'
#'
#' @param motif_kmers set of \emph{k}-mers that were used to compute
#' the \code{actual_mean}
#' @param random_permutations a set of random permutations of the
#' original data, used
#' to generate an empirical null distribution.
#' @param produce_plot if distribution plot should be part of the
#' returned list
#' @inheritParams estimate_significance_core
#'
#' @return A list with the following components:
#' \tabular{rl}{
#'   \code{p_value_estimate} \tab the estimated p-value of the
#'   observed mean\cr
#'   \code{conf_int} \tab the confidence interval around that estimate\cr
#'   \code{plot} \tab plot of the empirical distribution of
#'   geometric means of the
#'   enrichment values
#' }
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_histogram
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_density
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 ggtitle
#' @family \emph{k}-mer functions
#' @export
estimate_significance <-
    function(actual_mean,
             motif_kmers,
             random_permutations,
             alternative = c("two_sided", "less", "greater"),
             conf_level = 0.95,
             produce_plot = TRUE) {
        # avoid CRAN note
        means <- actual <- ..density.. <- NULL

        idx <- which(random_permutations$kmers %in% motif_kmers)
        random_means <-
            unlist(lapply(random_permutations$dfs, function(enrichment_values) {
                return(geometric_mean(enrichment_values[idx]))
            }))

        if (produce_plot) {
            df <-
                data.frame(
                    means = random_means,
                    actual = rep(actual_mean, length(random_means))
                )
            bin_width <- (max(df$means) - min(df$means)) / 40
            if (bin_width <= 0) {
                bin_width <- 0.001
            }
            dist_plot <-
                ggplot2::ggplot(df, ggplot2::aes(x = means)) +
                ggplot2::geom_histogram(
                    ggplot2::aes(y = ..density..),
                    color = "black",
                    fill = "white",
                    binwidth = bin_width
                ) +
                ggplot2::geom_density(alpha = 0.2, fill = "#66FF66") +
                ggplot2::geom_vline(
                    ggplot2::aes(xintercept = mean(actual)),
                    color = "red",
                    linetype = "dashed",
                    size = 1
                ) +
                ggplot2::xlab("mean enrichment value") +
                ggplot2::theme_bw(base_size = 16) +
                ggplot2::theme(
                    axis.line = ggplot2::element_line(color = "black"),
                    panel.grid.major = ggplot2::element_blank(),
                    panel.grid.minor = ggplot2::element_blank(),
                    panel.border = ggplot2::element_blank(),
                    panel.background = ggplot2::element_blank()
                ) +
                ggplot2::ggtitle("Monte Carlo sampling of distribution of
                                 mean enrichment values")
        } else {
            dist_plot <- NULL
        }

        empirical_p_value <-
            estimate_significance_core(random_means, actual_mean,
                                       alternative = alternative, conf_level
            )
        return(
            list(
                p_value_estimate = empirical_p_value$p_value_estimate,
                conf_int = empirical_p_value$conf_int,
                plot = dist_plot
            )
        )
    }

#' @title \emph{k}-mer Enrichment Volcano Plot
#'
#' @description
#' Uses a volcano plot to visualize \emph{k}-mer enrichment.
#' X-axis is \eqn{\log_2} enrichment value,
#' y-axis is \eqn{\log_10} significance, i.e., multiple testing
#' corrected p-value from
#' Fisher's exact test or Pearson's chi-squared test.
#'
#' @param kmers data frame with the following columns: kmer,
#' adj_p_value, enrichment
#' @param motif_kmers set of \emph{k}-mers that are associated with
#' a certain motif,
#' will be highlighted in volcano plot
#' @param motif_rbps name of RNA-binding proteins associated with
#' highlighted \emph{k}-mers
#' (character vector of length 1)
#' @param significance_threshold p-value threshold for significance,
#' e.g., \code{0.05} or
#' \code{0.01}
#' @param show_legend whether or not a legend should be shown
#'
#' @return volcano plot
#'
#' @examples
#' motif <- get_motif_by_id("951_12324455")
#' draw_volcano_plot(transite:::kmers_enrichment, get_hexamers(motif[[1]]),
#'   get_rbps(motif[[1]]))
#'
#' \dontrun{
#' foreground_set <- c("UGUGGG", "GUGGGG", "GUGUGG", "UGUGGU")
#' background_set <- unique(c(foreground_set, c(
#'   "CAACAGCCUUAAUU", "CAGUCAAGACUCC", "CUUUGGGGAAU",
#'   "UCAUUUUAUUAAA", "AAUUGGUGUCUGGAUACUUCCCUGUACAU",
#'   "AUCAAAUUA", "AGAU", "GACACUUAAAGAUCCU",
#'   "UAGCAUUAACUUAAUG", "AUGGA", "GAAGAGUGCUCA",
#'   "AUAGAC", "AGUUC", "CCAGUAA",
#'   "CCACACAC", "CUCAUUGGAG", "ACUUUCCCACA", "CAGGUCAGCA",
#'   "CCACACCAG", "CCACACAUCAGU", "CACACACUCC", "CAGCCCCCCACAGGCA"
#' )))
#'
#' motif <- get_motif_by_id("M178_0.6")
#' results <- run_kmer_tsma(list(foreground_set), background_set,
#'                        motifs = motif)
#' draw_volcano_plot(results[[1]]$motif_kmers_dfs[[1]],
#'     get_hexamers(motif[[1]]), "test RBP")}
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_size_manual
#' @importFrom ggplot2 scale_alpha_manual
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom dplyr filter
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 element_blank
#' @family TSMA functions
#' @family \emph{k}-mer functions
#' @export
draw_volcano_plot <-
    function(kmers,
             motif_kmers,
             motif_rbps,
             significance_threshold = 0.01,
             show_legend = TRUE) {
        # avoid CRAN note
        log_enrichment <-
            log_adj_p_value <- group <- motifidx <- NULL

        indices <- as.character(kmers$kmer) %in% motif_kmers

        kmers$sig <- kmers$adj_p_value < significance_threshold
        kmers$log_adj_p_value <-
            -log10(as.numeric(kmers$adj_p_value))
        kmers$enrichment[is.na(kmers$enrichment)] <- 1
        kmers$log_enrichment <- log2(as.numeric(kmers$enrichment))
        kmers$group[kmers$adj_p_value >= significance_threshold] <-
            "n"
        kmers$group[kmers$enrichment > 1 &
                        kmers$adj_p_value <= significance_threshold] <-
            "e"
        kmers$group[kmers$enrichment < 1 &
                        kmers$adj_p_value <= significance_threshold] <-
            "d"
        kmers$group[indices == TRUE] <- "m"
        kmers$motifidx <- indices

        volcano_plot <-
            ggplot2::ggplot(
                kmers,
                ggplot2::aes(
                    x = log_enrichment,
                    y = log_adj_p_value,
                    color = group,
                    alpha = group,
                    size = motifidx
                )
            ) +
            ggplot2::geom_point() +
            ggplot2::scale_color_manual(
                values = c(
                    "n" = "#000000",
                    "m" = "yellow",
                    "e" = "#832424FF",
                    "d" = "#004F91FF"
                ),
                name = "",
                breaks = c("d", "e", "m", "n"),
                labels = c(
                    "depleted",
                    "enriched",
                    paste0(motif_rbps, " motif kmer"),
                    "not significant"
                )
            ) +
            ggplot2::scale_size_manual(values = c(4, 8), guide = FALSE) +
            ggplot2::scale_alpha_manual(values = c(0.5, 0.5, 0.7, 0.5),
                                        guide = FALSE) +
            ggplot2::geom_hline(yintercept = -log10(significance_threshold)) +
            ggplot2::ylab(expression(-log[10](q))) +
            ggplot2::xlab(expression(log[2](FC))) +
            ggplot2::geom_point(
                data = dplyr::filter(kmers, group == "m"),
                ggplot2::aes(log_enrichment, log_adj_p_value),
                shape = 21,
                color = "#000000",
                fill = "yellow"
            ) +
            ggplot2::theme_bw(base_size = 16) +
            ggplot2::theme(
                axis.line = ggplot2::element_line(colour = "black"),
                panel.grid.major = ggplot2::element_blank(),
                panel.grid.minor = ggplot2::element_blank(),
                panel.border = ggplot2::element_blank(),
                panel.background = ggplot2::element_blank(),
                legend.title = ggplot2::element_blank(),
                legend.position = "top"
            )
        if (!show_legend) {
            volcano_plot <-
                volcano_plot + ggplot2::theme(legend.position = "none")
        }

        return(volcano_plot)
    }

#' @title Check Validity of Set of \emph{k}-mers
#'
#' @description
#' Checks if the provided set of \emph{k}-mers is valid. A valid set
#' of \emph{k}-mers is
#' (1) non-empty, (2) contains either only hexamers or only heptamers,
#'  and (3) contains
#' only characters from the RNA alphabet (A, C, G, U)
#'
#' @param kmers set of \emph{k}-mers
#'
#' @return \code{TRUE} if set of \emph{k}-mers is valid
#'
#' @examples
#' # valid set
#' check_kmers(c("ACGCUC", "AAACCC", "UUUACA"))
#'
#' # invalid set (contains hexamers and heptamers)
#' check_kmers(c("ACGCUC", "AAACCC", "UUUACAA"))
#' @family \emph{k}-mer functions
#' @export
check_kmers <- function(kmers) {
    if (length(kmers) == 0) {
        return(FALSE)
    }
    # all k-mers have to be of the same length
    if (length(unique(nchar(kmers))) != 1) {
        return(FALSE)
    }
    # only hexamers or heptamers allowed
    if (!(unique(nchar(kmers)) %in% c(6, 7))) {
        return(FALSE)
    }
    invalid_characters <-
        gsub("A", "", gsub("C", "", gsub("G", "", gsub("U", "", kmers))))
    if (sum(nchar(invalid_characters)) > 0) {
        return(FALSE)
    }
    return(TRUE)
}
