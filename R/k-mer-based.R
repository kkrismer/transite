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
#' @param is.rna if \code{sequences} are RNA sequences, this
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
homopolymerCorrection <-
    function(sequences, k, kmers, is.rna = FALSE) {
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

        if (is.rna) {
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
#' Calls \code{\link{computeKmerEnrichment}} to compute \emph{k}-mer
#' enrichment values
#' for multiple foregrounds. Calculates enrichment for foreground sets in
#' parallel.
#'
#' @param foreground.sets list of foreground sets; a foreground set is a
#' character vector of
#' DNA or RNA sequences (not both) and a strict subset of the
#' \code{background.set}
#' @param background.set character vector of DNA or RNA sequences that
#' constitute the
#' background set
#' @param n.cores number of computing cores to use
#' @inheritParams generateKmers
#' @inheritParams computeKmerEnrichment
#'
#' @return A list with two entries:
#'
#' (1) dfs: a list of data frames with results from
#' \code{\link{computeKmerEnrichment}} for each of the foreground sets
#' (2) kmers: a character vector of all k-mers
#'
#' @examples
#' # define simple sequence sets for foreground and background
#' foreground.set1 <- c(
#'   "CAACAGCCUUAAUU", "CAGUCAAGACUCC", "CUUUGGGGAAU",
#'   "UCAUUUUAUUAAA", "AAUUGGUGUCUGGAUACUUCCCUGUACAU",
#'   "AUCAAAUUA", "AGAU", "GACACUUAAAGAUCCU",
#'   "UAGCAUUAACUUAAUG", "AUGGA", "GAAGAGUGCUCA",
#'   "AUAGAC", "AGUUC", "CCAGUAA"
#' )
#' foreground.set2 <- c("UUAUUUA", "AUCCUUUACA", "UUUUUUU", "UUUCAUCAUU")
#' foreground.sets <- list(foreground.set1, foreground.set2)
#' background.set <- c(foreground.set1, foreground.set2,
#'                     "CCACACAC", "CUCAUUGGAG", "ACUUUGGGACA", "CAGGUCAGCA")
#'
#' kmer.enrichment.values <- calculateKmerEnrichment(foreground.sets,
#'   background.set, 6)
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel parLapply
#' @family \emph{k}-mer functions
#' @export
calculateKmerEnrichment <-
    function(foreground.sets,
             background.set,
             k,
             permutation = FALSE,
             chisq.p.value.threshold = 0.05,
             p.adjust.method = "BH",
             n.cores = 4) {
        if (length(foreground.sets) > 0) {
            background.kmers <- generateKmers(background.set, k)

            n.cores <- min(n.cores, length(foreground.sets))
            if (n.cores == 1) {
                enrichment.dfs <- lapply(foreground.sets,
                                         function(foreground.set) {
                    foreground.kmers <- generateKmers(foreground.set, k)
                    return(
                        computeKmerEnrichment(
                            foreground.kmers,
                            background.kmers,
                            permutation = permutation,
                            chisq.p.value.threshold = chisq.p.value.threshold,
                            p.adjust.method = p.adjust.method
                        )
                    )
                })
            } else {
                cluster <-
                    parallel::makeCluster(mc <-
                                              getOption("cl.cores", n.cores))
                parallel::clusterExport(
                    cl = cluster,
                    varlist = c(
                        "k",
                        "background.kmers",
                        "permutation",
                        "chisq.p.value.threshold",
                        "p.adjust.method",
                        "computeKmerEnrichment",
                        "generateKmers"
                    ),
                    envir = environment()
                )
                enrichment.dfs <-
                    parallel::parLapply(
                        cl = cluster, foreground.sets,
                        function(foreground.set) {
                            foreground.kmers <- generateKmers(foreground.set, k)
                            return(
                                computeKmerEnrichment(
                                    foreground.kmers,
                                    background.kmers,
                                    permutation = permutation,
                                    chisq.p.value.threshold = chisq.p.value.threshold,
                                    p.adjust.method = p.adjust.method
                                )
                            )
                        }
                    )
            }
            return(list(
                dfs = enrichment.dfs,
                kmers = gsub("T", "U", as.character(names(
                    background.kmers
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
#' @inheritParams homopolymerCorrection
#'
#' @return Returns a named numeric vector, where the elements are
#' \emph{k}-mer counts and the
#' names are DNA \emph{k}-mers.
#'
#' @examples
#' # count hexamers in set of RNA sequences
#' rna.sequences <- c(
#'   "CAACAGCCUUAAUU", "CAGUCAAGACUCC", "CUUUGGGGAAU",
#'   "UCAUUUUAUUAAA", "AAUUGGUGUCUGGAUACUUCCCUGUACAU",
#'   "AUCAAAUUA", "AGAU", "GACACUUAAAGAUCCU",
#'   "UAGCAUUAACUUAAUG", "AUGGA", "GAAGAGUGCUCA",
#'   "AUAGAC", "AGUUC", "CCAGUAA",
#'   "UUAUUUA", "AUCCUUUACA", "UUUUUUU", "UUUCAUCAUU",
#'   "CCACACAC", "CUCAUUGGAG", "ACUUUGGGACA", "CAGGUCAGCA"
#' )
#' hexamer.counts <- generateKmers(rna.sequences, 6)
#'
#'
#' # count heptamers in set of DNA sequences
#' dna.sequences <- c(
#'   "CAACAGCCTTAATT", "CAGTCAAGACTCC", "CTTTGGGGAAT",
#'   "TCATTTTATTAAA", "AATTGGTGTCTGGATACTTCCCTGTACAT",
#'   "ATCAAATTA", "AGAT", "GACACTTAAAGATCCT",
#'   "TAGCATTAACTTAATG", "ATGGA", "GAAGAGTGCTCA",
#'   "ATAGAC", "AGTTC", "CCAGTAA",
#'   "TTATTTA", "ATCCTTTACA", "TTTTTTT", "TTTCATCATT",
#'   "CCACACAC", "CTCATTGGAG", "ACTTTGGGACA", "CAGGTCAGCA"
#' )
#' hexamer.counts <- generateKmers(dna.sequences, 7)
#' @section Warning:
#' \code{generateKmers} always returns DNA \emph{k}-mers, even if
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
generateKmers <- function(sequences, k) {
    kmers <- tryCatch({
        kmer.counts <-
            Biostrings::oligonucleotideFrequency(Biostrings::DNAStringSet(sequences), k)
        homopolymerCorrection(sequences, k, colSums(kmer.counts))
    }, error = function(e) {
        kmer.counts <-
            Biostrings::oligonucleotideFrequency(Biostrings::RNAStringSet(sequences), k)
        rna.kmers <-
            homopolymerCorrection(sequences, k, colSums(kmer.counts), is.rna = TRUE)
        names(rna.kmers) <- gsub("U", "T", names(rna.kmers))
        return(rna.kmers)
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
#' @param foreground.kmers \emph{k}-mer counts of the foreground set
#' (generated by \code{\link{generateKmers}})
#' @param background.kmers \emph{k}-mer counts of the background set
#' (generated by \code{\link{generateKmers}})
#' @param permutation if \code{TRUE}, only the enrichment value is returned
#' (efficiency mode
#' used for permutation testing)
#' @param chisq.p.value.threshold threshold below which Fisher's exact test
#' is used instead of Pearson's chi-squared test
#' @param p.adjust.method see \code{\link[stats]{p.adjust}}
#'
#' @return enrichment of \emph{k}-mers in specified foreground sequences.
#' A data frame with the following columns is returned:
#' \tabular{rl}{
#'   \code{foreground.count} \tab foreground counts for each \emph{k}-mer\cr
#'   \code{background.count} \tab background counts for each \emph{k}-mer\cr
#'   \code{enrichment} \tab \emph{k}-mer enrichment\cr
#'   \code{p.value} \tab p-value of \emph{k}-mer enrichment (either from
#'   Fisher's exact test or Pearson's chi-squared test.\cr
#'   \code{adj.p.value} \tab multiple testing corrected p-value\cr
#' }
#'
#' @details
#' Usually uses Pearson's chi-squared test, but recalculates p-values
#' with Fisher's exact test
#' for Pearson's chi-squared test p-values \code{<= chisq.p.value.threshold}.
#' The reason this is done is
#' computational efficiency. Fisher's exact tests are computationally
#' demanding and are only
#' performed in situations, where exact p-values are preferred, e.g.,
#' if expected < 5 or
#' significant p-values.
#'
#' @examples
#' # define simple sequence sets for foreground and background
#' foreground.set <- c(
#'   "CAACAGCCUUAAUU", "CAGUCAAGACUCC", "CUUUGGGGAAU",
#'   "UCAUUUUAUUAAA", "AAUUGGUGUCUGGAUACUUCCCUGUACAU",
#'   "AUCAAAUUA", "AGAU", "GACACUUAAAGAUCCU",
#'   "UAGCAUUAACUUAAUG", "AUGGA", "GAAGAGUGCUCA",
#'   "AUAGAC", "AGUUC", "CCAGUAA"
#' )
#' background.set <- c(
#'   "CAACAGCCUUAAUU", "CAGUCAAGACUCC", "CUUUGGGGAAU",
#'   "UCAUUUUAUUAAA", "AAUUGGUGUCUGGAUACUUCCCUGUACAU",
#'   "AUCAAAUUA", "AGAU", "GACACUUAAAGAUCCU",
#'   "UAGCAUUAACUUAAUG", "AUGGA", "GAAGAGUGCUCA",
#'   "AUAGAC", "AGUUC", "CCAGUAA",
#'   "UUAUUUA", "AUCCUUUACA", "UUUUUUU", "UUUCAUCAUU",
#'   "CCACACAC", "CUCAUUGGAG", "ACUUUGGGACA", "CAGGUCAGCA"
#' )
#' foreground.kmers <- generateKmers(foreground.set, 6)
#' background.kmers <- generateKmers(background.set, 6)
#'
#'
#' kmer.enrichment.values <- computeKmerEnrichment(foreground.kmers,
#'   background.kmers)
#' @importFrom stats chisq.test
#' @importFrom stats p.adjust
#' @importFrom stats fisher.test
#' @family \emph{k}-mer functions
#' @export
computeKmerEnrichment <-
    function(foreground.kmers,
             background.kmers,
             permutation = FALSE,
             chisq.p.value.threshold = 0.05,
             p.adjust.method = "BH") {
        background.kmers.sum <- sum(background.kmers)
        other.background.kmers <-
            background.kmers.sum - background.kmers

        foreground.kmers.sum <- sum(foreground.kmers)
        other.foreground.kmers <-
            foreground.kmers.sum - foreground.kmers

        # calculate enrichment value
        if (foreground.kmers.sum == 0 ||
            background.kmers.sum == 0) {
            enrichment <- as.numeric(NA)
        } else {
            enrichment <-
                (foreground.kmers / foreground.kmers.sum) / (background.kmers / background.kmers.sum)
        }

        if (permutation) {
            return(enrichment)
        } else {
            # Pearson's Chi-squared test
            chisq.p.values <-
                vapply(seq_len(length(foreground.kmers)), function(i) {
                    cont <- matrix(
                        c(
                            foreground.kmers[i],
                            other.foreground.kmers[i],
                            background.kmers[i],
                            other.background.kmers[i]
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
            # p-values <= chisq.p.value.threshold
            idx <-
                which(
                    stats::p.adjust(chisq.p.values, method = p.adjust.method) <= chisq.p.value.threshold
                )

            if (length(idx) > 0) {
                # Fisher's exact test
                fisher.p.values <- vapply(idx, function(i) {
                    cont <- matrix(
                        c(
                            foreground.kmers[i],
                            other.foreground.kmers[i],
                            background.kmers[i],
                            other.background.kmers[i]
                        ),
                        nrow = 2
                    )
                    return(stats::fisher.test(cont)$p.value)
                }, numeric(1))

                # replace Chi-squared test p-values with Fisher's
                # exact test p-values
                chisq.p.values[idx] <- fisher.p.values
            }

            # adjust for multiple testing
            adj.p.value <-
                stats::p.adjust(chisq.p.values, method = p.adjust.method)
            df <- data.frame(
                foreground.count = foreground.kmers,
                background.count = background.kmers,
                enrichment = enrichment,
                p.value = chisq.p.values,
                adj.p.value = adj.p.value,
                stringsAsFactors = FALSE
            )
            rownames(df) <- NULL
            return(df)
        }
    }

#' @title Significance of Observed Mean
#'
#' @description
#' \code{empiricalEnrichmentMeanCDF} returns an estimate of the significance
#' of the observed
#' mean, given a vector of means based on random permutations of the data.
#'
#' @param random.means numeric vector of means based on random permutations
#' of the data
#' (empirical null distribution)
#' @param actual.mean observed mean
#' @param alternative side of the test, one of the following:
#' \code{"two.sided"},
#' \code{"less"}, \code{"greater"}
#' @inheritParams stats::binom.test
#'
#' @return A list with the following components:
#' \tabular{rl}{
#'   \code{p.value.estimate} \tab the estimated p-value of the observed mean\cr
#'   \code{conf.int} \tab the confidence interval around that estimate
#' }
#'
#' @examples
#' test.sd <- 1.0
#' test.null.distribution <- rnorm(n = 10000, mean = 1.0, sd = test.sd)
#'
#' empiricalEnrichmentMeanCDF(test.null.distribution, test.sd * 2, "greater")
#' @importFrom stats binom.test
#' @family \emph{k}-mer functions
#' @export
empiricalEnrichmentMeanCDF <- function(random.means,
                                       actual.mean,
                                       alternative = c("two.sided",
                                                       "less", "greater"),
                                       conf.level = 0.95) {
    alternative <- match.arg(alternative,
                             choices = c("two.sided", "less", "greater"))
    if (alternative == "greater") {
        # upper tail probability
        k <- sum(random.means >= actual.mean)
    } else if (alternative == "less") {
        # lower tail probability
        k <- sum(random.means <= actual.mean)
    } else {
        # two-tailed probability
        k <- sum(abs(log2(random.means)) >= abs(log2(actual.mean)))
    }
    p.value <- (k + 1) / (length(random.means) + 1)

    conf.int <- stats::binom.test(k,
                                  length(random.means),
                                  p = 0.5,
                                  conf.level = conf.level
    )$conf.int[seq_len(2)]
    return(list(p.value.estimate = p.value, conf.int = conf.int))
}

#' @title Geometric Mean
#'
#' @description
#' Calculates the geometric mean of the specified values.
#'
#' @param x numeric vector of values for which the geometric mean
#' will be computed
#' @inheritParams base::sum
#'
#' @return Geometric mean of \code{x} or \code{1} if length of \code{x} is 0
#'
#' @examples
#' geometricMean(c(0.123, 0.441, 0.83))
#' @export
geometricMean <- function(x, na.rm = TRUE) {
    if (length(x) > 0) {
        return(exp(sum(log(x), na.rm = na.rm) / length(x)))
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
#' @param n.transcripts.foreground number of transcripts in the original
#' foreground set
#' @param n.permutations number of permutations to perform
#' @inheritParams calculateKmerEnrichment
#'
#' @return The result of \code{\link{calculateKmerEnrichment}} for the
#' random foreground sets.
#'
#' @family \emph{k}-mer functions
#' @export
generatePermutedEnrichments <-
    function(n.transcripts.foreground,
             background.set,
             k,
             n.permutations = 1000,
             n.cores = 4) {
        background.seq.n <- length(background.set)
        random.foreground.sets <-
            lapply(seq_len(n.permutations), function(i) {
                return(background.set[sample.int(background.seq.n,
                                                 n.transcripts.foreground)])
            })
        return(
            calculateKmerEnrichment(
                random.foreground.sets,
                background.set,
                k,
                permutation = TRUE,
                n.cores = n.cores
            )
        )
    }

#' @title Permutation Test Based Significance of Observed Mean
#'
#' @description
#' \code{permTestGeometricMean} returns an estimate of the significance
#' of the observed
#' mean, given a set of random permutations of the data.
#'
#'
#' @param motif.kmers set of \emph{k}-mers that were used to compute
#' the \code{actual.mean}
#' @param random.permutations a set of random permutations of the
#' original data, used
#' to generate an empirical null distribution.
#' @param produce.plot if distribution plot should be part of the
#' returned list
#' @inheritParams empiricalEnrichmentMeanCDF
#'
#' @return A list with the following components:
#' \tabular{rl}{
#'   \code{p.value.estimate} \tab the estimated p-value of the
#'   observed mean\cr
#'   \code{conf.int} \tab the confidence interval around that estimate\cr
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
permTestGeometricMean <-
    function(actual.mean,
             motif.kmers,
             random.permutations,
             alternative = c("two.sided", "less", "greater"),
             conf.level = 0.95,
             produce.plot = TRUE) {
        # avoid CRAN note
        means <- actual <- ..density.. <- NULL

        idx <- which(random.permutations$kmers %in% motif.kmers)
        random.means <-
            unlist(lapply(random.permutations$dfs, function(enrichment.values) {
                return(geometricMean(enrichment.values[idx]))
            }))

        if (produce.plot) {
            df <-
                data.frame(
                    means = random.means,
                    actual = rep(actual.mean, length(random.means))
                )
            bin.width <- (max(df$means) - min(df$means)) / 40
            if (bin.width <= 0) {
                bin.width <- 0.001
            }
            dist.plot <-
                ggplot2::ggplot(df, ggplot2::aes(x = means)) +
                ggplot2::geom_histogram(
                    ggplot2::aes(y = ..density..),
                    colour = "black",
                    fill = "white",
                    binwidth = bin.width
                ) +
                ggplot2::geom_density(alpha = .2, fill = "#66FF66") +
                ggplot2::geom_vline(
                    ggplot2::aes(xintercept = mean(actual)),
                    color = "red",
                    linetype = "dashed",
                    size = 1
                ) +
                ggplot2::xlab("mean enrichment value") +
                ggplot2::theme_bw(base_size = 16) +
                ggplot2::theme(
                    axis.line = ggplot2::element_line(colour = "black"),
                    panel.grid.major = ggplot2::element_blank(),
                    panel.grid.minor = ggplot2::element_blank(),
                    panel.border = ggplot2::element_blank(),
                    panel.background = ggplot2::element_blank()
                ) +
                ggplot2::ggtitle("Monte Carlo sampling of distribution of
                                 mean enrichment values")
        } else {
            dist.plot <- NULL
        }

        empirical.p.value <-
            empiricalEnrichmentMeanCDF(random.means, actual.mean,
                                       alternative = alternative, conf.level
            )
        return(
            list(
                p.value.estimate = empirical.p.value$p.value.estimate,
                conf.int = empirical.p.value$conf.int,
                plot = dist.plot
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
#' adj.p.value, enrichment
#' @param motif.kmers set of \emph{k}-mers that are associated with
#' a certain motif,
#' will be highlighted in volcano plot
#' @param motif.rbps name of RNA-binding proteins associated with
#' highlighted \emph{k}-mers
#' (character vector of length 1)
#' @param significance.threshold p-value threshold for significance,
#' e.g., \code{0.05} or
#' \code{0.01}
#' @param show.legend whether or not a legend should be shown
#'
#' @return volcano plot
#'
#' @examples
#' motif <- getMotifById("951_12324455")
#' drawVolcanoPlot(transite:::kmers.enrichment, motifHexamers(motif[[1]]),
#'   motifRbps(motif[[1]]))
#'
#' \dontrun{
#' foreground.set <- c("UGUGGG", "GUGGGG", "GUGUGG", "UGUGGU")
#' background.set <- unique(c(foreground.set, c(
#'   "CAACAGCCUUAAUU", "CAGUCAAGACUCC", "CUUUGGGGAAU",
#'   "UCAUUUUAUUAAA", "AAUUGGUGUCUGGAUACUUCCCUGUACAU",
#'   "AUCAAAUUA", "AGAU", "GACACUUAAAGAUCCU",
#'   "UAGCAUUAACUUAAUG", "AUGGA", "GAAGAGUGCUCA",
#'   "AUAGAC", "AGUUC", "CCAGUAA",
#'   "CCACACAC", "CUCAUUGGAG", "ACUUUCCCACA", "CAGGUCAGCA",
#'   "CCACACCAG", "CCACACAUCAGU", "CACACACUCC", "CAGCCCCCCACAGGCA"
#' )))
#'
#' motif <- getMotifById("M178_0.6")
#' results <- runKmerTSMA(list(foreground.set), background.set,
#'                        motifs = motif)
#' drawVolcanoPlot(results[[1]]$motif.kmers.dfs[[1]],
#'     motifHexamers(motif[[1]]), "test RBP")}
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
drawVolcanoPlot <-
    function(kmers,
             motif.kmers,
             motif.rbps,
             significance.threshold = 0.01,
             show.legend = TRUE) {
        # avoid CRAN note
        log.enrichment <-
            log.adj.p.value <- group <- motifidx <- NULL

        indices <- as.character(kmers$kmer) %in% motif.kmers

        kmers$sig <- kmers$adj.p.value < significance.threshold
        kmers$log.adj.p.value <-
            -log10(as.numeric(kmers$adj.p.value))
        kmers$enrichment[is.na(kmers$enrichment)] <- 1
        kmers$log.enrichment <- log2(as.numeric(kmers$enrichment))
        kmers$group[kmers$adj.p.value >= significance.threshold] <-
            "n"
        kmers$group[kmers$enrichment > 1 &
                        kmers$adj.p.value <= significance.threshold] <-
            "e"
        kmers$group[kmers$enrichment < 1 &
                        kmers$adj.p.value <= significance.threshold] <-
            "d"
        kmers$group[indices == TRUE] <- "m"
        kmers$motifidx <- indices

        volcano.plot <-
            ggplot2::ggplot(
                kmers,
                ggplot2::aes(
                    log.enrichment,
                    log.adj.p.value,
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
                    paste0(motif.rbps, " motif kmer"),
                    "not significant"
                )
            ) +
            ggplot2::scale_size_manual(values = c(4, 8), guide = FALSE) +
            ggplot2::scale_alpha_manual(values = c(0.5, 0.5, 0.7, 0.5),
                                        guide = FALSE) +
            ggplot2::geom_hline(yintercept = -log10(significance.threshold)) +
            ggplot2::ylab(expression(-log[10](q))) +
            ggplot2::xlab(expression(log[2](FC))) +
            ggplot2::geom_point(
                data = dplyr::filter(kmers, group == "m"),
                ggplot2::aes(log.enrichment, log.adj.p.value),
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
        if (!show.legend) {
            volcano.plot <-
                volcano.plot + ggplot2::theme(legend.position = "none")
        }

        return(volcano.plot)
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
#' checkKmers(c("ACGCUC", "AAACCC", "UUUACA"))
#'
#' # invalid set (contains hexamers and heptamers)
#' checkKmers(c("ACGCUC", "AAACCC", "UUUACAA"))
#' @family \emph{k}-mer functions
#' @export
checkKmers <- function(kmers) {
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
    invalid.characters <-
        gsub("A", "", gsub("C", "", gsub("G", "", gsub("U", "", kmers))))
    if (sum(nchar(invalid.characters)) > 0) {
        return(FALSE)
    }
    return(TRUE)
}
