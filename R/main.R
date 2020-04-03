#' @title Matrix-based Transcript Set Motif Analysis
#'
#' @description
#' Calculates motif enrichment in foreground sets versus a background
#' set using position
#' weight matrices to identify putative binding sites
#'
#' @param foreground_sets a list of named character vectors of
#' foreground sequences
#' (only containing upper case characters A, C, G, T), where the
#' names are RefSeq identifiers
#' and sequence
#' type qualifiers (\code{"3UTR"}, \code{"5UTR"}, \code{"mRNA"}), e.g.
#' \code{"NM_010356|3UTR"}. Names are only used to cache results.
#' @param background_set a named character vector of background
#' sequences (naming follows same
#' rules as foreground set sequences)
#' @inheritParams score_transcripts
#' @inheritParams calculate_motif_enrichment
#'
#' @return A list with the following components:
#' \tabular{rl}{
#'   \code{foreground_scores} \tab the result of \code{\link{score_transcripts}}
#'   for the foreground
#'   sets\cr
#'   \code{background_scores} \tab the result of \code{\link{score_transcripts}}
#'   for the background
#'   set\cr
#'   \code{enrichment_dfs} \tab a list of data frames, returned by
#'   \code{\link{calculate_motif_enrichment}}
#' }
#'
#' @details
#' Motif transcript set analysis can be used to identify RNA binding proteins,
#' whose targets are
#' significantly overrepresented or underrepresented in certain sets of
#' transcripts.
#'
#' The aim of Transcript Set Motif Analysis (TSMA) is to identify the
#' overrepresentation
#' and underrepresentation of potential RBP targets (binding sites)
#' in a set (or sets) of
#' sequences, i.e., the foreground set, relative to the entire population
#' of sequences.
#' The latter is called background set, which can be composed of all
#' sequences of the genes
#' of a microarray platform or all sequences of an organism or any
#' other meaningful
#' superset of the foreground sets.
#'
#' The matrix-based approach skips the \emph{k}-merization step of
#' the \emph{k}-mer-based approach
#' and instead scores the transcript sequence as a whole with a
#' position specific scoring matrix.
#'
#' For each sequence in foreground and background sets and each
#' sequence motif,
#' the scoring algorithm evaluates the score for each sequence position.
#' Positions with
#' a relative score greater than a certain threshold are considered hits, i.e.,
#' putative binding sites.
#'
#' By scoring all sequences in foreground and background sets, a hit count
#' for each motif and
#' each set is obtained, which is used to calculate enrichment values and
#' associated p-values
#' in the same way in which motif-compatible hexamer enrichment values are
#' calculated in
#' the k -mer-based approach. P-values are adjusted with one of the available
#' adjustment methods.
#'
#' An advantage of the matrix-based approach is the possibility of detecting
#' clusters of
#' binding sites. This can be done by counting regions with many hits using
#' positional
#' hit information or by simply applying a hit count threshold per sequence,
#' e.g., only
#' sequences with more than some number of hits are considered. Homotypic
#' clusters of RBP
#' binding sites may play a similar role as clusters of transcription factors.
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
#' names(foreground_set1) <- c(
#'   "NM_1_DUMMY|3UTR", "NM_2_DUMMY|3UTR", "NM_3_DUMMY|3UTR",
#'   "NM_4_DUMMY|3UTR", "NM_5_DUMMY|3UTR", "NM_6_DUMMY|3UTR",
#'   "NM_7_DUMMY|3UTR",
#'   "NM_8_DUMMY|3UTR", "NM_9_DUMMY|3UTR", "NM_10_DUMMY|3UTR",
#'   "NM_11_DUMMY|3UTR",
#'   "NM_12_DUMMY|3UTR", "NM_13_DUMMY|3UTR", "NM_14_DUMMY|3UTR"
#' )
#'
#' foreground_set2 <- c("UUAUUUA", "AUCCUUUACA", "UUUUUUU", "UUUCAUCAUU")
#' names(foreground_set2) <- c(
#'   "NM_15_DUMMY|3UTR", "NM_16_DUMMY|3UTR", "NM_17_DUMMY|3UTR",
#'   "NM_18_DUMMY|3UTR"
#' )
#'
#' foreground_sets <- list(foreground_set1, foreground_set2)
#'
#' background_set <- c(
#'   "CAACAGCCUUAAUU", "CAGUCAAGACUCC", "CUUUGGGGAAU",
#'   "UCAUUUUAUUAAA", "AAUUGGUGUCUGGAUACUUCCCUGUACAU",
#'   "AUCAAAUUA", "AGAU", "GACACUUAAAGAUCCU",
#'   "UAGCAUUAACUUAAUG", "AUGGA", "GAAGAGUGCUCA",
#'   "AUAGAC", "AGUUC", "CCAGUAA",
#'   "UUAUUUA", "AUCCUUUACA", "UUUUUUU", "UUUCAUCAUU",
#'   "CCACACAC", "CUCAUUGGAG", "ACUUUGGGACA", "CAGGUCAGCA"
#' )
#' names(background_set) <- c(
#'   "NM_1_DUMMY|3UTR", "NM_2_DUMMY|3UTR", "NM_3_DUMMY|3UTR",
#'   "NM_4_DUMMY|3UTR", "NM_5_DUMMY|3UTR", "NM_6_DUMMY|3UTR",
#'   "NM_7_DUMMY|3UTR",
#'   "NM_8_DUMMY|3UTR", "NM_9_DUMMY|3UTR", "NM_10_DUMMY|3UTR",
#'   "NM_11_DUMMY|3UTR",
#'   "NM_12_DUMMY|3UTR", "NM_13_DUMMY|3UTR", "NM_14_DUMMY|3UTR",
#'   "NM_15_DUMMY|3UTR",
#'   "NM_16_DUMMY|3UTR", "NM_17_DUMMY|3UTR", "NM_18_DUMMY|3UTR",
#'   "NM_19_DUMMY|3UTR",
#'   "NM_20_DUMMY|3UTR", "NM_21_DUMMY|3UTR", "NM_22_DUMMY|3UTR"
#' )
#'
#' # run cached version of TSMA with all Transite motifs (recommended):
#' # results <- run_matrix_tsma(foreground_sets, background_set)
#'
#' # run uncached version with one motif:
#' motif_db <- get_motif_by_id("M178_0.6")
#' results <- run_matrix_tsma(foreground_sets, background_set, motifs = motif_db,
#' cache = FALSE)
#'
#' \dontrun{
#' # define example sequence sets for foreground and background
#' foreground1_df <- transite:::ge$foreground1_df
#' foreground_set1 <- gsub("T", "U", foreground1_df$seq)
#' names(foreground_set1) <- paste0(foreground1_df$refseq, "|",
#'   foreground1_df$seq_type)
#'
#' foreground2_df <- transite:::ge$foreground2_df
#' foreground_set2 <- gsub("T", "U", foreground2_df$seq)
#' names(foreground_set2) <- paste0(foreground2_df$refseq, "|",
#'   foreground2_df$seq_type)
#'
#' foreground_sets <- list(foreground_set1, foreground_set2)
#'
#' background_df <- transite:::ge$background_df
#' background_set <- gsub("T", "U", background_df$seq)
#' names(background_set) <- paste0(background_df$refseq, "|",
#'   background_df$seq_type)
#'
#' # run cached version of TSMA with all Transite motifs (recommended)
#' results <- run_matrix_tsma(foreground_sets, background_set)
#'
#' # run uncached version of TSMA with all Transite motifs
#' results <- run_matrix_tsma(foreground_sets, background_set, cache = FALSE)
#'
#' # run TSMA with a subset of Transite motifs
#' results <- run_matrix_tsma(foreground_sets, background_set,
#'   motifs = get_motif_by_rbp("ELAVL1"))
#'
#' # run TSMA with user-defined motif
#' toy_motif <- create_matrix_motif(
#'   "toy_motif", "example RBP", toy_motif_matrix,
#'   "example type", "example species", "user"
#' )
#' results <- run_matrix_tsma(foreground_sets, background_set,
#'   motifs = list(toy_motif))
#' }
#'
#' @family TSMA functions
#' @family matrix functions
#' @importFrom dplyr arrange
#' @export
run_matrix_tsma <- function(foreground_sets,
                            background_set,
                            motifs = NULL,
                            max_hits = 5,
                            threshold_method = "p_value",
                            threshold_value = 0.25^6,
                            max_fg_permutations = 1000000,
                            min_fg_permutations = 1000,
                            e = 5,
                            p_adjust_method = "BH",
                            n_cores = 1,
                            cache = paste0(tempdir(), "/sc/")) {
    # avoid CRAN note
    adj_p_value <- p_value <- NULL

    i <- 1
    foreground_scores <-
        lapply(foreground_sets, function(foreground_set) {
            result <-
                score_transcripts(
                    foreground_set,
                    motifs = motifs,
                    max_hits = max_hits,
                    threshold_method = threshold_method,
                    threshold_value = threshold_value,
                    n_cores = n_cores,
                    cache = cache
                )
            message(paste0("scored transcripts in foreground set ", i))
            i <<- i + 1
            return(result)
        })

    background_scores <-
        score_transcripts(
            background_set,
            motifs = motifs,
            max_hits = max_hits,
            threshold_method = threshold_method,
            threshold_value = threshold_value,
            n_cores = n_cores,
            cache = cache
        )
    message("scored transcripts in background set")

    i <- 1
    enrichment_dfs <-
        lapply(foreground_scores, function(scores_per_condition) {
            enrichment_df <-
                calculate_motif_enrichment(
                    scores_per_condition$df,
                    background_scores$df,
                    background_scores$total_sites,
                    background_scores$absolute_hits,
                    length(foreground_sets[[i]]),
                    max_fg_permutations = max_fg_permutations,
                    min_fg_permutations = min_fg_permutations,
                    e = e,
                    p_adjust_method = p_adjust_method
                )
            enrichment_df <-
                dplyr::arrange(enrichment_df, adj_p_value, p_value)
            message(paste0("calculated enrichment for foreground set ", i))
            i <<- i + 1
            return(enrichment_df)
        })
    return(
        list(
            foreground_scores = foreground_scores,
            background_scores = background_scores,
            enrichment_dfs = enrichment_dfs
        )
    )
}

#' @title Matrix-based Spectrum Motif Analysis
#'
#' @description
#' SPMA helps to illuminate the relationship between RBP binding
#' evidence and the transcript
#' sorting criterion, e.g., fold change between treatment and control samples.
#'
#' @param background_set named character vector of ranked sequences
#' (only containing upper case characters A, C, G, T), where the
#' names are RefSeq identifiers
#' and sequence
#' type qualifiers (\code{"3UTR"}, \code{"5UTR"} or \code{"mRNA"}), separated by
#' \code{"|"}, e.g.
#' \code{"NM_010356|3UTR"}. Names are only used to cache results.
#' The sequences in \code{background_set} must be ranked (i.e., sorted).
#' Commonly used sorting criteria are measures of differential expression, such
#' as fold change or signal-to-noise ratio (e.g., between treatment and control
#' samples in gene expression profiling experiments).
#' @inheritParams run_matrix_tsma
#' @inheritParams subdivide_data
#' @inheritParams score_spectrum
#'
#' @return A list with the following components:
#' \tabular{rl}{
#'   \code{foreground_scores} \tab the result of \code{\link{score_transcripts}}
#'   for the foreground
#'   sets (the bins)\cr
#'   \code{background_scores} \tab the result of \code{\link{score_transcripts}}
#'   for the background
#'   set\cr
#'   \code{enrichment_dfs} \tab a list of data frames, returned by
#'   \code{\link{calculate_motif_enrichment}}\cr
#'   \code{spectrum_info_df} \tab a data frame with the SPMA results\cr
#'   \code{spectrum_plots} \tab a list of spectrum plots, as generated by
#'   \code{\link{score_spectrum}}\cr
#'   \code{classifier_scores} \tab a list of classifier scores, as returned by
#'   \code{\link{spectrum_classifier}}
#' }
#'
#' @details
#' In order to investigate how motif targets are distributed across a
#' spectrum of
#' transcripts (e.g., all transcripts of a platform, ordered by fold change),
#' Spectrum Motif Analysis visualizes the gradient of RBP binding evidence
#' across all transcripts.
#'
#' The matrix-based approach skips the \emph{k}-merization step of the
#' \emph{k}-mer-based approach
#' and instead scores the transcript sequence as a whole with a position
#' specific scoring matrix.
#'
#' For each sequence in foreground and background sets and each sequence motif,
#' the scoring algorithm evaluates the score for each sequence position.
#' Positions with
#' a relative score greater than a certain threshold are considered hits, i.e.,
#' putative binding sites.
#'
#' By scoring all sequences in foreground and background sets, a hit count
#' for each motif and
#' each set is obtained, which is used to calculate enrichment values and
#' associated p-values
#' in the same way in which motif-compatible hexamer enrichment values are
#' calculated in
#' the \emph{k}-mer-based approach. P-values are adjusted with one of the
#' available adjustment methods.
#'
#' An advantage of the matrix-based approach is the possibility of detecting
#' clusters of
#' binding sites. This can be done by counting regions with many hits using
#' positional
#' hit information or by simply applying a hit count threshold per
#' sequence, e.g., only
#' sequences with more than some number of hits are considered. Homotypic
#' clusters of RBP
#' binding sites may play a similar role as clusters of transcription factors.
#'
#' @examples
#' # example data set
#' background_df <- transite:::ge$background_df
#' # sort sequences by signal-to-noise ratio
#' background_df <- dplyr::arrange(background_df, value)
#' # character vector of named and ranked (by signal-to-noise ratio) sequences
#' background_set <- gsub("T", "U", background_df$seq)
#' names(background_set) <- paste0(background_df$refseq, "|",
#'   background_df$seq_type)
#'
#' results <- run_matrix_spma(background_set,
#'                          motifs = get_motif_by_id("M178_0.6"),
#'                          n_bins = 20,
#'                          max_fg_permutations = 10000)
#'
#' \dontrun{
#' results <- run_matrix_spma(background_set) }
#'
#' @family SPMA functions
#' @family matrix functions
#' @importFrom stats p.adjust
#' @importFrom dplyr filter
#' @export
run_matrix_spma <- function(background_set,
                            motifs = NULL,
                            n_bins = 40,
                            max_model_degree = 1,
                            max_cs_permutations = 10000000,
                            min_cs_permutations = 5000,
                            max_hits = 5,
                            threshold_method = "p_value",
                            threshold_value = 0.25^6,
                            max_fg_permutations = 1000000,
                            min_fg_permutations = 1000,
                            e = 5,
                            p_adjust_method = "BH",
                            n_cores = 1,
                            cache = paste0(tempdir(), "/sc/")) {
    # avoid CRAN note
    motif_id <- motif_rbps <- adj_r_squared <- degree <-
        residuals <- slope <- NULL
    f_statistic <- f_statistic_p_value <- f_statistic_adj_p_value <- NULL
    consistency_score <-
        consistency_score_p_value <-
        consistency_score_adj_p_value <-
        consistency_score_n <- NULL
    n_significant <-
        n_very_significant <-
        n_extremely_significant <-
        aggregate_classifier_score <- NULL

    foreground_sets <- subdivide_data(background_set, n_bins)

    results <- run_matrix_tsma(
        foreground_sets,
        background_set,
        motifs = motifs,
        max_hits = max_hits,
        threshold_method = threshold_method,
        threshold_value = threshold_value,
        max_fg_permutations = max_fg_permutations,
        min_fg_permutations = min_fg_permutations,
        e = e,
        p_adjust_method = p_adjust_method,
        n_cores = n_cores,
        cache = cache
    )

    if (length(results$enrichment_dfs) > 0) {
        enrichment_df <- do.call("rbind", results$enrichment_dfs)
        enrichment_df$adj_p_value <-
            stats::p.adjust(enrichment_df$p_value, method = p_adjust_method)

        if (is.null(motifs)) {
            motifs <- get_motifs()
        }
        spectrum_info <- lapply(motifs, function(motif) {
            motif_data_df <- dplyr::filter(enrichment_df,
                                           motif_id == motifId(motif))
            values <- motif_data_df$enrichment
            values[values == 0] <-
                0.01 # avoid -Inf after taking the log
            score <-
                score_spectrum(
                    log(values),
                    motif_data_df$adj_p_value,
                    max_model_degree = max_model_degree,
                    max_cs_permutations = max_cs_permutations,
                    min_cs_permutations = min_cs_permutations
                )

            n_significant <-
                sum(motif_data_df$adj_p_value <= 0.05)
            classifier_score <-
                spectrum_classifier(
                    spectrumAdjRSquared(score),
                    spectrumDegree(score),
                    spectrumSlope(score),
                    spectrumConsistencyScoreN(score),
                    n_significant,
                    n_bins
                )

            return(
                list(
                    info = list(
                        motif_id = motifId(motif),
                        motif_rbps = paste(motifRbps(motif), collapse = ", "),
                        adj_r_squared = spectrumAdjRSquared(score),
                        degree = spectrumDegree(score),
                        residuals = spectrumResiduals(score),
                        slope = spectrumSlope(score),
                        f_statistic = spectrumFStatistic(score),
                        f_statistic_p_value = spectrumFStatisticPValue(score),
                        consistency_score = spectrumConsistencyScore(score),
                        consistency_score_p_value = spectrumConsistencyScorePValue(score),
                        consistency_score_n = spectrumConsistencyScoreN(score),
                        n_significant = n_significant,
                        n_very_significant = sum(motif_data_df$adj_p_value <= 0.01),
                        n_extremely_significant = sum(motif_data_df$adj_p_value <= 0.001),
                        aggregate_classifier_score = sum(classifier_score)
                    ),
                    spectrum_plot = score@plot,
                    classifier_score = classifier_score
                )
            )
        })
        spectrum_info_df <-
            as.data.frame(do.call("rbind", lapply(spectrum_info, function(x)
                x$info)),
                stringsAsFactors = FALSE
            )
        spectrum_info_df$motif_id <-
            as.character(spectrum_info_df$motif_id)
        spectrum_info_df$motif_rbps <-
            as.character(spectrum_info_df$motif_rbps)
        spectrum_info_df$adj_r_squared <-
            as.numeric(spectrum_info_df$adj_r_squared)
        spectrum_info_df$degree <-
            as.integer(spectrum_info_df$degree)
        spectrum_info_df$residuals <-
            as.numeric(spectrum_info_df$residuals)
        spectrum_info_df$slope <-
            as.numeric(spectrum_info_df$slope)
        spectrum_info_df$f_statistic <-
            as.numeric(spectrum_info_df$f_statistic)
        spectrum_info_df$f_statistic_p_value <-
            as.numeric(spectrum_info_df$f_statistic_p_value)
        spectrum_info_df$consistency_score <-
            as.numeric(spectrum_info_df$consistency_score)
        spectrum_info_df$consistency_score_p_value <-
            as.numeric(spectrum_info_df$consistency_score_p_value)
        spectrum_info_df$consistency_score_n <-
            as.integer(spectrum_info_df$consistency_score_n)
        spectrum_info_df$n_significant <-
            as.integer(spectrum_info_df$n_significant)
        spectrum_info_df$n_very_significant <-
            as.integer(spectrum_info_df$n_very_significant)
        spectrum_info_df$n_extremely_significant <-
            as.integer(spectrum_info_df$n_extremely_significant)
        spectrum_info_df$aggregate_classifier_score <-
            as.integer(spectrum_info_df$aggregate_classifier_score)

        spectrum_info_df$f_statistic_adj_p_value <-
            stats::p.adjust(spectrum_info_df$f_statistic_p_value,
                            method = p_adjust_method)
        spectrum_info_df$consistency_score_adj_p_value <-
            stats::p.adjust(spectrum_info_df$consistency_score_p_value,
                            method = p_adjust_method)

        spectrum_info_df <-
            dplyr::select(
                spectrum_info_df,
                motif_id,
                motif_rbps,
                adj_r_squared,
                degree,
                residuals,
                slope,
                f_statistic,
                f_statistic_p_value,
                f_statistic_adj_p_value,
                consistency_score,
                consistency_score_p_value,
                consistency_score_adj_p_value,
                consistency_score_n,
                n_significant,
                n_very_significant,
                n_extremely_significant,
                aggregate_classifier_score
            )

        spectrum_plots <-
            lapply(spectrum_info, function(x)
                x$spectrum_plot)
        classifier_scores <-
            lapply(spectrum_info, function(x)
                x$classifier_score)
    } else {
        spectrum_info_df <- data.frame(
            motif_id = character(0),
            motif_rbps = character(0),
            adj_r_squared = numeric(0),
            degree = integer(0),
            residuals = numeric(0),
            slope = numeric(0),
            f_statistic = numeric(0),
            f_statistic_p_value = numeric(0),
            f_statistic_adj_p_value = numeric(0),
            consistency_score = numeric(0),
            consistency_score_p_value = numeric(0),
            consistency_score_adj_p_value = numeric(0),
            consistency_score_n = integer(0),
            n_significant = integer(0),
            n_very_significant = integer(0),
            n_extremely_significant = integer(0),
            aggregate_classifier_score = integer(0)
        )
        spectrum_plots <- NULL
        warning("no sequences in any condition")
    }

    return(
        list(
            foreground_scores = results$foreground_scores,
            background_scores = results$background_scores,
            enrichment_dfs = results$enrichment_dfs,
            spectrum_info_df = spectrum_info_df,
            spectrum_plots = spectrum_plots,
            classifier_scores = classifier_scores
        )
    )
}

#' @title \emph{k}-mer-based Transcript Set Motif Analysis
#'
#' @description
#' Calculates the enrichment of putative binding sites in foreground sets
#' versus a background set
#' using \emph{k}-mers to identify putative binding sites
#'
#' @param motifs a list of motifs that is used to score the specified sequences.
#' If \code{is.null(motifs)} then all Transite motifs are used.
#' @param fg_permutations numer of foreground permutations
#' @param kmer_significance_threshold p-value threshold for significance,
#' e.g., \code{0.05} or
#' \code{0.01} (used for volcano plots)
#' @param produce_plot if \code{TRUE} volcano plots and distribution plots
#' are created
#' @param p_combining_method one of the following: Fisher (1932)
#' (\code{"fisher"}), Stouffer (1949),
#' Liptak (1958) (\code{"SL"}), Mudholkar and George (1979)
#' (\code{"MG"}), and Tippett (1931)
#' (\code{"tippett"}) (see \code{\link{p_combine}})
#' @inheritParams calculate_kmer_enrichment
#'
#' @return A list of lists (one for each transcript set) with the
#' following components:
#' \tabular{rl}{
#'   \code{enrichment_df} \tab the result of
#'   \code{\link{compute_kmer_enrichment}} \cr
#'   \code{motif_df} \tab \cr
#'   \code{motif_kmers_dfs} \tab \cr
#'   \code{volcano_plots} \tab volcano plots for each
#'   motif (see \code{\link{draw_volcano_plot}}) \cr
#'   \code{perm_test_plots} \tab plots of the empirical distribution of
#'   \emph{k}-mer enrichment values for each motif \cr
#'   \code{enriched_kmers_combined_p_values} \tab \cr
#'   \code{depleted_kmers_combined_p_values} \tab
#' }
#'
#' @details
#' Motif transcript set analysis can be used to identify RNA binding
#' proteins, whose targets are
#' significantly overrepresented or underrepresented in certain sets
#' of transcripts.
#'
#' The aim of Transcript Set Motif Analysis (TSMA) is to identify the
#' overrepresentation
#' and underrepresentation of potential RBP targets (binding sites)
#' in a set (or sets) of
#' sequences, i.e., the foreground set, relative to the entire population
#' of sequences.
#' The latter is called background set, which can be composed of all
#' sequences of the genes
#' of a microarray platform or all sequences of an organism or any
#' other meaningful
#' superset of the foreground sets.
#'
#' The \emph{k}-mer-based approach breaks the sequences of foreground
#' and background sets into
#' \emph{k}-mers and calculates the enrichment on a \emph{k}-mer level.
#' In this case, motifs are
#' not represented as position weight matrices, but as lists of \emph{k}-mers.
#'
#' Statistically significantly enriched or depleted \emph{k}-mers
#' are then used to
#' calculate a score for each RNA-binding protein, which quantifies its
#' target overrepresentation.
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
#' background_set <- unique(c(foreground_set1, foreground_set2, c(
#'   "CCACACAC", "CUCAUUGGAG", "ACUUUGGGACA", "CAGGUCAGCA",
#'   "CCACACCGG", "GUCAUCAGU", "GUCAGUCC", "CAGGUCAGGGGCA"
#' )))
#'
#' # run k-mer based TSMA with all Transite motifs (recommended):
#' # results <- run_kmer_tsma(foreground_sets, background_set)
#'
#' # run TSMA with one motif:
#' motif_db <- get_motif_by_id("M178_0.6")
#' results <- run_kmer_tsma(foreground_sets, background_set, motifs = motif_db)
#' \dontrun{
#' # define example sequence sets for foreground and background
#' foreground_set1 <- gsub("T", "U", transite:::ge$foreground1_df$seq)
#' foreground_set2 <- gsub("T", "U", transite:::ge$foreground2_df$seq)
#' foreground_sets <- list(foreground_set1, foreground_set2)
#' background_set <- gsub("T", "U", transite:::ge$background_df$seq)
#'
#' # run TSMA with all Transite motifs
#' results <- run_kmer_tsma(foreground_sets, background_set)
#'
#' # run TSMA with a subset of Transite motifs
#' results <- run_kmer_tsma(foreground_sets, background_set,
#'   motifs = get_motif_by_rbp("ELAVL1"))
#'
#' # run TSMA with user-defined motif
#' toy_motif <- create_kmer_motif(
#'   "toy_motif", "example RBP",
#'   c("AACCGG", "AAAACG", "AACACG"), "example type", "example species", "user"
#' )
#' results <- run_matrix_tsma(foreground_sets, background_set,
#'   motifs = list(toy_motif))
#' }
#'
#' @family TSMA functions
#' @family \emph{k}-mer functions
#' @importFrom stats p.adjust
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @export
run_kmer_tsma <- function(foreground_sets,
                          background_set,
                          motifs = NULL,
                          k = 6,
                          fg_permutations = 5000,
                          kmer_significance_threshold = 0.01,
                          produce_plot = TRUE,
                          p_adjust_method = "BH",
                          p_combining_method = "fisher",
                          n_cores = 1) {
    # avoid CRAN note
    kmer <- enrichment <- p_value <- adj_p_value <- NULL

    if (is.null(motifs)) {
        motifs <- get_motifs()
    }

    enrichment_result <-
        calculate_kmer_enrichment(
            foreground_sets,
            background_set,
            k,
            p_adjust_method = p_adjust_method,
            n_cores = n_cores
        )
    message("calculated enrichment for all foreground sets")

    i <- 1
    if (!is.null(enrichment_result)) {
        set_sizes <- unique(unlist(lapply(foreground_sets, length)))
        random_enrichments <- new.env()
        for (set_size in set_sizes) {
            if (set_size < 10) {
                adapted_fg_permutations <- min(fg_permutations, 100)
            } else if (set_size < 30) {
                adapted_fg_permutations <- min(fg_permutations, 500)
            } else {
                adapted_fg_permutations <- fg_permutations
            }
            assign(
                as.character(set_size),
                generate_permuted_enrichments(
                    set_size,
                    background_set,
                    k,
                    n_permutations = adapted_fg_permutations,
                    n_cores = n_cores
                ),
                envir = random_enrichments
            )
        }
        message("calculated enrichment for all permuted sets")

        motif_result <-
            lapply(seq_len(length(foreground_sets)), function(i) {
                kmers_df <- enrichment_result$dfs[[i]]
                kmers_df$kmer <- enrichment_result$kmers

                foreground_result <-
                    lapply(motifs, function(motif) {
                        if (k == 6) {
                            rbp_kmers <- motifHexamers(motif)
                        } else if (k == 7) {
                            rbp_kmers <- motifHeptamers(motif)
                        }

                        idx <-
                            which(enrichment_result$kmers %in% rbp_kmers)
                        motif_kmers_df <- kmers_df[idx, ]
                        motif_kmers_df <-
                            dplyr::select(
                                motif_kmers_df,
                                kmer,
                                enrichment,
                                p_value,
                                adj_p_value
                            )

                        geo_mean <-
                            geometric_mean(motif_kmers_df$enrichment)
                        perm_test <-
                            perm_test_geometric_mean(
                                geo_mean,
                                rbp_kmers,
                                get(
                                    as.character(length(foreground_sets[[i]])),
                                    envir = random_enrichments,
                                    inherits = FALSE
                                ),
                                alternative = "two_sided",
                                conf_level = 0.95,
                                produce_plot = produce_plot
                            )

                        motif_kmers_enriched_df <-
                            dplyr::filter(motif_kmers_df, enrichment > 1)
                        enriched_kmers_combined_p_value <-
                            p_combine(motif_kmers_enriched_df$adj_p_value,
                                      method = p_combining_method
                            )

                        motif_kmers_depleted_df <-
                            dplyr::filter(motif_kmers_df, enrichment < 1)
                        depleted_kmers_combined_p_value <-
                            p_combine(motif_kmers_depleted_df$adj_p_value,
                                      method = p_combining_method
                            )

                        if (produce_plot) {
                            volcano_plot <-
                                draw_volcano_plot(
                                    kmers_df,
                                    rbp_kmers,
                                    paste(motifRbps(motif), collapse = ", "),
                                    kmer_significance_threshold
                                )
                        } else {
                            volcano_plot <- NULL
                        }

                        return(
                            list(
                                df = list(
                                    motif_id = motifId(motif),
                                    motif_rbps = paste0(motifRbps(motif),
                                                        collapse = ", "),
                                    geo_mean_enrichment = geo_mean,
                                    p_value_estimate = perm_test$p_value_estimate,
                                    cont_int_low = perm_test$conf_int[1],
                                    cont_int_high = perm_test$conf_int[2],
                                    enriched_kmers_combined_p_value = enriched_kmers_combined_p_value$p_value,
                                    depleted_kmers_combined_p_value = depleted_kmers_combined_p_value$p_value
                                ),
                                motif_kmers_df = motif_kmers_df,
                                volcano_plot = volcano_plot,
                                perm_test_plot = perm_test$plot,
                                enriched_kmers_combined_p_value = enriched_kmers_combined_p_value,
                                depleted_kmers_combined_p_value = depleted_kmers_combined_p_value
                            )
                        )
                    })

                df <-
                    as.data.frame(do.call(
                        "rbind",
                        lapply(foreground_result, function(x)
                            x$df)
                    ), stringsAsFactors = FALSE)
                df$motif_id <- as.character(df$motif_id)
                df$motif_rbps <- as.character(df$motif_rbps)
                df$geo_mean_enrichment <-
                    as.numeric(df$geo_mean_enrichment)
                df$p_value_estimate <-
                    as.numeric(df$p_value_estimate)
                df$cont_int_low <- as.numeric(df$cont_int_low)
                df$cont_int_high <- as.numeric(df$cont_int_high)
                df$enriched_kmers_combined_p_value <-
                    as.numeric(df$enriched_kmers_combined_p_value)
                df$depleted_kmers_combined_p_value <-
                    as.numeric(df$depleted_kmers_combined_p_value)

                df$adj_p_value_estimate <-
                    stats::p.adjust(df$p_value_estimate,
                                    method = p_adjust_method)
                df$adj_enriched_kmers_combined_p_value <-
                    stats::p.adjust(df$enriched_kmers_combined_p_value,
                                    method = p_adjust_method)
                df$adj_depleted_kmers_combined_p_value <-
                    stats::p.adjust(df$depleted_kmers_combined_p_value,
                                    method = p_adjust_method)

                motif_kmers_dfs <-
                    lapply(foreground_result, function(x)
                        x$motif_kmers_df)
                volcano_plots <-
                    lapply(foreground_result, function(x)
                        x$volcano_plot)
                perm_test_plots <-
                    lapply(foreground_result, function(x)
                        x$perm_test_plot)
                enriched_kmers_combined_p_values <-
                    lapply(foreground_result, function(x)
                        x$enriched_kmers_combined_p_value)
                depleted_kmers_combined_p_values <-
                    lapply(foreground_result, function(x)
                        x$depleted_kmers_combined_p_value)

                return(
                    list(
                        df = df,
                        motif_kmers_dfs = motif_kmers_dfs,
                        volcano_plots = volcano_plots,
                        perm_test_plots = perm_test_plots,
                        enriched_kmers_combined_p_values = enriched_kmers_combined_p_values,
                        depleted_kmers_combined_p_values = depleted_kmers_combined_p_values
                    )
                )
            })

        result <-
            lapply(seq_len(length(foreground_sets)), function(i) {
                enrichment_df <- enrichment_result$dfs[[i]]
                enrichment_df$kmer <- enrichment_result$kmers
                return(
                    list(
                        enrichment_df = enrichment_df,
                        motif_df = motif_result[[i]]$df,
                        motif_kmers_dfs = motif_result[[i]]$motif_kmers_dfs,
                        volcano_plots = motif_result[[i]]$volcano_plots,
                        perm_test_plots = motif_result[[i]]$perm_test_plots,
                        enriched_kmers_combined_p_values = motif_result[[i]]$enriched_kmers_combined_p_values,
                        depleted_kmers_combined_p_values = motif_result[[i]]$depleted_kmers_combined_p_values
                    )
                )
            })
        return(result)
    } else {
        return(NULL)
    }
}



#' @title \emph{k}-mer-based Spectrum Motif Analysis
#'
#' @description
#' SPMA helps to illuminate the relationship between RBP binding evidence
#' and the transcript
#' sorting criterion, e.g., fold change between treatment and control samples.
#'
#' @param background_set character vector of ranked sequences, either DNA
#' (only containing upper case characters A, C, G, T) or RNA (A, C, G, U).
#' The sequences in \code{background.set} must be ranked (i.e., sorted).
#' Commonly used sorting criteria are measures of differential expression, such
#' as fold change or signal-to-noise ratio (e.g., between treatment and control
#' samples in gene expression profiling experiments).
#'
#' @inheritParams run_kmer_tsma
#' @inheritParams subdivide_data
#' @inheritParams score_spectrum
#'
#' @return A list with the following components:
#' \tabular{rl}{
#'   \code{foreground_scores} \tab the result of \code{\link{run_kmer_tsma}}
#'   for the binned data\cr
#'   \code{spectrum_info_df} \tab a data frame with the SPMA results\cr
#'   \code{spectrum_plots} \tab a list of spectrum plots, as generated by
#'   \code{\link{score_spectrum}}\cr
#'   \code{classifier_scores} \tab a list of classifier scores, as returned by
#'   \code{\link{spectrum_classifier}}
#' }
#'
#' @details
#' In order to investigate how motif targets are distributed across a
#' spectrum of
#' transcripts (e.g., all transcripts of a platform, ordered by fold change),
#' Spectrum Motif Analysis visualizes the gradient of RBP binding evidence
#' across all transcripts.
#'
#' The \emph{k}-mer-based approach differs from the matrix-based approach by
#' how the sequences are
#' scored. Here, sequences are broken into \emph{k}-mers, i.e.,
#' oligonucleotide sequences of
#' \emph{k} bases.
#' And only statistically significantly enriched or depleted \emph{k}-mers
#' are then used to
#' calculate a score for each RNA-binding protein, which quantifies its
#' target overrepresentation.
#'
#' @examples
#' # example data set
#' background_df <- transite:::ge$background_df
#' # sort sequences by signal-to-noise ratio
#' background_df <- dplyr::arrange(background_df, value)
#' # character vector of named and ranked (by signal-to-noise ratio) sequences
#' background_set <- gsub("T", "U", background_df$seq)
#' names(background_set) <- paste0(background_df$refseq, "|",
#'   background_df$seq_type)
#'
#' results <- run_kmer_spma(background_set,
#'                        motifs = get_motif_by_id("M178_0.6"),
#'                        n_bins = 20,
#'                        fg_permutations = 10)
#'
#' \dontrun{
#' results <- run_kmer_spma(background_set)}
#'
#' @family SPMA functions
#' @family \emph{k}-mer functions
#' @importFrom stats p.adjust
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @export
run_kmer_spma <- function(background_set,
                          motifs = NULL,
                          k = 6,
                          n_bins = 40,
                          max_model_degree = 1,
                          max_cs_permutations = 10000000,
                          min_cs_permutations = 5000,
                          fg_permutations = 5000,
                          p_adjust_method = "BH",
                          p_combining_method = "fisher",
                          n_cores = 1) {
    # avoid CRAN note
    motif_id <-
        motif_rbps <-
        adj_r_squared <- degree <- residuals <- slope <- NULL
    f_statistic <-
        f_statistic_p_value <- f_statistic_adj_p_value <- NULL
    consistency_score <-
        consistency_score_p_value <-
        consistency_score_adj_p_value <- NULL
    consistency_score_n <- aggregate_classifier_score <- NULL
    n_significant <-
        n_very_significant <- n_extremely_significant <- NULL

    foreground_sets <- subdivide_data(background_set, n_bins)

    results <- run_kmer_tsma(
        foreground_sets,
        background_set,
        motifs = motifs,
        k = k,
        fg_permutations = fg_permutations,
        kmer_significance_threshold = 0.01,
        produce_plot = FALSE,
        p_adjust_method = p_adjust_method,
        p_combining_method = p_combining_method,
        n_cores = n_cores
    )

    if (length(results) > 0) {
        dfs <- lapply(results, function(result) {
            result$motif_df
        })
        enrichment_df <- do.call("rbind", dfs)
        enrichment_df$adj_p_value_estimate <-
            stats::p.adjust(enrichment_df$p_value_estimate,
                            method = p_adjust_method
            )

        if (is.null(motifs)) {
            motifs <- get_motifs()
        }
        spectrum_info <- lapply(motifs, function(motif) {
            motif_data_df <- dplyr::filter(enrichment_df,
                                           motif_id == motifId(motif))
            values <- motif_data_df$geo_mean_enrichment
            values[values == 0] <-
                0.01 # avoid -Inf after taking the log
            score <-
                score_spectrum(
                    log(values),
                    motif_data_df$adj_p_value_estimate,
                    max_model_degree = max_model_degree,
                    max_cs_permutations = max_cs_permutations,
                    min_cs_permutations = min_cs_permutations
                )

            n_significant <-
                sum(motif_data_df$adj_p_value_estimate <= 0.05)
            classifier_score <-
                spectrum_classifier(
                    spectrumAdjRSquared(score),
                    spectrumDegree(score),
                    spectrumSlope(score),
                    spectrumConsistencyScoreN(score),
                    n_significant,
                    n_bins
                )

            return(
                list(
                    info = list(
                        motif_id = motifId(motif),
                        motif_rbps = paste(motifRbps(motif), collapse = ", "),
                        adj_r_squared = spectrumAdjRSquared(score),
                        degree = spectrumDegree(score),
                        residuals = spectrumResiduals(score),
                        slope = spectrumSlope(score),
                        f_statistic = spectrumFStatistic(score),
                        f_statistic_p_value = spectrumFStatisticPValue(score),
                        consistency_score = spectrumConsistencyScore(score),
                        consistency_score_p_value = spectrumConsistencyScorePValue(score),
                        consistency_score_n = spectrumConsistencyScoreN(score),
                        n_significant = n_significant,
                        n_very_significant = sum(motif_data_df$adj_p_value <= 0.01),
                        n_extremely_significant = sum(motif_data_df$adj_p_value <= 0.001),
                        aggregate_classifier_score = sum(classifier_score)
                    ),
                    spectrum_plot = score@plot,
                    classifier_score = classifier_score
                )
            )
        })
        spectrum_info_df <-
            as.data.frame(do.call("rbind", lapply(spectrum_info, function(x)
                x$info)),
                stringsAsFactors = FALSE
            )
        spectrum_info_df$motif_id <-
            as.character(spectrum_info_df$motif_id)
        spectrum_info_df$motif_rbps <-
            as.character(spectrum_info_df$motif_rbps)
        spectrum_info_df$adj_r_squared <-
            as.numeric(spectrum_info_df$adj_r_squared)
        spectrum_info_df$degree <-
            as.integer(spectrum_info_df$degree)
        spectrum_info_df$residuals <-
            as.numeric(spectrum_info_df$residuals)
        spectrum_info_df$slope <-
            as.numeric(spectrum_info_df$slope)
        spectrum_info_df$f_statistic <-
            as.numeric(spectrum_info_df$f_statistic)
        spectrum_info_df$f_statistic_p_value <-
            as.numeric(spectrum_info_df$f_statistic_p_value)
        spectrum_info_df$consistency_score <-
            as.numeric(spectrum_info_df$consistency_score)
        spectrum_info_df$consistency_score_p_value <-
            as.numeric(spectrum_info_df$consistency_score_p_value)
        spectrum_info_df$consistency_score_n <-
            as.integer(spectrum_info_df$consistency_score_n)
        spectrum_info_df$n_significant <-
            as.integer(spectrum_info_df$n_significant)
        spectrum_info_df$n_very_significant <-
            as.integer(spectrum_info_df$n_very_significant)
        spectrum_info_df$n_extremely_significant <-
            as.integer(spectrum_info_df$n_extremely_significant)
        spectrum_info_df$aggregate_classifier_score <-
            as.integer(spectrum_info_df$aggregate_classifier_score)

        spectrum_info_df$f_statistic_adj_p_value <-
            stats::p.adjust(spectrum_info_df$f_statistic_p_value,
                            method = p_adjust_method)
        spectrum_info_df$consistency_score_adj_p_value <-
            stats::p.adjust(spectrum_info_df$consistency_score_p_value,
                            method = p_adjust_method)

        spectrum_info_df <-
            dplyr::select(
                spectrum_info_df,
                motif_id,
                motif_rbps,
                adj_r_squared,
                degree,
                residuals,
                slope,
                f_statistic,
                f_statistic_p_value,
                f_statistic_adj_p_value,
                consistency_score,
                consistency_score_p_value,
                consistency_score_adj_p_value,
                consistency_score_n,
                n_significant,
                n_very_significant,
                n_extremely_significant,
                aggregate_classifier_score
            )

        spectrum_plots <-
            lapply(spectrum_info, function(x)
                x$spectrum_plot)
        classifier_scores <-
            lapply(spectrum_info, function(x)
                x$classifier_score)
    } else {
        spectrum_info_df <- data.frame(
            motif_id = character(0),
            motif_rbps = character(0),
            adj_r_squared = numeric(0),
            degree = integer(0),
            residuals = numeric(0),
            slope = numeric(0),
            f_statistic = numeric(0),
            f_statistic_p_value = numeric(0),
            f_statistic_adj_p_value = numeric(0),
            consistency_score = numeric(0),
            consistency_score_p_value = numeric(0),
            consistency_score_adj_p_value = numeric(0),
            consistency_score_n = integer(0),
            n_significant = integer(0),
            n_very_significant = integer(0),
            n_extremely_significant = integer(0),
            aggregate_classifier_score = integer(0)
        )
        spectrum_plots <- NULL
        warning("no sequences in any condition")
    }

    return(
        list(
            foreground_scores = results,
            spectrum_info_df = spectrum_info_df,
            spectrum_plots = spectrum_plots,
            classifier_scores = classifier_scores
        )
    )
}
