# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @title Score Sequences with PWM
#'
#' @description
#' C++ implementation of PWM scoring algorithm
#'
#' @param sequences list of sequences
#' @param pwm position weight matrix
#'
#' @return list of PWM scores for each sequence
#' @examples
#' motif <- get_motif_by_id("M178_0.6")[[1]]
#' sequences <- c("CAACAGCCUUAAUU", "CAGUCAAGACUCC", "CUUUGGGGAAU",
#'                "UCAUUUUAUUAAA", "AAUUGGUGUCUGGAUACUUCCCUGUACAU",
#'                "AUCAAAUUA", "UGUGGGG", "GACACUUAAAGAUCCU",
#'                "UAGCAUUAACUUAAUG", "AUGGA", "GAAGAGUGCUCA", "AUAGAC",
#'                "AGUUC", "CCAGUAA")
#' seq_char_vectors <- lapply(sequences, function(seq) {
#'   unlist(strsplit(seq, ""))
#' })
#' score_sequences(seq_char_vectors, as.matrix(get_motif_matrix(motif)))
#'
#' @export
score_sequences <- function(sequences, pwm) {
    .Call('_transite_score_sequences', PACKAGE = 'transite', sequences, pwm)
}

#' @title Local Consistency Score
#'
#' @description
#' C++ implementation of Local Consistency Score algorithm.
#'
#' @param x numeric vector that contains values for shuffling
#' @param numPermutations maximum number of permutations performed in
#' Monte Carlo test
#' for consistency score
#' @param minPermutations minimum number of permutations performed in
#' Monte Carlo test
#' for consistency score
#' @param e stop criterion for consistency score Monte Carlo test:
#' aborting permutation
#' process after observing \code{e} random consistency values with
#' more extreme values
#' than the actual consistency value
#' @return list with \code{score}, \code{p_value}, and \code{n} components,
#' where \code{score} is the raw local consistency score (usually not used),
#' \code{p_value} is the associated p-value for that score, obtained by
#' Monte Carlo testing, and \code{n} is the number of permutations performed
#' in the Monte Carlo test (the higher, the more significant)
#'
#' @examples
#' poor_enrichment_spectrum <- c(0.1, 0.5, 0.6, 0.4,
#'   0.7, 0.6, 1.2, 1.1, 1.8, 1.6)
#' local_consistency <- calculate_local_consistency(poor_enrichment_spectrum,
#'   1000000, 1000, 5)
#'
#' enrichment_spectrum <- c(0.1, 0.3, 0.6, 0.7, 0.8,
#'   0.9, 1.2, 1.4, 1.6, 1.4)
#' local_consistency <- calculate_local_consistency(enrichment_spectrum,
#'   1000000, 1000, 5)
#' @export
calculate_local_consistency <- function(x, numPermutations, minPermutations, e) {
    .Call('_transite_calculate_local_consistency', PACKAGE = 'transite', x, numPermutations, minPermutations, e)
}

#' @title Motif Enrichment calculation
#'
#' @description
#' C++ implementation of Motif Enrichment calculation
#'
#' @param absoluteHits number of putative binding sites per sequence
#' (returned by \code{\link{score_transcripts}})
#' @param totalSites number of potential binding sites per sequence
#' (returned by \code{\link{score_transcripts}})
#' @param relHitsForeground relative number of hits in foreground set
#' @param n number of sequences in the foreground set
#' @param maxPermutations maximum number of foreground permutations
#' performed in
#' Monte Carlo test for enrichment score
#' @param minPermutations minimum number of foreground permutations
#' performed in
#' Monte Carlo test for enrichment score
#' @param e stop criterion for enrichment score Monte Carlo test:
#' aborting permutation process
#' after observing \code{e} random enrichment values with more extreme
#' values than the actual
#' enrichment value
#'
#' @return list with p-value and number of iterations of Monte Carlo sampling
#' for foreground enrichment
#'
#' @examples
#' foreground_seqs <- c("CAGUCAAGACUCC", "AAUUGGUUGUGGGGCUUCCCUGUACAU",
#'                      "AGAU", "CCAGUAA", "UGUGGGG")
#' background_seqs <- c(foreground_seqs, "CAACAGCCUUAAUU", "CUUUGGGGAAU",
#'                      "UCAUUUUAUUAAA", "AUCAAAUUA", "GACACUUAAAGAUCCU",
#'                      "UAGCAUUAACUUAAUG", "AUGGA", "GAAGAGUGCUCA",
#'                      "AUAGAC", "AGUUC")
#' motif_db <- get_motif_by_id("M178_0.6")
#' fg <- score_transcripts(foreground_seqs, cache = FALSE,
#'   motifs = motif_db)
#' bg <- score_transcripts(background_seqs, cache = FALSE,
#'   motifs = motif_db)
#'
#' mc_result <- calculate_transcript_mc(unlist(bg$absolute_hits),
#'  unlist(bg$total_sites),
#'  fg$df$absolute_hits / fg$df$total_sites,
#'  length(foreground_seqs), 1000, 500, 5)
#' @export
calculate_transcript_mc <- function(absoluteHits, totalSites, relHitsForeground, n, maxPermutations, minPermutations, e) {
    .Call('_transite_calculate_transcript_mc', PACKAGE = 'transite', absoluteHits, totalSites, relHitsForeground, n, maxPermutations, minPermutations, e)
}

