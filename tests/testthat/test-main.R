context("main functions")

test_that("run_matrix_tsma", {
    foreground_set1 <- c(
        "CAACAGCCUUAAUU", "CAGUCAAGACUCC", "CUUUGGGGAAU",
        "UCAUUUUAUUAAA", "AAUUGGUGUCUGGAUACUUCCCUGUACAU",
        "AUCAAAUUA", "AGAU", "GACACUUAAAGAUCCU",
        "UAGCAUUAACUUAAUG", "AUGGA", "GAAGAGUGCUCA",
        "AUAGAC", "AGUUC", "CCAGUAA"
    )
    names(foreground_set1) <- c(
        "NM_1_DUMMY|3UTR", "NM_2_DUMMY|3UTR", "NM_3_DUMMY|3UTR",
        "NM_4_DUMMY|3UTR", "NM_5_DUMMY|3UTR", "NM_6_DUMMY|3UTR",
        "NM_7_DUMMY|3UTR",
        "NM_8_DUMMY|3UTR", "NM_9_DUMMY|3UTR", "NM_10_DUMMY|3UTR",
        "NM_11_DUMMY|3UTR",
        "NM_12_DUMMY|3UTR", "NM_13_DUMMY|3UTR", "NM_14_DUMMY|3UTR"
    )

    foreground_set2 <- c("UUAUUUA", "AUCCUUUACA", "UUUUUUU", "UUUCAUCAUU")
    names(foreground_set2) <- c(
        "NM_15_DUMMY|3UTR", "NM_16_DUMMY|3UTR", "NM_17_DUMMY|3UTR",
        "NM_18_DUMMY|3UTR"
    )

    foreground_sets <- list(foreground_set1, foreground_set2)

    background_set <- c(
        "CAACAGCCUUAAUU", "CAGUCAAGACUCC", "CUUUGGGGAAU",
        "UCAUUUUAUUAAA", "AAUUGGUGUCUGGAUACUUCCCUGUACAU",
        "AUCAAAUUA", "AGAU", "GACACUUAAAGAUCCU",
        "UAGCAUUAACUUAAUG", "AUGGA", "GAAGAGUGCUCA",
        "AUAGAC", "AGUUC", "CCAGUAA",
        "UUAUUUA", "AUCCUUUACA", "UUUUUUU", "UUUCAUCAUU",
        "CCACACAC", "CUCAUUGGAG", "ACUUUGGGACA", "CAGGUCAGCA"
    )
    names(background_set) <- c(
        "NM_1_DUMMY|3UTR", "NM_2_DUMMY|3UTR", "NM_3_DUMMY|3UTR",
        "NM_4_DUMMY|3UTR", "NM_5_DUMMY|3UTR", "NM_6_DUMMY|3UTR",
        "NM_7_DUMMY|3UTR",
        "NM_8_DUMMY|3UTR", "NM_9_DUMMY|3UTR", "NM_10_DUMMY|3UTR",
        "NM_11_DUMMY|3UTR",
        "NM_12_DUMMY|3UTR", "NM_13_DUMMY|3UTR", "NM_14_DUMMY|3UTR",
        "NM_15_DUMMY|3UTR",
        "NM_16_DUMMY|3UTR", "NM_17_DUMMY|3UTR", "NM_18_DUMMY|3UTR",
        "NM_19_DUMMY|3UTR",
        "NM_20_DUMMY|3UTR", "NM_21_DUMMY|3UTR", "NM_22_DUMMY|3UTR"
    )

    motif_db <- get_motif_by_id("M178_0.6")
    results <- run_matrix_tsma(foreground_sets, background_set,
                               motifs = motif_db, cache = FALSE)

    expect_equal(length(results$foreground_scores), 2)
    expect_equal(length(results$foreground_scores[[1]]$total_sites[[1]]),
                 length(foreground_set1))
    expect_equal(length(results$foreground_scores[[2]]$total_sites[[1]]),
                 length(foreground_set2))
    expect_equal(as.numeric(results$foreground_scores[[2]]$total_sites[[1]]),
                 c(1, 4, 1, 4))
    expect_equal(as.numeric(results$foreground_scores[[2]]$df$total_sites),
                 sum(results$foreground_scores[[2]]$total_sites[[1]]))
    expect_equal(length(results$background_scores$total_sites[[1]]),
                 length(background_set))
    expect_equal(length(results$foreground_scores[[1]]$absolute_hits[[1]]),
                 length(foreground_set1))
    expect_equal(length(results$foreground_scores[[2]]$absolute_hits[[1]]),
                 length(foreground_set2))
    expect_equal(length(results$background_scores$absolute_hits[[1]]),
                 length(background_set))

    results <- run_matrix_tsma(foreground_sets, background_set,
                               motifs = motif_db, cache = TRUE)
    expect_equal(length(results$foreground_scores), 2)
    expect_equal(length(results$foreground_scores[[1]]$total_sites[[1]]),
                 length(foreground_set1))
    expect_equal(length(results$foreground_scores[[2]]$total_sites[[1]]),
                 length(foreground_set2))
    expect_equal(as.numeric(results$foreground_scores[[2]]$total_sites[[1]]),
                 c(1, 4, 1, 4))
    expect_equal(as.numeric(results$foreground_scores[[2]]$df$total_sites),
                 sum(results$foreground_scores[[2]]$total_sites[[1]]))
    expect_equal(length(results$background_scores$total_sites[[1]]),
                 length(background_set))
    expect_equal(length(results$foreground_scores[[1]]$absolute_hits[[1]]),
                 length(foreground_set1))
    expect_equal(length(results$foreground_scores[[2]]$absolute_hits[[1]]),
                 length(foreground_set2))
    expect_equal(length(results$background_scores$absolute_hits[[1]]),
                 length(background_set))
})

test_that("run_matrix_spma", {
    background_df <- transite:::ge$background_df
    background_df <- dplyr::arrange(background_df, value)
    background_seqs <- gsub("T", "U", background_df$seq)
    names(background_seqs) <- paste0(background_df$refseq, "|",
                                     background_df$seq_type)

    n_bins <- 10
    results <- run_matrix_spma(background_seqs,
                               motifs = get_motif_by_id("M178_0.6"),
                               n_bins = n_bins,
                               max_fg_permutations = 10000, cache = FALSE)

    expect_equal(length(results$foreground_scores), n_bins)
    expect_equal(length(results$background_scores$total_sites[[1]]),
                 length(background_seqs))
    expect_equal(length(results$background_scores$absolute_hits[[1]]),
                 length(background_seqs))

    results <- run_matrix_spma(background_seqs,
                               sorted_transcript_values = background_df$value,
                               motifs = get_motif_by_id("M043_0.6"),
                               n_bins = n_bins,
                               max_fg_permutations = 10000, cache = FALSE)

    expect_equal(length(results$foreground_scores), n_bins)
    expect_equal(length(results$background_scores$total_sites[[1]]),
                 length(background_seqs))
    expect_equal(length(results$background_scores$absolute_hits[[1]]),
                 length(background_seqs))

    results <- run_matrix_spma(background_seqs,
                               motifs = get_motif_by_id("M178_0.6"),
                               n_bins = n_bins,
                               max_fg_permutations = 10000, cache = TRUE)

    expect_equal(length(results$foreground_scores), n_bins)
    expect_equal(length(results$background_scores$total_sites[[1]]),
                 length(background_seqs))
    expect_equal(length(results$background_scores$absolute_hits[[1]]),
                 length(background_seqs))
})

test_that("run_kmer_tsma", {
    foreground_set1 <- c(
        "CAACAGCUGUGUGCUUAAUU", "CAGUCAAGUGUUGGACUCC",
        "UGUGUUCUUUGGGGAAUUGUGUU",
        "UCAUUUUAGUGGUGUUAAA", "AAUUGGUGUCUGGAUACUUCCUGUGUUCUGUACAUUGUGUU",
        "AUCAAAUUA", "AGAU", "GACACUUAAAGAUCCU",
        "UAGCAUUAACUUAAUGUGUGGU", "UGUGGUAUGGA", "GAAGAGUGCUCA",
        "UGUGUGAUAGAC", "AGUUC", "CCAGUAAUGUGGU"
    )
    foreground_set2 <- c("UUAUUUA", "AUCCUUUACAUGUGGG", "UUUGUGGGUUUUU",
                         "UUUCAUGUGGGUCAUUGUGGGG")
    foreground_sets <- list(foreground_set1, foreground_set2)
    background_set <- unique(c(foreground_set1, foreground_set2, c(
        "CCACACAC", "CUCAUUGGAG", "ACUUUGGGACA", "CAGGUCAGCA",
        "CCACACCGG", "GUCAUCAGU", "GUCAGUCC", "CAGGUCAGGGGCA"
    )))

    motif_db <- get_motif_by_id("M178_0.6")
    results <- run_kmer_tsma(foreground_sets, background_set,
                             motifs = motif_db)

    expect_equal(nrow(results[[1]]$enrichment_df), 4^6)
    expect_equal(nrow(results[[1]]$motif_df), 1)
    expect_equal(nrow(results[[1]]$motif_kmers_dfs[[1]]),
                 length(get_hexamers(motif_db[[1]])))
})

test_that("run_kmer_spma", {
    background_df <- transite:::ge$background_df
    background_df <- dplyr::arrange(background_df, value)
    background_seqs <- gsub("T", "U", background_df$seq)
    names(background_seqs) <- paste0(background_df$refseq, "|",
                                     background_df$seq_type)

    n_bins <- 10
    expect_warning(results <- run_kmer_spma(background_seqs,
                                            motifs = get_motif_by_id("M178_0.6"),
                                            n_bins = n_bins,
                                            fg_permutations = 10))

    expect_equal(length(results$foreground_scores), n_bins)
    expect_equal(nrow(results$spectrum_info_df), 1)
    expect_equal(ncol(results$spectrum_info_df), 17)
    expect_equal(length(results$spectrum_plots), 1)
    expect_equal(length(results$classifier_scores), 1)
    expect_equal(length(results$classifier_scores[[1]]), 3)
    expect_equal(sum(results$classifier_scores[[1]]), 0)
})
