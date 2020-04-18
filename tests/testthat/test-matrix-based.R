context("matrix-based functions")

test_that("score_transcripts", {
    foreground_set <- c(
        "CAACAGCCUUAAUU", "CAGUCAAGACUCC", "CUUUGGGGAAU",
        "UCAUUUUAUUAAA", "AAUUGGUGUCUGGAUACUUCCCUGUACAU",
        "AUCAAAUUA", "AGAU", "GACACUUAAAGAUCCU",
        "UAGCAUUAACUUAAUG", "AUGGA", "GAAGAGUGCUCA",
        "AUAGAC", "AGUUC", "CCAGUAA"
    )
    expect_error(score_transcripts(foreground_set,
                                   cache = paste0(tempdir(), "/cache/")))
    names(foreground_set) <- c(
        "NM_1_DUMMY|3UTR", "NM_2_DUMMY|3UTR", "NM_3_DUMMY|3UTR",
        "NM_4_DUMMY|3UTR", "NM_5_DUMMY|3UTR", "NM_6_DUMMY|3UTR",
        "NM_7_DUMMY|3UTR", "NM_8_DUMMY|3UTR", "NM_9_DUMMY|3UTR",
        "NM_10_DUMMY|3UTR", "NM_11_DUMMY|3UTR", "NM_12_DUMMY|3UTR",
        "NM_13_DUMMY|3UTR", "NM_14_DUMMY|3UTR"
    )
    expect_error(score_transcripts(foreground_set, max_hits = 0))

    scores_multi_threaded <- score_transcripts(foreground_set,
                                               max_hits = 10,
                                               threshold_method = "relative",
                                               threshold_value = 0.9,
                                               n_cores = 2)
    scores_single_threaded <- score_transcripts(foreground_set,
                                                max_hits = 10,
                                                threshold_method = "relative",
                                                threshold_value = 0.9,
                                                n_cores = 1)
    expect_equal(scores_single_threaded, scores_multi_threaded)
})
