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

    max_hits_3 <- score_transcripts(foreground_set, max_hits = 3)
    expect_equal(ncol(max_hits_3$df), 9)
    max_hits_4 <- score_transcripts(foreground_set, max_hits = 4)
    expect_equal(ncol(max_hits_4$df), 10)
    max_hits_5 <- score_transcripts(foreground_set, max_hits = 5)
    expect_equal(ncol(max_hits_5$df), 11)
    max_hits_7 <- score_transcripts(foreground_set, max_hits = 7)
    expect_equal(ncol(max_hits_7$df), 13)
    max_hits_8 <- score_transcripts(foreground_set, max_hits = 8)
    expect_equal(ncol(max_hits_8$df), 14)
    max_hits_9 <- score_transcripts(foreground_set, max_hits = 9)
    expect_equal(ncol(max_hits_9$df), 15)

    expect_warning(score_transcripts(foreground_set,
                                     max_hits = 2,
                                     threshold_method = "relative",
                                     threshold_value = -0.1),
                   "threshold score below zero")

    skip_if(Sys.getenv("R_ARCH") == "/i386", "multi-core not available on i386")
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
