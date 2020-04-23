context("spectrum plot functions")

test_that("score_spectrum", {
    set.seed(0)
    sp <- score_spectrum(runif(n = 40, min = -1, max = 1), max_model_degree = 1)
    expect_equal(transite::get_adj_r_squared(sp), 0)
    expect_equal(transite::get_model_degree(sp), 0)

    sp <- score_spectrum(runif(n = 40, min = -1, max = 1),
                         p_values = runif(n = 40, min = 0, max = 1),
                         max_model_degree = 1)
    expect_equal(transite::get_adj_r_squared(sp), 0)
    expect_equal(transite::get_model_degree(sp), 0)

    log_fold_change <- log(runif(n = 1000, min = 0, max = 1) /
                               runif(n = 1000, min = 0, max = 1))
    sp <- score_spectrum(runif(n = 40, min = -1, max = 1),
                   sorted_transcript_values = sort(log_fold_change),
                   max_model_degree = 1)
    expect_equal(transite::get_adj_r_squared(sp), 0)
    expect_equal(transite::get_model_degree(sp), 0)
})
