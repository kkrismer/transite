context("combine p-values")

test_that("p_combine", {
    expect_equal(p_combine(c(0.01, 0.05, 0.5), method = "fisher")$p_value,
                 0.01092242, tolerance = 0.0002)
    expect_equal(p_combine(c(0.01, 0.05, 0.5), method = "SL")$p_value,
                 0.0109305, tolerance = 0.0002)
    expect_equal(p_combine(c(0.01, 0.05, 0.5), method = "MG")$p_value,
                 0.01004926, tolerance = 0.0002)
    expect_equal(p_combine(c(0.01, 0.05, 0.5), method = "tippett")$p_value,
                 0.029701, tolerance = 0.0002)

    expect_equal(p_combine(rep(0.05, 5), method = "fisher")$p_value,
                 0.0008705185, tolerance = 0.0002)
    expect_equal(p_combine(rep(0.05, 5), method = "SL")$p_value,
                 0.0001175329, tolerance = 0.0002)
    expect_equal(p_combine(rep(0.05, 5), method = "MG")$p_value,
                 0.000380212, tolerance = 0.0002)
    expect_equal(p_combine(rep(0.05, 5), method = "tippett")$p_value,
                 0.2262191, tolerance = 0.0002)

    expect_error(p_combine(c(0.01, 1.1)))
    expect_error(p_combine(c(0.01, 0.23), method = "xyc"))
    expect_warning(p_combine(c(), method = "fisher"))
    expect_warning(p_combine(c(), method = "SL"))
    expect_warning(p_combine(c(), method = "MG"))
    expect_warning(p_combine(c(), method = "tippett"))

    expect_error(p_combine(c(0.01, 0.23), method = "SL", w = c(0.2)))
})
