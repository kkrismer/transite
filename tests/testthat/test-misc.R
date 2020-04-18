context("miscellaneous functions")

test_that("create_matrix_motif", {
    custom_motif <- create_matrix_motif("custom_motif", "RBP1",
                                        transite:::toy_motif_matrix,
                                        "HITS-CLIP", "Homo sapiens", "user")
    expect_equal(transite::get_id(custom_motif), "custom_motif")
    expect_equal(transite::get_rbps(custom_motif), "RBP1")
    expect_equal(transite::get_motif_matrix(custom_motif),
                 transite:::toy_motif_matrix)
    expect_equal(transite::get_type(custom_motif), "HITS-CLIP")
    expect_equal(transite::get_species(custom_motif), "Homo sapiens")
    expect_equal(transite::get_source(custom_motif), "user")
})

test_that("create_kmer_motif", {
    custom_motif <- create_kmer_motif("custom_motif", "RBP1",
                                      c("AAAAAAA", "CAAAAAA"), "HITS-CLIP",
                                      "Homo sapiens", "user")

    expect_equal(transite::get_id(custom_motif), "custom_motif")
    expect_equal(transite::get_rbps(custom_motif), "RBP1")
    expect_equal(transite::get_hexamers(custom_motif), c("AAAAAA", "CAAAAA"))
    expect_equal(transite::get_type(custom_motif), "HITS-CLIP")
    expect_equal(transite::get_species(custom_motif), "Homo sapiens")
    expect_equal(transite::get_source(custom_motif), "user")
})

test_that("get_motifs_meta_info", {
    df <- get_motifs_meta_info()

    expect_equal(nrow(df), 174)
    expect_equal(ncol(df), 7)
})

test_that("get_motif_by_rbp", {
    motifs <- get_motif_by_rbp(c("ELAVL1", "ELAVL2"))

    expect_equal(length(motifs), 8)
})

test_that("get_ppm", {
    ppm <- get_ppm(get_motif_by_id("M178_0.6")[[1]])

    expect_equal(colnames(ppm), c("A", "C", "G", "U"))
    expect_equal(nrow(ppm), 7)
    expect_equal(sum(ppm), 7, tolerance = 0.0002)
})

test_that("generate_iupac_by_matrix", {
    expect_equal(generate_iupac_by_matrix(
        get_motif_matrix(get_motif_by_id("M178_0.6")[[1]])),
                 "UGUGKKG")
})

test_that("generate_iupac_by_kmers", {
    expect_equal(generate_iupac_by_kmers(c("AACCAA", "AACCGG", "CACCGA")),
        "MACCRR")
    expect_equal(generate_iupac_by_kmers(c("AAAAAA", "AAAAAA", "AAAAAA")),
                 "AAAAAA")
})

test_that("generate_kmers_from_iupac", {
    expect_equal(generate_kmers_from_iupac(
        get_iupac(get_motif_by_id("M178_0.6")[[1]]), k = 7),
                 c("UGUGGGG", "UGUGUGG", "UGUGGUG", "UGUGUUG"))
    expect_equal(generate_kmers_from_iupac(
        get_iupac(get_motif_by_id("M178_0.6")[[1]]), k = 8),
                 c("AUGUGGGG", "CUGUGGGG", "GUGUGGGG", "UUGUGGGG", "AUGUGUGG",
                   "CUGUGUGG", "GUGUGUGG", "UUGUGUGG", "AUGUGGUG", "CUGUGGUG",
                   "GUGUGGUG", "UUGUGGUG", "AUGUGUUG", "CUGUGUUG", "GUGUGUUG",
                   "UUGUGUUG", "UGUGGGGA", "UGUGUGGA", "UGUGGUGA", "UGUGUUGA",
                   "UGUGGGGC", "UGUGUGGC", "UGUGGUGC", "UGUGUUGC", "UGUGGGGG",
                   "UGUGUGGG", "UGUGGUGG", "UGUGUUGG", "UGUGGGGU", "UGUGUGGU",
                   "UGUGGUGU", "UGUGUUGU"))
})
