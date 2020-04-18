context("k-mer-based functions")

test_that("generate_kmers", {
    seqs <- c(
        "CAACAGCTGTGTGCTTAATT", "CAGTCAAGTGTTGGACTCC",
        "TGTGTTCTTTGGGGAATTGTGTT",
        "TCATTTTAGTGGTGTTAAA", "AATTGGTGTCTGGATACTTCCTGTGTTCTGTACATTGTGTT",
        "ATCAAATTA", "AGAT", "GACACTTAAAGATCCT",
        "TAGCATTAACTTAATGTGTGGT", "TGTGGTATGGA", "GAAGAGTGCTCA",
        "TGTGTGATAGAC", "AGTTC", "CCAGTAATGTGGT"
    )
    kmer_counts <- generate_kmers(seqs, 6)
    expect_equal(length(kmer_counts), 4^6)
})

test_that("calculate_kmer_enrichment", {
    foreground_set1 <- c(
        "CAACAGCCUUAAUU", "CAGUCAAGACUCC", "CUUUGGGGAAU",
        "UCAUUUUAUUAAA", "AAUUGGUGUCUGGAUACUUCCCUGUACAU",
        "AUCAAAUUA", "AGAU", "GACACUUAAAGAUCCU",
        "UAGCAUUAACUUAAUG", "AUGGA", "GAAGAGUGCUCA",
        "AUAGAC", "AGUUC", "CCAGUAA"
    )
    foreground_set2 <- c("UUAUUUA", "AUCCUUUACA", "UUUUUUU", "UUUCAUCAUU")
    foreground_sets <- list(foreground_set1, foreground_set2)
    background_set <- c(foreground_set1, foreground_set2,
                        "CCACACAC", "CUCAUUGGAG", "ACUUUGGGACA", "CAGGUCAGCA")

    single_threaded_results <- calculate_kmer_enrichment(foreground_sets,
                                                         background_set, 6,
                                                         n_cores = 1)
    multi_threaded_results <- calculate_kmer_enrichment(foreground_sets,
                                                         background_set, 6,
                                                         n_cores = 2)
    expect_equal(single_threaded_results, multi_threaded_results)
})

test_that("draw_volcano_plot", {
    motif <- get_motif_by_id("951_12324455")
    p <- draw_volcano_plot(transite:::kmers_enrichment,
                           get_hexamers(motif[[1]]),
                           get_rbps(motif[[1]]), show_legend = FALSE)
    expect_equal(class(p), c("gg", "ggplot"))
})
