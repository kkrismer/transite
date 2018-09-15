#' @title Matrix-based Transcript Set Motif Analysis
#'
#' @description
#' Calculates motif enrichment in foreground sets versus a background set using position
#' weight matrices to identify putative binding sites
#'
#' @param foreground.sets a list of named character vectors of foreground sequences
#' (only containing upper case characters A, C, G, T), where the names are RefSeq identifiers
#' and sequence
#' type qualifiers (\code{"3UTR"}, \code{"5UTR"}, \code{"mRNA"}), e.g. \code{"NM_010356|3UTR"}
#' @param background.set a named character vector of background sequences (naming follows same
#' rules as foreground set sequences)
#' @inheritParams scoreTranscripts
#' @inheritParams calculateMotifEnrichment
#'
#' @return A list with the following components:
#' \tabular{rl}{
#'   \code{foreground.scores} \tab the result of \code{\link{scoreTranscripts}} for the foreground
#'   sets\cr
#'   \code{background.scores} \tab the result of \code{\link{scoreTranscripts}} for the background
#'   set\cr
#'   \code{enrichment.dfs} \tab a list of data frames, returned by
#'   \code{\link{calculateMotifEnrichment}}
#' }
#'
#' @details
#' Motif transcript set analysis can be used to identify RNA binding proteins, whose targets are
#' significantly overrepresented or underrepresented in certain sets of transcripts.
#'
#' The aim of Transcript Set Motif Analysis (TSMA) is to identify the overrepresentation
#' and underrepresentation of potential RBP targets (binding sites) in a set (or sets) of
#' sequences, i.e., the foreground set, relative to the entire population of sequences.
#' The latter is called background set, which can be composed of all sequences of the genes
#' of a microarray platform or all sequences of an organism or any other meaningful
#' superset of the foreground sets.
#'
#' The matrix-based approach skips the \emph{k}-merization step of the \emph{k}-mer-based approach
#' and instead scores the transcript sequence as a whole with a position specific scoring matrix.
#'
#' For each sequence in foreground and background sets and each sequence motif,
#' the scoring algorithm evaluates the score for each sequence position. Positions with
#' a relative score greater than a certain threshold are considered hits, i.e.,
#' putative binding sites.
#'
#' By scoring all sequences in foreground and background sets, a hit count for each motif and
#' each set is obtained, which is used to calculate enrichment values and associated p-values
#' in the same way in which motif-compatible hexamer enrichment values are calculated in
#' the k -mer-based approach. P-values are adjusted with one of the available adjustment methods.
#'
#' An advantage of the matrix-based approach is the possibility of detecting clusters of
#' binding sites. This can be done by counting regions with many hits using positional
#' hit information or by simply applying a hit count threshold per sequence, e.g., only
#' sequences with more than some number of hits are considered. Homotypic clusters of RBP
#' binding sites may play a similar role as clusters of transcription factors.
#'
#' @examples
#' # define simple sequence sets for foreground and background
#' foreground.set1 <- c(
#'   "CAACAGCCTTAATT", "CAGTCAAGACTCC", "CTTTGGGGAAT",
#'   "TCATTTTATTAAA", "AATTGGTGTCTGGATACTTCCCTGTACAT",
#'   "ATCAAATTA", "AGAT", "GACACTTAAAGATCCT",
#'   "TAGCATTAACTTAATG", "ATGGA", "GAAGAGTGCTCA",
#'   "ATAGAC", "AGTTC", "CCAGTAA"
#' )
#' names(foreground.set1) <- c(
#'   "NM_1_DUMMY|3UTR", "NM_2_DUMMY|3UTR", "NM_3_DUMMY|3UTR",
#'   "NM_4_DUMMY|3UTR", "NM_5_DUMMY|3UTR", "NM_6_DUMMY|3UTR", "NM_7_DUMMY|3UTR",
#'   "NM_8_DUMMY|3UTR", "NM_9_DUMMY|3UTR", "NM_10_DUMMY|3UTR", "NM_11_DUMMY|3UTR",
#'   "NM_12_DUMMY|3UTR", "NM_13_DUMMY|3UTR", "NM_14_DUMMY|3UTR"
#' )
#' 
#' foreground.set2 <- c("TTATTTA", "ATCCTTTACA", "TTTTTTT", "TTTCATCATT")
#' names(foreground.set2) <- c(
#'   "NM_15_DUMMY|3UTR", "NM_16_DUMMY|3UTR", "NM_17_DUMMY|3UTR",
#'   "NM_18_DUMMY|3UTR"
#' )
#' 
#' foreground.sets <- list(foreground.set1, foreground.set2)
#' 
#' background.set <- c(
#'   "CAACAGCCTTAATT", "CAGTCAAGACTCC", "CTTTGGGGAAT",
#'   "TCATTTTATTAAA", "AATTGGTGTCTGGATACTTCCCTGTACAT",
#'   "ATCAAATTA", "AGAT", "GACACTTAAAGATCCT",
#'   "TAGCATTAACTTAATG", "ATGGA", "GAAGAGTGCTCA",
#'   "ATAGAC", "AGTTC", "CCAGTAA",
#'   "TTATTTA", "ATCCTTTACA", "TTTTTTT", "TTTCATCATT",
#'   "CCACACAC", "CTCATTGGAG", "ACTTTGGGACA", "CAGGTCAGCA"
#' )
#' names(background.set) <- c(
#'   "NM_1_DUMMY|3UTR", "NM_2_DUMMY|3UTR", "NM_3_DUMMY|3UTR",
#'   "NM_4_DUMMY|3UTR", "NM_5_DUMMY|3UTR", "NM_6_DUMMY|3UTR", "NM_7_DUMMY|3UTR",
#'   "NM_8_DUMMY|3UTR", "NM_9_DUMMY|3UTR", "NM_10_DUMMY|3UTR", "NM_11_DUMMY|3UTR",
#'   "NM_12_DUMMY|3UTR", "NM_13_DUMMY|3UTR", "NM_14_DUMMY|3UTR", "NM_15_DUMMY|3UTR",
#'   "NM_16_DUMMY|3UTR", "NM_17_DUMMY|3UTR", "NM_18_DUMMY|3UTR", "NM_19_DUMMY|3UTR",
#'   "NM_20_DUMMY|3UTR", "NM_21_DUMMY|3UTR", "NM_22_DUMMY|3UTR"
#' )
#' 
#' # run cached version of TSMA with all Transite motifs (recommended)
#' results <- runMatrixTSMA(foreground.sets, background.set)
#' \dontrun{
#' # define exemplary sequence sets for foreground and background
#' foreground1.df <- ge$foreground1
#' foreground.set1 <- foreground1.df$seq
#' names(foreground.set1) <- paste0(foreground1.df$refseq, "|", foreground1.df$seq.type)
#' 
#' foreground2.df <- ge$foreground2
#' foreground.set2 <- foreground2.df$seq
#' names(foreground.set2) <- paste0(foreground2.df$refseq, "|", foreground2.df$seq.type)
#' 
#' foreground.sets <- list(foreground.set1, foreground.set2)
#' 
#' background.df <- ge$background
#' background.set <- background.df$seq
#' names(background.set) <- paste0(background.df$refseq, "|", background.df$seq.type)
#' 
#' # run cached version of TSMA with all Transite motifs (recommended)
#' results <- runMatrixTSMA(foreground.sets, background.set)
#' 
#' # run uncached version of TSMA with all Transite motifs
#' results <- runMatrixTSMA(foreground.sets, background.set, cache = FALSE)
#' 
#' # run TSMA with a subset of Transite motifs
#' results <- runMatrixTSMA(foreground.sets, background.set, motifs = getMotifByRBP("ELAVL1"))
#' 
#' # run TSMA with user-defined motif
#' toy.motif <- createMatrixMotif(
#'   "toy.motif", "example RBP", toy.motif.matrix,
#'   "example type", "example species", "user"
#' )
#' results <- runMatrixTSMA(foreground.sets, background.set, motifs = list(toy.motif))
#' }
#' 
#' @family TSMA functions
#' @family matrix functions
#' @importFrom dplyr arrange
#' @export
runMatrixTSMA <-
  function(foreground.sets,
             background.set,
             motifs = NULL,
             max.hits = 5,
             threshold.method = "p.value",
             threshold.value = 0.25^6,
             max.fg.permutations = 1000000,
             min.fg.permutations = 1000,
             e = 5,
             p.adjust.method = "BH",
             n.cores = 1,
             cache = paste0(getwd(), "/sc/")) {
    # avoid CRAN note
    adj.p.value <- p.value <- NULL

    i <- 1
    foreground.scores <-
      lapply(foreground.sets, function(foreground.set) {
        result <-
          scoreTranscripts(
            foreground.set,
            motifs = motifs,
            max.hits = max.hits,
            threshold.method = threshold.method,
            threshold.value = threshold.value,
            n.cores = n.cores,
            cache = cache
          )
        message(paste0("scored transcripts in foreground set ", i))
        i <<- i + 1
        return(result)
      })

    background.scores <-
      scoreTranscripts(
        background.set,
        motifs = motifs,
        max.hits = max.hits,
        threshold.method = threshold.method,
        threshold.value = threshold.value,
        n.cores = n.cores,
        cache = cache
      )
    message("scored transcripts in background set")

    i <- 1
    enrichment.dfs <-
      lapply(foreground.scores, function(scores.per.condition) {
        enrichment.df <-
          calculateMotifEnrichment(
            scores.per.condition$df,
            background.scores$df,
            background.scores$total.sites,
            background.scores$absolute.hits,
            length(foreground.sets[[i]]),
            max.fg.permutations = max.fg.permutations,
            min.fg.permutations = min.fg.permutations,
            e = e,
            p.adjust.method = p.adjust.method
          )
        enrichment.df <-
          dplyr::arrange(enrichment.df, adj.p.value, p.value)
        message(paste0("calculated enrichment for foreground set ", i))
        i <<- i + 1
        return(enrichment.df)
      })
    return(
      list(
        foreground.scores = foreground.scores,
        background.scores = background.scores,
        enrichment.dfs = enrichment.dfs
      )
    )
  }

#' @title Matrix-based Spectrum Motif Analysis
#'
#' @description
#' SPMA helps to illuminate the relationship between RBP binding evidence and the transcript
#' sorting criterion, e.g., fold change between treatment and control samples.
#'
#' @inheritParams runMatrixTSMA
#' @inheritParams subdivideData
#' @inheritParams scoreSpectrum
#'
#' @return A list with the following components:
#' \tabular{rl}{
#'   \code{foreground.scores} \tab the result of \code{\link{scoreTranscripts}} for the foreground
#'   sets (the bins)\cr
#'   \code{background.scores} \tab the result of \code{\link{scoreTranscripts}} for the background
#'   set\cr
#'   \code{enrichment.dfs} \tab a list of data frames, returned by
#'   \code{\link{calculateMotifEnrichment}}\cr
#'   \code{spectrum.info.df} \tab a data frame with the SPMA results\cr
#'   \code{spectrum.plots} \tab a list of spectrum plots, as generated by \code{\link{scoreSpectrum}}\cr
#'   \code{classifier.scores} \tab a list of classifier scores, as returned by
#'   \code{\link{spectrumClassifier}}
#' }
#'
#' @details
#' In order to investigate how motif targets are distributed across a spectrum of
#' transcripts (e.g., all transcripts of a platform, ordered by fold change),
#' Spectrum Motif Analysis visualizes the gradient of RBP binding evidence
#' across all transcripts.
#'
#' The matrix-based approach skips the \emph{k}-merization step of the \emph{k}-mer-based approach
#' and instead scores the transcript sequence as a whole with a position specific scoring matrix.
#'
#' For each sequence in foreground and background sets and each sequence motif,
#' the scoring algorithm evaluates the score for each sequence position. Positions with
#' a relative score greater than a certain threshold are considered hits, i.e.,
#' putative binding sites.
#'
#' By scoring all sequences in foreground and background sets, a hit count for each motif and
#' each set is obtained, which is used to calculate enrichment values and associated p-values
#' in the same way in which motif-compatible hexamer enrichment values are calculated in
#' the \emph{k}-mer-based approach. P-values are adjusted with one of the
#' available adjustment methods.
#'
#' An advantage of the matrix-based approach is the possibility of detecting clusters of
#' binding sites. This can be done by counting regions with many hits using positional
#' hit information or by simply applying a hit count threshold per sequence, e.g., only
#' sequences with more than some number of hits are considered. Homotypic clusters of RBP
#' binding sites may play a similar role as clusters of transcription factors.
#'
#' @examples
#' \dontrun{
#' # exemplary data set
#' background.df <- ge$background
#' # sort sequences by signal-to-noise ratio
#' background.df <- dplyr::arrange(background.df, value)
#' # character vector of named sequences
#' background.set <- background.df$seq
#' names(background.set) <- paste0(background.df$refseq, "|", background.df$seq.type)
#' 
#' results <- runMatrixSPMA(background.set)
#' }
#' 
#' @family SPMA functions
#' @family matrix functions
#' @importFrom stats p.adjust
#' @importFrom dplyr filter
#' @export
runMatrixSPMA <-
  function(background.set,
             motifs = NULL,
             n.bins = 40,
             max.model.degree = 1,
             max.cs.permutations = 10000000,
             min.cs.permutations = 5000,
             max.hits = 5,
             threshold.method = "p.value",
             threshold.value = 0.25^6,
             max.fg.permutations = 1000000,
             min.fg.permutations = 1000,
             e = 5,
             p.adjust.method = "BH",
             n.cores = 1,
             cache = paste0(getwd(), "/sc/")) {
    # avoid CRAN note
    motif.id <-
      motif.rbps <-
      adj.r.squared <- degree <- residuals <- slope <- NULL
    f.statistic <-
      f.statistic.p.value <- f.statistic.adj.p.value <- NULL
    consistency.score <-
      consistency.score.p.value <-
      consistency.score.adj.p.value <-
      consistency.score.n <- NULL
    n.significant <-
      n.very.significant <-
      n.extremely.significant <-
      aggregate.classifier.score <- NULL

    foreground.sets <- subdivideData(background.set, n.bins)

    results <- runMatrixTSMA(
      foreground.sets,
      background.set,
      motifs = motifs,
      max.hits = max.hits,
      threshold.method = threshold.method,
      threshold.value = threshold.value,
      max.fg.permutations = max.fg.permutations,
      min.fg.permutations = min.fg.permutations,
      e = e,
      p.adjust.method = p.adjust.method,
      n.cores = n.cores,
      cache = cache
    )

    if (length(results$enrichment.dfs) > 0) {
      enrichment.df <- do.call("rbind", results$enrichment.dfs)
      enrichment.df$adj.p.value <-
        stats::p.adjust(enrichment.df$p.value, method = p.adjust.method)

      motifs <- getMotifs()
      spectrum.info <- lapply(motifs, function(motif) {
        motif.data.df <- dplyr::filter(enrichment.df, motif.id == motif$id)
        values <- motif.data.df$enrichment
        values[values == 0] <-
          0.01 # avoid -Inf after taking the log
        score <-
          scoreSpectrum(
            log(values),
            motif.data.df$adj.p.value,
            max.model.degree = max.model.degree,
            max.cs.permutations = max.cs.permutations,
            min.cs.permutations = min.cs.permutations
          )

        n.significant <-
          sum(motif.data.df$adj.p.value <= 0.05)
        classifier.score <-
          spectrumClassifier(
            score$adj.r.squared,
            score$degree,
            score$slope,
            score$consistency.score.n,
            n.significant,
            n.bins
          )

        return(
          list(
            info = list(
              motif.id = motif$id,
              motif.rbps = paste(motif$rbps, collapse = ", "),
              adj.r.squared = score$adj.r.squared,
              degree = score$degree,
              residuals = score$residuals,
              slope = score$slope,
              f.statistic = score$f.statistic,
              f.statistic.p.value = score$f.statistic.p.value,
              consistency.score = score$consistency.score,
              consistency.score.p.value = score$consistency.score.p.value,
              consistency.score.n = score$consistency.score.n,
              n.significant = n.significant,
              n.very.significant = sum(motif.data.df$adj.p.value <= 0.01),
              n.extremely.significant = sum(motif.data.df$adj.p.value <= 0.001),
              aggregate.classifier.score = sum(classifier.score)
            ),
            spectrum.plot = score$plot,
            classifier.score = classifier.score
          )
        )
      })
      spectrum.info.df <-
        as.data.frame(do.call("rbind", lapply(spectrum.info, function(x)
          x$info)),
        stringsAsFactors = FALSE
        )
      spectrum.info.df$motif.id <-
        as.character(spectrum.info.df$motif.id)
      spectrum.info.df$motif.rbps <-
        as.character(spectrum.info.df$motif.rbps)
      spectrum.info.df$adj.r.squared <-
        as.numeric(spectrum.info.df$adj.r.squared)
      spectrum.info.df$degree <-
        as.integer(spectrum.info.df$degree)
      spectrum.info.df$residuals <-
        as.numeric(spectrum.info.df$residuals)
      spectrum.info.df$slope <-
        as.numeric(spectrum.info.df$slope)
      spectrum.info.df$f.statistic <-
        as.numeric(spectrum.info.df$f.statistic)
      spectrum.info.df$f.statistic.p.value <-
        as.numeric(spectrum.info.df$f.statistic.p.value)
      spectrum.info.df$consistency.score <-
        as.numeric(spectrum.info.df$consistency.score)
      spectrum.info.df$consistency.score.p.value <-
        as.numeric(spectrum.info.df$consistency.score.p.value)
      spectrum.info.df$consistency.score.n <-
        as.integer(spectrum.info.df$consistency.score.n)
      spectrum.info.df$n.significant <-
        as.integer(spectrum.info.df$n.significant)
      spectrum.info.df$n.very.significant <-
        as.integer(spectrum.info.df$n.very.significant)
      spectrum.info.df$n.extremely.significant <-
        as.integer(spectrum.info.df$n.extremely.significant)
      spectrum.info.df$aggregate.classifier.score <-
        as.integer(spectrum.info.df$aggregate.classifier.score)

      spectrum.info.df$f.statistic.adj.p.value <-
        stats::p.adjust(spectrum.info.df$f.statistic.p.value, method = p.adjust.method)
      spectrum.info.df$consistency.score.adj.p.value <-
        stats::p.adjust(spectrum.info.df$consistency.score.p.value, method = p.adjust.method)

      spectrum.info.df <-
        dplyr::select(
          spectrum.info.df,
          motif.id,
          motif.rbps,
          adj.r.squared,
          degree,
          residuals,
          slope,
          f.statistic,
          f.statistic.p.value,
          f.statistic.adj.p.value,
          consistency.score,
          consistency.score.p.value,
          consistency.score.adj.p.value,
          consistency.score.n,
          n.significant,
          n.very.significant,
          n.extremely.significant,
          aggregate.classifier.score
        )

      spectrum.plots <-
        lapply(spectrum.info, function(x)
          x$spectrum.plot)
      classifier.scores <-
        lapply(spectrum.info, function(x)
          x$classifier.score)
    } else {
      spectrum.info.df <- data.frame(
        motif.id = character(0),
        motif.rbps = character(0),
        adj.r.squared = numeric(0),
        degree = integer(0),
        residuals = numeric(0),
        slope = numeric(0),
        f.statistic = numeric(0),
        f.statistic.p.value = numeric(0),
        f.statistic.adj.p.value = numeric(0),
        consistency.score = numeric(0),
        consistency.score.p.value = numeric(0),
        consistency.score.adj.p.value = numeric(0),
        consistency.score.n,
        n.significant = integer(0),
        n.very.significant = integer(0),
        n.extremely.significant = integer(0),
        aggregate.classifier.score = integer(0)
      )
      spectrum.plots <- NULL
      warning("no sequences in any condition")
    }

    return(
      list(
        foreground.scores = results$foreground.scores,
        background.scores = results$background.scores,
        enrichment.dfs = results$enrichment.dfs,
        spectrum.info.df = spectrum.info.df,
        spectrum.plots = spectrum.plots,
        classifier.scores = classifier.scores
      )
    )
  }

#' @title \emph{k}-mer-based Transcript Set Motif Analysis
#'
#' @description
#' Calculates the enrichment of putative binding sites in foreground sets versus a background set
#' using \emph{k}-mers to identify putative binding sites
#'
#' @param motifs a list of motifs that is used to score the specified sequences.
#' If \code{is.null(motifs)} then all Transite motifs are used.
#' @param fg.permutations numer of foreground permutations
#' @param kmer.significance.threshold p-value threshold for significance, e.g., \code{0.05} or
#' \code{0.01} (used for volcano plots)
#' @param produce.plot if \code{TRUE} volcano plots and distribution plots are created
#' @param p.combining.method one of the following: Fisher (1932) (\code{"fisher"}), Stouffer (1949),
#' Liptak (1958) (\code{"SL"}), Mudholkar and George (1979) (\code{"MG"}), and Tippett (1931)
#' (\code{"tippett"}) (see \code{\link{pCombine}})
#' @inheritParams calculateKmerEnrichment
#'
#' @return A list of lists with the following components:
#' \tabular{rl}{
#'   \code{enrichment.df} \tab \cr
#'   \code{motif.df} \tab \cr
#'   \code{motif.kmers.dfs} \tab \cr
#'   \code{volcano.plots} \tab \cr
#'   \code{perm.test.plots} \tab \cr
#'   \code{enriched.kmers.combined.p.values} \tab \cr
#'   \code{depleted.kmers.combined.p.values} \tab
#' }
#'
#' @details
#' Motif transcript set analysis can be used to identify RNA binding proteins, whose targets are
#' significantly overrepresented or underrepresented in certain sets of transcripts.
#'
#' The aim of Transcript Set Motif Analysis (TSMA) is to identify the overrepresentation
#' and underrepresentation of potential RBP targets (binding sites) in a set (or sets) of
#' sequences, i.e., the foreground set, relative to the entire population of sequences.
#' The latter is called background set, which can be composed of all sequences of the genes
#' of a microarray platform or all sequences of an organism or any other meaningful
#' superset of the foreground sets.
#'
#' The \emph{k}-mer-based approach breaks the sequences of foreground and background sets into
#' \emph{k}-mers and calculates the enrichment on a \emph{k}-mer level. In this case, motifs are
#' not represented as position weight matrices, but as lists of \emph{k}-mers.
#'
#' Statistically significantly enriched or depleted \emph{k}-mers are then used to
#' calculate a score for each RNA-binding protein, which quantifies its
#' target overrepresentation.
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
#' background.set <- c(
#'   "CAACAGCCUUAAUU", "CAGUCAAGACUCC", "CUUUGGGGAAU",
#'   "UCAUUUUAUUAAA", "AAUUGGUGUCUGGAUACUUCCCUGUACAU",
#'   "AUCAAAUUA", "AGAU", "GACACUUAAAGAUCCU",
#'   "UAGCAUUAACUUAAUG", "AUGGA", "GAAGAGUGCUCA",
#'   "AUAGAC", "AGUUC", "CCAGUAA",
#'   "UUAUUUA", "AUCCUUUACA", "UUUUUUU", "UUUCAUCAUU",
#'   "CCACACAC", "CUCAUUGGAG", "ACUUUGGGACA", "CAGGUCAGCA"
#' )
#' 
#' # run k-mer based TSMA with all Transite motifs
#' results <- runKmerTSMA(foreground.sets, background.set)
#' \dontrun{
#' # define exemplary sequence sets for foreground and background
#' foreground.set1 <- gsub("T", "U", ge$foreground1$seq)
#' foreground.set2 <- gsub("T", "U", ge$foreground2$seq)
#' foreground.sets <- list(foreground.set1, foreground.set2)
#' background.set <- gsub("T", "U", ge$background$seq)
#' 
#' # run TSMA with all Transite motifs
#' results <- runKmerTSMA(foreground.sets, background.set)
#' 
#' # run TSMA with a subset of Transite motifs
#' results <- runKmerTSMA(foreground.sets, background.set, motifs = getMotifByRBP("ELAVL1"))
#' 
#' # run TSMA with user-defined motif
#' toy.motif <- createKmerMotif(
#'   "toy.motif", "example RBP",
#'   c("AACCGG", "AAAACG", "AACACG"), "example type", "example species", "user"
#' )
#' results <- runMatrixTSMA(foreground.sets, background.set, motifs = list(toy.motif))
#' }
#' 
#' @family TSMA functions
#' @family \emph{k}-mer functions
#' @importFrom stats p.adjust
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @export
runKmerTSMA <-
  function(foreground.sets,
             background.set,
             motifs = NULL,
             k = 6,
             fg.permutations = 5000,
             kmer.significance.threshold = 0.01,
             produce.plot = TRUE,
             p.adjust.method = "BH",
             p.combining.method = "fisher",
             n.cores = 1) {
    # avoid CRAN note
    kmer <- enrichment <- p.value <- adj.p.value <- NULL

    if (is.null(motifs)) {
      motifs <- getMotifs()
    }

    enrichment.result <-
      calculateKmerEnrichment(
        foreground.sets,
        background.set,
        k,
        p.adjust.method = p.adjust.method,
        n.cores = n.cores
      )
    message("calculated enrichment for all foreground sets")

    i <- 1
    if (!is.null(enrichment.result)) {
      set.sizes <- unique(unlist(lapply(foreground.sets, length)))
      random.enrichments <- new.env()
      for (set.size in set.sizes) {
        if (set.size < 10) {
          adapted.fg.permutations <- min(fg.permutations, 100)
        } else if (set.size < 30) {
          adapted.fg.permutations <- min(fg.permutations, 500)
        } else {
          adapted.fg.permutations <- fg.permutations
        }
        assign(
          as.character(set.size),
          generatePermutedEnrichments(
            set.size,
            background.set,
            k,
            n.permutations = adapted.fg.permutations,
            n.cores = n.cores
          ),
          envir = random.enrichments
        )
      }
      message("calculated enrichment for all permuted sets")

      motif.result <-
        lapply(seq_len(length(foreground.sets)), function(i) {
          kmers.df <- enrichment.result$dfs[[i]]
          kmers.df$kmer <- enrichment.result$kmers

          foreground.result <-
            lapply(motifs, function(motif) {
              if (k == 6) {
                rbp.kmers <- motif$hexamers
              } else if (k == 7) {
                rbp.kmers <- motif$heptamers
              }

              idx <-
                which(enrichment.result$kmers %in% rbp.kmers)
              motif.kmers.df <- kmers.df[idx, ]
              motif.kmers.df <-
                dplyr::select(
                  motif.kmers.df,
                  kmer,
                  enrichment,
                  p.value,
                  adj.p.value
                )

              geo.mean <-
                geometricMean(motif.kmers.df$enrichment)
              perm.test <-
                permTestGeometricMean(
                  geo.mean,
                  rbp.kmers,
                  get(
                    as.character(length(foreground.sets[[i]])),
                    envir = random.enrichments,
                    inherits = FALSE
                  ),
                  alternative = "two.sided",
                  conf.level = 0.95,
                  produce.plot = produce.plot
                )

              motif.kmers.enriched.df <-
                dplyr::filter(motif.kmers.df, enrichment > 1)
              enriched.kmers.combined.p.value <-
                pCombine(motif.kmers.enriched.df$adj.p.value,
                  method = p.combining.method
                )

              motif.kmers.depleted.df <-
                dplyr::filter(motif.kmers.df, enrichment < 1)
              depleted.kmers.combined.p.value <-
                pCombine(motif.kmers.depleted.df$adj.p.value,
                  method = p.combining.method
                )

              if (produce.plot) {
                volcano.plot <-
                  drawVolcanoPlot(
                    kmers.df,
                    rbp.kmers,
                    paste(motif$rbps, collapse = ", "),
                    kmer.significance.threshold
                  )
              } else {
                volcano.plot <- NULL
              }

              return(
                list(
                  df = list(
                    motif.id = motif$id,
                    motif.rbps = paste0(motif$rbps, collapse = ", "),
                    geo.mean.enrichment = geo.mean,
                    p.value.estimate = perm.test$p.value.estimate,
                    cont.int.low = perm.test$conf.int[1],
                    cont.int.high = perm.test$conf.int[2],
                    enriched.kmers.combined.p.value = enriched.kmers.combined.p.value$p.value,
                    depleted.kmers.combined.p.value = depleted.kmers.combined.p.value$p.value
                  ),
                  motif.kmers.df = motif.kmers.df,
                  volcano.plot = volcano.plot,
                  perm.test.plot = perm.test$plot,
                  enriched.kmers.combined.p.value = enriched.kmers.combined.p.value,
                  depleted.kmers.combined.p.value = depleted.kmers.combined.p.value
                )
              )
            })

          df <-
            as.data.frame(do.call(
              "rbind",
              lapply(foreground.result, function(x)
                x$df)
            ), stringsAsFactors = FALSE)
          df$motif.id <- as.character(df$motif.id)
          df$motif.rbps <- as.character(df$motif.rbps)
          df$geo.mean.enrichment <-
            as.numeric(df$geo.mean.enrichment)
          df$p.value.estimate <-
            as.numeric(df$p.value.estimate)
          df$cont.int.low <- as.numeric(df$cont.int.low)
          df$cont.int.high <- as.numeric(df$cont.int.high)
          df$enriched.kmers.combined.p.value <-
            as.numeric(df$enriched.kmers.combined.p.value)
          df$depleted.kmers.combined.p.value <-
            as.numeric(df$depleted.kmers.combined.p.value)

          df$adj.p.value.estimate <-
            stats::p.adjust(df$p.value.estimate, method = p.adjust.method)
          df$adj.enriched.kmers.combined.p.value <-
            stats::p.adjust(df$enriched.kmers.combined.p.value, method = p.adjust.method)
          df$adj.depleted.kmers.combined.p.value <-
            stats::p.adjust(df$depleted.kmers.combined.p.value, method = p.adjust.method)

          motif.kmers.dfs <-
            lapply(foreground.result, function(x)
              x$motif.kmers.df)
          volcano.plots <-
            lapply(foreground.result, function(x)
              x$volcano.plot)
          perm.test.plots <-
            lapply(foreground.result, function(x)
              x$perm.test.plot)
          enriched.kmers.combined.p.values <-
            lapply(foreground.result, function(x)
              x$enriched.kmers.combined.p.value)
          depleted.kmers.combined.p.values <-
            lapply(foreground.result, function(x)
              x$depleted.kmers.combined.p.value)

          return(
            list(
              df = df,
              motif.kmers.dfs = motif.kmers.dfs,
              volcano.plots = volcano.plots,
              perm.test.plots = perm.test.plots,
              enriched.kmers.combined.p.values = enriched.kmers.combined.p.values,
              depleted.kmers.combined.p.values = depleted.kmers.combined.p.values
            )
          )
        })

      result <-
        lapply(seq_len(length(foreground.sets)), function(i) {
          enrichment.df <- enrichment.result$dfs[[i]]
          enrichment.df$kmer <- enrichment.result$kmers
          return(
            list(
              enrichment.df = enrichment.df,
              motif.df = motif.result[[i]]$df,
              motif.kmers.dfs = motif.result[[i]]$motif.kmers.dfs,
              volcano.plots = motif.result[[i]]$volcano.plots,
              perm.test.plots = motif.result[[i]]$perm.test.plots,
              enriched.kmers.combined.p.values = motif.result[[i]]$enriched.kmers.combined.p.values,
              depleted.kmers.combined.p.values = motif.result[[i]]$depleted.kmers.combined.p.values
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
#' SPMA helps to illuminate the relationship between RBP binding evidence and the transcript
#' sorting criterion, e.g., fold change between treatment and control samples.
#'
#' @inheritParams runKmerTSMA
#' @inheritParams subdivideData
#' @inheritParams scoreSpectrum
#'
#' @return A list with the following components:
#' \tabular{rl}{
#'   \code{foreground.scores} \tab the result of \code{\link{runKmerTSMA}} for the binned data\cr
#'   \code{spectrum.info.df} \tab a data frame with the SPMA results\cr
#'   \code{spectrum.plots} \tab a list of spectrum plots, as generated by \code{\link{scoreSpectrum}}\cr
#'   \code{classifier.scores} \tab a list of classifier scores, as returned by
#'   \code{\link{spectrumClassifier}}
#' }
#'
#' @details
#' In order to investigate how motif targets are distributed across a spectrum of
#' transcripts (e.g., all transcripts of a platform, ordered by fold change),
#' Spectrum Motif Analysis visualizes the gradient of RBP binding evidence
#' across all transcripts.
#'
#' The \emph{k}-mer-based approach differs from the matrix-based approach by how the sequences are
#' scored. Here, sequences are broken into \emph{k}-mers, i.e., oligonucleotide sequences of
#' \emph{k} bases.
#' And only statistically significantly enriched or depleted \emph{k}-mers are then used to
#' calculate a score for each RNA-binding protein, which quantifies its
#' target overrepresentation.
#'
#' @examples
#' \dontrun{
#' # exemplary data set
#' background.df <- ge$background
#' # sort sequences by signal-to-noise ratio
#' background.df <- dplyr::arrange(background.df, value)
#' # character vector of named sequences
#' background.set <- background.df$seq
#' names(background.set) <- paste0(background.df$refseq, "|", background.df$seq.type)
#' 
#' results <- runKmerSPMA(background.set)
#' }
#' 
#' @family SPMA functions
#' @family \emph{k}-mer functions
#' @importFrom stats p.adjust
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @export
runKmerSPMA <-
  function(background.set,
             motifs = NULL,
             k = 6,
             n.bins = 40,
             max.model.degree = 1,
             max.cs.permutations = 10000000,
             min.cs.permutations = 5000,
             fg.permutations = 5000,
             p.adjust.method = "BH",
             p.combining.method = "fisher",
             n.cores = 1) {
    # avoid CRAN note
    motif.id <-
      motif.rbps <-
      adj.r.squared <- degree <- residuals <- slope <- NULL
    f.statistic <-
      f.statistic.p.value <- f.statistic.adj.p.value <- NULL
    consistency.score <-
      consistency.score.p.value <-
      consistency.score.adj.p.value <- NULL
    consistency.score.n <- aggregate.classifier.score <- NULL
    n.significant <-
      n.very.significant <- n.extremely.significant <- NULL

    foreground.sets <- subdivideData(background.set, n.bins)

    results <- runKmerTSMA(
      foreground.sets,
      background.set,
      motifs = motifs,
      k = k,
      fg.permutations = fg.permutations,
      kmer.significance.threshold = 0.01,
      produce.plot = FALSE,
      p.adjust.method = p.adjust.method,
      p.combining.method = p.combining.method,
      n.cores = n.cores
    )

    if (length(results) > 0) {
      dfs <- lapply(results, function(result) {
        result$motif.df
      })
      enrichment.df <- do.call("rbind", dfs)
      enrichment.df$adj.p.value.estimate <-
        stats::p.adjust(enrichment.df$p.value.estimate,
          method = p.adjust.method
        )

      motifs <- getMotifs()
      spectrum.info <- lapply(motifs, function(motif) {
        motif.data.df <- dplyr::filter(enrichment.df, motif.id == motif$id)
        values <- motif.data.df$geo.mean.enrichment
        values[values == 0] <-
          0.01 # avoid -Inf after taking the log
        score <-
          scoreSpectrum(
            log(values),
            motif.data.df$adj.p.value.estimate,
            max.model.degree = max.model.degree,
            max.cs.permutations = max.cs.permutations,
            min.cs.permutations = min.cs.permutations
          )

        n.significant <-
          sum(motif.data.df$adj.p.value.estimate <= 0.05)
        classifier.score <-
          spectrumClassifier(
            score$adj.r.squared,
            score$degree,
            score$slope,
            score$consistency.score.n,
            n.significant,
            n.bins
          )

        return(
          list(
            info = list(
              motif.id = motif$id,
              motif.rbps = paste(motif$rbps, collapse = ", "),
              adj.r.squared = score$adj.r.squared,
              degree = score$degree,
              residuals = score$residuals,
              slope = score$slope,
              f.statistic = score$f.statistic,
              f.statistic.p.value = score$f.statistic.p.value,
              consistency.score = score$consistency.score,
              consistency.score.p.value = score$consistency.score.p.value,
              consistency.score.n = score$consistency.score.n,
              n.significant = n.significant,
              n.very.significant = sum(motif.data.df$adj.p.value <= 0.01),
              n.extremely.significant = sum(motif.data.df$adj.p.value <= 0.001),
              aggregate.classifier.score = sum(classifier.score)
            ),
            spectrum.plot = score$plot,
            classifier.score = classifier.score
          )
        )
      })
      spectrum.info.df <-
        as.data.frame(do.call("rbind", lapply(spectrum.info, function(x)
          x$info)),
        stringsAsFactors = FALSE
        )
      spectrum.info.df$motif.id <-
        as.character(spectrum.info.df$motif.id)
      spectrum.info.df$motif.rbps <-
        as.character(spectrum.info.df$motif.rbps)
      spectrum.info.df$adj.r.squared <-
        as.numeric(spectrum.info.df$adj.r.squared)
      spectrum.info.df$degree <-
        as.integer(spectrum.info.df$degree)
      spectrum.info.df$residuals <-
        as.numeric(spectrum.info.df$residuals)
      spectrum.info.df$slope <-
        as.numeric(spectrum.info.df$slope)
      spectrum.info.df$f.statistic <-
        as.numeric(spectrum.info.df$f.statistic)
      spectrum.info.df$f.statistic.p.value <-
        as.numeric(spectrum.info.df$f.statistic.p.value)
      spectrum.info.df$consistency.score <-
        as.numeric(spectrum.info.df$consistency.score)
      spectrum.info.df$consistency.score.p.value <-
        as.numeric(spectrum.info.df$consistency.score.p.value)
      spectrum.info.df$consistency.score.n <-
        as.integer(spectrum.info.df$consistency.score.n)
      spectrum.info.df$n.significant <-
        as.integer(spectrum.info.df$n.significant)
      spectrum.info.df$n.very.significant <-
        as.integer(spectrum.info.df$n.very.significant)
      spectrum.info.df$n.extremely.significant <-
        as.integer(spectrum.info.df$n.extremely.significant)
      spectrum.info.df$aggregate.classifier.score <-
        as.integer(spectrum.info.df$aggregate.classifier.score)

      spectrum.info.df$f.statistic.adj.p.value <-
        stats::p.adjust(spectrum.info.df$f.statistic.p.value, method = p.adjust.method)
      spectrum.info.df$consistency.score.adj.p.value <-
        stats::p.adjust(spectrum.info.df$consistency.score.p.value, method = p.adjust.method)

      spectrum.info.df <-
        dplyr::select(
          spectrum.info.df,
          motif.id,
          motif.rbps,
          adj.r.squared,
          degree,
          residuals,
          slope,
          f.statistic,
          f.statistic.p.value,
          f.statistic.adj.p.value,
          consistency.score,
          consistency.score.p.value,
          consistency.score.adj.p.value,
          consistency.score.n,
          n.significant,
          n.very.significant,
          n.extremely.significant,
          aggregate.classifier.score
        )

      spectrum.plots <-
        lapply(spectrum.info, function(x)
          x$spectrum.plot)
      classifier.scores <-
        lapply(spectrum.info, function(x)
          x$classifier.score)
    } else {
      spectrum.info.df <- data.frame(
        motif.id = character(0),
        motif.rbps = character(0),
        adj.r.squared = numeric(0),
        degree = integer(0),
        residuals = numeric(0),
        slope = numeric(0),
        f.statistic = numeric(0),
        f.statistic.p.value = numeric(0),
        f.statistic.adj.p.value = numeric(0),
        consistency.score = numeric(0),
        consistency.score.p.value = numeric(0),
        consistency.score.adj.p.value = numeric(0),
        consistency.score.n,
        n.significant = integer(0),
        n.very.significant = integer(0),
        n.extremely.significant = integer(0),
        aggregate.classifier.score = integer(0)
      )
      spectrum.plots <- NULL
      warning("no sequences in any condition")
    }

    return(
      list(
        foreground.scores = results,
        spectrum.info.df = spectrum.info.df,
        spectrum.plots = spectrum.plots,
        classifier.scores = classifier.scores
      )
    )
  }
