// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <random>
#include <algorithm>
#include <chrono>
#include <unordered_set>

//' @title Score Sequences with PWM
//'
//' @description
//' C++ implementation of PWM scoring algorithm
//'
//' @param sequences list of sequences
//' @param pwm position weight matrix
//'
//' @return list of PWM scores for each sequence
//' @examples
//' motif <- get_motif_by_id("M178_0.6")[[1]]
//' sequences <- c("CAACAGCCUUAAUU", "CAGUCAAGACUCC", "CUUUGGGGAAU",
//'                "UCAUUUUAUUAAA", "AAUUGGUGUCUGGAUACUUCCCUGUACAU",
//'                "AUCAAAUUA", "UGUGGGG", "GACACUUAAAGAUCCU",
//'                "UAGCAUUAACUUAAUG", "AUGGA", "GAAGAGUGCUCA", "AUAGAC",
//'                "AGUUC", "CCAGUAA")
//' seq_char_vectors <- lapply(sequences, function(seq) {
//'   unlist(strsplit(seq, ""))
//' })
//' score_sequences(seq_char_vectors, as.matrix(get_motif_matrix(motif)))
//'
//' @export
// [[Rcpp::export]]
SEXP score_sequences(Rcpp::List sequences, Rcpp::NumericMatrix pwm) {
    std::vector<Rcpp::NumericVector> scores;
    for(int i(0); i < sequences.size(); ++i) {
        Rcpp::CharacterVector seq(sequences[i]);

        if(seq.size() >= pwm.nrow()) {
            Rcpp::NumericVector positionalScores(seq.size() - pwm.nrow() + 1);
            for (int j(0); j < (seq.size() - pwm.nrow() + 1); ++j) {
                double sum = 0;
                for(int k(0); k < pwm.nrow(); ++k) {
                    Rcpp::String pos = seq[j + k];
                    if(pos == "A") {
                        sum += pwm(k, 0);
                    } else if(pos == "C") {
                        sum += pwm(k, 1);
                    } else if(pos == "G") {
                        sum += pwm(k, 2);
                    } else if(pos == "U") {
                        sum += pwm(k, 3);
                    } else if(pos == "T") {
                        sum += pwm(k, 3);
                    } else {
                        throw std::invalid_argument(std::string("invalid character in RNA sequence: ") + std::string(pos) + std::string(" (valid characters: A, C, G, U)"));
                    }
                }
                positionalScores[j] = sum;
            }
            scores.push_back(positionalScores);
        } else {
            Rcpp::NumericVector positionalScores(0);
            scores.push_back(positionalScores);
        }
    }
    return wrap(scores);
}

double calculate_consistency_score(Rcpp::NumericVector x) {
    double score(0.0);
    for(int i(0); i < x.size() - 2; ++i) {
        score += std::abs(((x[i] + x[i + 2]) / 2) - x[i + 1]);
    }
    return score / ((double)x.size());
}


//' @title Local Consistency Score
//'
//' @description
//' C++ implementation of Local Consistency Score algorithm.
//'
//' @param x numeric vector that contains values for shuffling
//' @param numPermutations maximum number of permutations performed in
//' Monte Carlo test
//' for consistency score
//' @param minPermutations minimum number of permutations performed in
//' Monte Carlo test
//' for consistency score
//' @param e stop criterion for consistency score Monte Carlo test:
//' aborting permutation
//' process after observing \code{e} random consistency values with
//' more extreme values
//' than the actual consistency value
//' @return list with \code{score}, \code{p_value}, and \code{n} components,
//' where \code{score} is the raw local consistency score (usually not used),
//' \code{p_value} is the associated p-value for that score, obtained by
//' Monte Carlo testing, and \code{n} is the number of permutations performed
//' in the Monte Carlo test (the higher, the more significant)
//'
//' @examples
//' poor_enrichment_spectrum <- c(0.1, 0.5, 0.6, 0.4,
//'   0.7, 0.6, 1.2, 1.1, 1.8, 1.6)
//' local_consistency <- calculate_local_consistency(poor_enrichment_spectrum,
//'   1000000, 1000, 5)
//'
//' enrichment_spectrum <- c(0.1, 0.3, 0.6, 0.7, 0.8,
//'   0.9, 1.2, 1.4, 1.6, 1.4)
//' local_consistency <- calculate_local_consistency(enrichment_spectrum,
//'   1000000, 1000, 5)
//' @export
// [[Rcpp::export]]
Rcpp::List calculate_local_consistency(Rcpp::NumericVector x, int numPermutations,
                               int minPermutations, int e) {
    double score(calculate_consistency_score(x));
    int k(0);
    int i(1);
    std::mt19937 gen(std::chrono::high_resolution_clock::now().time_since_epoch().count());

    while(i <= numPermutations && (i < minPermutations || k < e)) {
        // shuffle spectrum
        Rcpp::NumericVector shuffledSpectrum = clone(x);
        std::shuffle(shuffledSpectrum.begin(), shuffledSpectrum.end(), gen);

        // calculate consistency score
        double shuffledScore(calculate_consistency_score(shuffledSpectrum));
        // lower tail probability
        if(shuffledScore <= score) {
            ++k;
        }

        ++i;
    }

    double pValue(((double)k + 1.0) / ((double)i + 1.0));

    return Rcpp::List::create(Rcpp::Named("score") = score,
                              Rcpp::Named("p_value") = pValue,
                              Rcpp::Named("n") = i);
}

//' @title Motif Enrichment calculation
//'
//' @description
//' C++ implementation of Motif Enrichment calculation
//'
//' @param absoluteHits number of putative binding sites per sequence
//' (returned by \code{\link{score_transcripts}})
//' @param totalSites number of potential binding sites per sequence
//' (returned by \code{\link{score_transcripts}})
//' @param relHitsForeground relative number of hits in foreground set
//' @param n number of sequences in the foreground set
//' @param maxPermutations maximum number of foreground permutations
//' performed in
//' Monte Carlo test for enrichment score
//' @param minPermutations minimum number of foreground permutations
//' performed in
//' Monte Carlo test for enrichment score
//' @param e stop criterion for enrichment score Monte Carlo test:
//' aborting permutation process
//' after observing \code{e} random enrichment values with more extreme
//' values than the actual
//' enrichment value
//'
//' @return list with p-value and number of iterations of Monte Carlo sampling
//' for foreground enrichment
//'
//' @examples
//' foreground_seqs <- c("CAGUCAAGACUCC", "AAUUGGUUGUGGGGCUUCCCUGUACAU",
//'                      "AGAU", "CCAGUAA", "UGUGGGG")
//' background_seqs <- c(foreground_seqs, "CAACAGCCUUAAUU", "CUUUGGGGAAU",
//'                      "UCAUUUUAUUAAA", "AUCAAAUUA", "GACACUUAAAGAUCCU",
//'                      "UAGCAUUAACUUAAUG", "AUGGA", "GAAGAGUGCUCA",
//'                      "AUAGAC", "AGUUC")
//' motif_db <- get_motif_by_id("M178_0.6")
//' fg <- score_transcripts(foreground_seqs, cache = FALSE,
//'   motifs = motif_db)
//' bg <- score_transcripts(background_seqs, cache = FALSE,
//'   motifs = motif_db)
//'
//' mc_result <- calculate_transcript_mc(unlist(bg$absolute_hits),
//'  unlist(bg$total_sites),
//'  fg$df$absolute_hits / fg$df$total_sites,
//'  length(foreground_seqs), 1000, 500, 5)
//' @export
// [[Rcpp::export]]
Rcpp::List calculate_transcript_mc(Rcpp::NumericVector absoluteHits, Rcpp::NumericVector totalSites,
                           double relHitsForeground,
                           int n, int maxPermutations, int minPermutations, int e) {
    double relHitsBackground(sum(absoluteHits) / sum(totalSites));
    double actualScore(std::abs(relHitsForeground - relHitsBackground));
    int k(0);

    // random-number engine used (Mersenne-Twister) and init seed
    std::mt19937 gen(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    std::uniform_int_distribution<unsigned int> dist(0, absoluteHits.length() - 1);

    int i(1);
    unsigned int idx(0);
    while(i <= maxPermutations && (i < minPermutations || k < e)) {
        // select n transcripts randomly
        int randomAbsoluteHits(0);
        int randomTotalSites(0);
        std::unordered_set<unsigned int> indices;

        for(int j(0); j < n; ++j) {
            // sampling with replacement, i.e., bootstrap
            // idx = dist(gen);

            // sampling without replacement, i.e., permutation test
            // keep sampling until new sample not already in set
            do {
                idx = dist(gen);
            } while (indices.find(idx) != indices.end());
            indices.insert(idx);

            randomAbsoluteHits += absoluteHits[idx];
            randomTotalSites += totalSites[idx];
        }

        // calculate random score
        double randomScore(((double)randomAbsoluteHits) / ((double)randomTotalSites) - relHitsBackground);
        // two-tailed probability
        if(std::abs(randomScore) >= actualScore) {
            ++k;
        }

        ++i;
    }

    double pValue(((double)k + 1.0) / ((double)i + 1.0));

    return Rcpp::List::create(Rcpp::Named("p_value") = pValue,
                              Rcpp::Named("n") = i);
}
