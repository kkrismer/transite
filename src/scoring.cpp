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
//' motif <- getMotifById("M178_0.6")[[1]]
//' sequences <- c("CAACAGCCUUAAUU", "CAGUCAAGACUCC", "CUUUGGGGAAU",
//'                "UCAUUUUAUUAAA", "AAUUGGUGUCUGGAUACUUCCCUGUACAU",
//'                "AUCAAAUUA", "UGUGGGG", "GACACUUAAAGAUCCU",
//'                "UAGCAUUAACUUAAUG", "AUGGA", "GAAGAGUGCUCA", "AUAGAC",
//'                "AGUUC", "CCAGUAA")
//' seq_char_vectors <- lapply(sequences, function(seq) {
//'   unlist(strsplit(seq, ""))
//' })
//' scoreSequences(seq_char_vectors, as.matrix(motifMatrix(motif)))
//'
//' @export
// [[Rcpp::export]]
SEXP scoreSequences(Rcpp::List sequences, Rcpp::NumericMatrix pwm) {
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

//' @title \emph{k}-mer Score Calculation
//'
//' @description
//' C++ implementation of \emph{k}-mer score calculation
//'
//' @param kmers list of \emph{k}-mers
//' @param pwm position weight matrix
//'
//' @return list of PWM scores for the specified \emph{k}-mers
//'
//' @examples
//' motif <- getMotifById("M178_0.6")[[1]]
//' kmers <- c("AAAAAA", "CAAAAA", "GAAAAA")
//' calculateKmerScores(kmers, as.matrix(motifMatrix(motif)))
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector calculateKmerScores(Rcpp::List kmers, Rcpp::NumericMatrix pwm) {
    Rcpp::NumericVector scores(kmers.size());
    for(int i(0); i < kmers.size(); ++i) {
        Rcpp::CharacterVector seq(kmers[i]);

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
                    } else {
                        sum += pwm(k, 3);
                    }
                }
                positionalScores[j] = sum;
            }
            scores[i] = max(positionalScores);
        } else {
            Rcpp::NumericVector positionalScores(pwm.nrow() - seq.size() + 1);
            for (int j(0); j < (pwm.nrow() - seq.size() + 1); ++j) {
                double sum = 0;
                for(int k(0); k < seq.size(); ++k) {
                    Rcpp::String pos = seq[k];
                    if(pos == "A") {
                        sum += pwm(j + k, 0);
                    } else if(pos == "C") {
                        sum += pwm(j + k, 1);
                    } else if(pos == "G") {
                        sum += pwm(j + k, 2);
                    } else {
                        sum += pwm(j + k, 3);
                    }
                }
                positionalScores[j] = sum;
            }
            scores[i] = max(positionalScores);
        }
    }
    return scores;
}

//' @title \emph{k}-mer Score Lookup Table Access Function
//'
//' @description
//' C++ implementation of \emph{k}-mer score hash table lookup.
//'
//' @param kmers list of \emph{k}-mers
//' @param kmerScores position weight matrix
//'
//' @return numeric vector of \emph{k}-mer scores
// [[Rcpp::export]]
Rcpp::NumericVector lookupKmerScores(Rcpp::List kmers, Rcpp::Environment kmerScores) {
    Rcpp::NumericVector scores(kmers.size());
    for(int i(0); i < kmers.size(); ++i) {
        std::string kmer = kmers[i];
        scores[i] = kmerScores[kmer];
    }
    return scores;
}

//' @title Motif Score Algorithm
//'
//' @description
//' C++ implementation of motif score algorithm.
//'
//' @param kmers list of \emph{k}-mers
//' @return data frame with columns \code{score}, \code{top_kmer},
//' and \code{top_kmer_enrichment}
// [[Rcpp::export]]
Rcpp::DataFrame computeMotifScore(Rcpp::List kmers) {
    // kmers is a list of data frames (sorted desc(score), filtered score > 0)
    Rcpp::NumericVector scores(kmers.size());
    Rcpp::CharacterVector topKmers(kmers.size());
    Rcpp::NumericVector topKmerEnrichments(kmers.size());

    for(int i(0); i < kmers.size(); ++i) {
        Rcpp::DataFrame df = Rcpp::as<Rcpp::DataFrame>(kmers[i]);

        if(df.nrows() == 0) {
            scores[i] = 0;
            topKmers[i] = "";
            topKmerEnrichments[i] = 1;
        } else {
            Rcpp::NumericVector dfScore = df["score"];
            Rcpp::CharacterVector dfKmer = df["kmer"];
            Rcpp::NumericVector dfEnrichment = df["enrichment"];
            double enriched(0.0);
            double depleted(0.0);
            for(int j(0); j < dfScore.size(); ++j) {
                // requirement: dfScore sorted in descending order
                double localScore(dfScore[j]);
                if(dfEnrichment[j] < 1) {
                    depleted += localScore / (j + 1);
                } else {
                    enriched += localScore / (j + 1);
                }
            }
            scores[i] = enriched - depleted;
            topKmers[i] = dfKmer[0];
            topKmerEnrichments[i] = dfEnrichment[0];
        }
    }
    return Rcpp::DataFrame::create(Rcpp::_["score"] = scores,
                                   Rcpp::_["top_kmer"] = topKmers,
                                   Rcpp::_["top_kmer_enrichment"] = topKmerEnrichments);
}

double calculateConsistencyScore(Rcpp::NumericVector x) {
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
//' local_consistency <- calculateLocalConsistency(poor_enrichment_spectrum,
//'   1000000, 1000, 5)
//'
//' enrichment_spectrum <- c(0.1, 0.3, 0.6, 0.7, 0.8,
//'   0.9, 1.2, 1.4, 1.6, 1.4)
//' local_consistency <- calculateLocalConsistency(enrichment_spectrum,
//'   1000000, 1000, 5)
//' @export
// [[Rcpp::export]]
Rcpp::List calculateLocalConsistency(Rcpp::NumericVector x, int numPermutations,
                               int minPermutations, int e) {
    double score(calculateConsistencyScore(x));
    int k(0);
    int i(1);
    std::mt19937 gen(std::chrono::high_resolution_clock::now().time_since_epoch().count());

    while(i <= numPermutations && (i < minPermutations || k < e)) {
        // shuffle spectrum
        Rcpp::NumericVector shuffledSpectrum = clone(x);
        std::shuffle(shuffledSpectrum.begin(), shuffledSpectrum.end(), gen);

        // calculate consistency score
        double shuffledScore(calculateConsistencyScore(shuffledSpectrum));
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
//' (returned by \code{\link{scoreTranscripts}})
//' @param totalSites number of potential binding sites per sequence
//' (returned by \code{\link{scoreTranscripts}})
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
//' motif_db <- getMotifById("M178_0.6")
//' fg <- scoreTranscripts(foreground_seqs, cache = FALSE,
//'   motifs = motif_db)
//' bg <- scoreTranscripts(background_seqs, cache = FALSE,
//'   motifs = motif_db)
//'
//' mc_result <- calculateTranscriptMC(unlist(bg$absolute_hits),
//'  unlist(bg$total_sites),
//'  fg$df$absolute_hits / fg$df$total_sites,
//'  length(foreground_seqs), 1000, 500, 5)
//' @export
// [[Rcpp::export]]
Rcpp::List calculateTranscriptMC(Rcpp::NumericVector absoluteHits, Rcpp::NumericVector totalSites,
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

