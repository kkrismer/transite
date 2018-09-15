// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <random>
#include <algorithm>
#include <chrono>
#include <unordered_set>

using namespace Rcpp;

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
//' sequences <- c("CAACAGCCTTAATT", "CAGTCAAGACTCC", "CTTTGGGGAAT", "TCATTTTATTAAA",
//'   "AATTGGTGTCTGGATACTTCCCTGTACAT", "ATCAAATTA", "TGTGGGG", "GACACTTAAAGATCCT",
//'   "TAGCATTAACTTAATG", "ATGGA", "GAAGAGTGCTCA", "ATAGAC", "AGTTC", "CCAGTAA")
//' seq.char.vectors <- lapply(sequences, function(seq) {
//'   unlist(strsplit(seq, ""))
//' })
//' scoreSequences(seq.char.vectors, as.matrix(motif$matrix))
//'
//' @export
// [[Rcpp::export]]
SEXP scoreSequences(List sequences, NumericMatrix pwm) {
  std::vector<NumericVector> scores;
  for(int i(0); i < sequences.size(); ++i) {
    CharacterVector seq(sequences[i]);

    if(seq.size() >= pwm.nrow()) {
      NumericVector positionalScores(seq.size() - pwm.nrow() + 1);
      for (int j(0); j < (seq.size() - pwm.nrow() + 1); ++j) {
        double sum = 0;
        for(int k(0); k < pwm.nrow(); ++k) {
          String pos = seq[j + k];
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
      scores.push_back(positionalScores);
    } else {
      NumericVector positionalScores(0);
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
//' calculateKmerScores(kmers, as.matrix(motif$matrix))
//'
//' @export
// [[Rcpp::export]]
NumericVector calculateKmerScores(List kmers, NumericMatrix pwm) {
  NumericVector scores(kmers.size());
  for(int i(0); i < kmers.size(); ++i) {
    CharacterVector seq(kmers[i]);

    if(seq.size() >= pwm.nrow()) {
      NumericVector positionalScores(seq.size() - pwm.nrow() + 1);
      for (int j(0); j < (seq.size() - pwm.nrow() + 1); ++j) {
        double sum = 0;
        for(int k(0); k < pwm.nrow(); ++k) {
          String pos = seq[j + k];
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
      NumericVector positionalScores(pwm.nrow() - seq.size() + 1);
      for (int j(0); j < (pwm.nrow() - seq.size() + 1); ++j) {
        double sum = 0;
        for(int k(0); k < seq.size(); ++k) {
          String pos = seq[k];
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
NumericVector lookupKmerScores(List kmers, Environment kmerScores) {
  NumericVector scores(kmers.size());
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
//' @return data frame with columns \code{score}, \code{top.kmer},
//' and \code{top.kmer.enrichment}
// [[Rcpp::export]]
DataFrame computeMotifScore(List kmers) {
// kmers is a list of data frames (sorted desc(score), filtered score > 0)
  NumericVector scores(kmers.size());
  CharacterVector topKmers(kmers.size());
  NumericVector topKmerEnrichments(kmers.size());

  for(int i(0); i < kmers.size(); ++i) {
    DataFrame df = kmers[i];

    if(df.nrows() == 0) {
      scores[i] = 0;
      topKmers[i] = "";
      topKmerEnrichments[i] = 1;
    } else {
      NumericVector dfScore = df["score"];
      CharacterVector dfKmer = df["kmer"];
      NumericVector dfEnrichment = df["enrichment"];
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
  return DataFrame::create(_["score"]= scores, _["top.kmer"]= topKmers, _["top.kmer.enrichment"]= topKmerEnrichments);
}

double calculateConsistencyScore(NumericVector x) {
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
//' @param numPermutations maximum number of permutations performed in Monte Carlo test
//' for consistency score
//' @param minPermutations minimum number of permutations performed in Monte Carlo test
//' for consistency score
//' @param e stop criterion for consistency score Monte Carlo test: aborting permutation
//' process after observing \code{e} random consistency values with more extreme values
//' than the actual consistency value
//' @return list with \code{score}, \code{p.value}, and \code{n} components
//'
//' @examples
//' poor.enrichment.spectrum <- c(0.1, 0.5, 0.6, 0.4, 0.7, 0.6, 1.2, 1.1, 1.8, 1.6)
//' local.consistency <- calculateLocalConsistency(enrichment.values, 1000000, 1000, 5)
//'
//' enrichment.spectrum <- c(0.1, 0.3, 0.6, 0.7, 0.8, 0.9, 1.2, 1.4, 1.6, 1.4)
//' local.consistency <- calculateLocalConsistency(enrichment.values, 1000000, 1000, 5)
//' @export
// [[Rcpp::export]]
List calculateLocalConsistency(NumericVector x, int numPermutations, int minPermutations, int e) {
  double score(calculateConsistencyScore(x));
  int k(0);
  int i(1);
  std::mt19937 gen(std::chrono::high_resolution_clock::now().time_since_epoch().count());

  while(i <= numPermutations && (i < minPermutations || k < e)) {
    // shuffle spectrum
    NumericVector shuffledSpectrum = clone(x);
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

  return List::create(Named("score") = score, Named("p.value") = pValue, Named("n") = i);
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
//' @param maxPermutations maximum number of foreground permutations performed in
//' Monte Carlo test for enrichment score
//' @param minPermutations minimum number of foreground permutations performed in
//' Monte Carlo test for enrichment score
//' @param e stop criterion for enrichment score Monte Carlo test: aborting permutation process
//' after observing \code{e} random enrichment values with more extreme values than the actual
//' enrichment value
//'
//' @return list with p-value and number of iterations of Monte Carlo sampling
//' for local consistency score
//'
//' @examples
//' foreground.seqs <- c("CAGTCAAGACTCC", "AATTGGTGTCTGGATACTTCCCTGTACAT", "AGAT", "CCAGTAA")
//' background.seqs <- c("CAACAGCCTTAATT", "CAGTCAAGACTCC", "CTTTGGGGAAT",
//'                      "TCATTTTATTAAA", "AATTGGTGTCTGGATACTTCCCTGTACAT",
//'                      "ATCAAATTA", "AGAT", "GACACTTAAAGATCCT",
//'                      "TAGCATTAACTTAATG", "ATGGA", "GAAGAGTGCTCA",
//'                      "ATAGAC", "AGTTC", "CCAGTAA")
//' foreground.scores <- scoreTranscripts(foreground.seqs, cache = FALSE)
//' background.scores <- scoreTranscripts(background.seqs, cache = FALSE)
//'
//' fg <- dplyr::filter(foreground.scores$df, motif.id == "M178_0.6")
//' bg <- dplyr::filter(background.scores$df, motif.id == "M178_0.6")
//'
//' mc.result <- calculateTranscriptMC(bg$absolute.hits, bg$total.sites,
//'                                    fg$absolute.hits / fg$total.sites,
//'                                    length(foreground.seqs), 10000, 5000, 5)
//'
// [[Rcpp::export]]
List calculateTranscriptMC(NumericVector absoluteHits, NumericVector totalSites, double relHitsForeground,
                             int n, int maxPermutations, int minPermutations, int e) {
  double relHitsBackground(sum(absoluteHits) / sum(totalSites));
  double actualScore(std::abs(relHitsForeground - relHitsBackground));
  int k(0);

  // random-number engine used (Mersenne-Twister) and init seed (alternatively: std::mt19937 gen(std::time(0));)
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

  return List::create(Named("p.value") = pValue, Named("n") = i);
}

