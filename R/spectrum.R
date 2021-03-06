#' @title Subdivides Sequences into \emph{n} Bins
#'
#' @description
#' Preprocessing function for SPMA, divides transcript sequences
#' into \emph{n} bins.
#'
#' @param sorted_transcript_sequences character vector of named sequences (names are
#' usually RefSeq
#' identifiers and sequence region labels,
#'  e.g., "NM_1_DUMMY|3UTR"). It is important that the sequences are
#'  already sorted by fold change,
#'  signal-to-noise ratio or any other meaningful measure.
#' @param n_bins specifies the number of bins in which the sequences
#' will be divided,
#' valid values are between 7 and 100
#' @return An array of \code{n_bins} length, containing the binned sequences
#' @examples
#' # toy example
#' toy_seqs <- c(
#'   "CAACAGCCUUAAUU", "CAGUCAAGACUCC", "CUUUGGGGAAU", "UCAUUUUAUUAAA",
#'   "AAUUGGUGUCUGGAUACUUCCCUGUACAU", "AUCAAAUUA", "AGAU", "GACACUUAAAGAUCCU",
#'   "UAGCAUUAACUUAAUG", "AUGGA", "GAAGAGUGCUCA", "AUAGAC", "AGUUC", "CCAGUAA"
#' )
#' # names are used as keys in the hash table (cached version only)
#' # ideally sequence identifiers (e.g., RefSeq ids) and
#' # sequence region labels (e.g., 3UTR for 3'-UTR)
#' names(toy_seqs) <- c(
#'   "NM_1_DUMMY|3UTR", "NM_2_DUMMY|3UTR", "NM_3_DUMMY|3UTR",
#'   "NM_4_DUMMY|3UTR", "NM_5_DUMMY|3UTR", "NM_6_DUMMY|3UTR",
#'   "NM_7_DUMMY|3UTR",
#'   "NM_8_DUMMY|3UTR", "NM_9_DUMMY|3UTR", "NM_10_DUMMY|3UTR",
#'   "NM_11_DUMMY|3UTR",
#'   "NM_12_DUMMY|3UTR", "NM_13_DUMMY|3UTR", "NM_14_DUMMY|3UTR"
#' )
#'
#' foreground_sets <- subdivide_data(toy_seqs, n_bins = 7)
#'
#' # example data set
#' background_df <- transite:::ge$background_df
#' # sort sequences by signal-to-noise ratio
#' background_df <- dplyr::arrange(background_df, value)
#' # character vector of named sequences
#' background_seqs <- background_df$seq
#' names(background_seqs) <- paste0(background_df$refseq, "|",
#'   background_df$seq_type)
#'
#' foreground_sets <- subdivide_data(background_seqs)
#' @family SPMA functions
#' @export
subdivide_data <- function(sorted_transcript_sequences, n_bins = 40) {
    if (n_bins < 7 || n_bins > 100) {
        stop("value of bin_num is invalid, valid values are integers between 7 and 100")
    }
    if (n_bins > length(sorted_transcript_sequences)) {
        stop("more bins than transcripts")
    }

    bin_size <- floor(length(sorted_transcript_sequences) / n_bins)
    cuts <- seq_len(n_bins - 1) * bin_size
    cuts <- c(0, cuts, length(sorted_transcript_sequences))
    foreground_sets <- tapply(sorted_transcript_sequences,
                              cut(seq_len(length(sorted_transcript_sequences)),
                                  breaks = cuts), identity)
    return(foreground_sets)
}

#' @title Calculates spectrum scores and creates spectrum plots
#'
#' @description
#' Spectrum scores are a means to evaluate if a spectrum has a meaningful
#' (i.e., biologically relevant) or a random pattern.
#'
#' @details
#' One way to quantify the meaningfulness of a spectrum is to calculate the
#' deviance
#' between the linear interpolation of the scores of two adjoining bins and
#' the score
#' of the middle bin, for each position in the spectrum. The lower the score,
#' the more consistent
#' the trend in the spectrum plot. Formally, the local consistency score
#' \eqn{x_c} is defined as
#' \deqn{x_c = \frac{1}{n} \sum_{i = 1}^{n - 2}{|\frac{s_i + s_{i + 2}}{2} - s_{i + 1}|}.}
#' In order to obtain an estimate of the significance of a particular
#' score \eqn{x_c'},
#' Monte Carlo sampling
#' is performed by randomly permuting the coordinates of the scores
#' vector \eqn{s}
#' and recomputing \eqn{x_c}.
#' The probability estimate \eqn{\hat{p}} is given by the lower tail
#' version of the cumulative
#' distribution function
#' \deqn{\hat{Pr}(T(x)) = \frac{\sum_{i = 1}^n 1(T(y_i) \le T(x)) + 1}{n + 1},}
#' where \eqn{1} is the indicator function, \eqn{n} is the
#' sample size, i.e.,
#' the number of performed permutations, and \eqn{T} equals \eqn{x_c}
#' in the above equation.
#'
#' An alternative approach to assess the consistency of a spectrum plot is
#' via polynomial regression.
#' In a first step, polynomial regression models of various degrees are
#' fitted to the data, i.e.,
#' the dependent variable \eqn{s} (vector of scores), and
#' orthogonal
#' polynomials of the independent variable \eqn{b} (vector of
#' bin numbers).
#' Secondly, the model that reflects best the true nature of the data is
#' selected by
#' means of the F-test. And lastly, the adjusted \eqn{R^2} and the sum of
#' squared residuals
#' are calculated to indicate how well the model fits the data. These
#' statistics are used as
#' scores to rank the spectrum plots.
#' In general, the polynomial regression equation is
#' \deqn{y_i = \beta_0 + \beta_1 x_i + \beta_2 x_i^2 + \cdots + \beta_m x_i^m + \epsilon_i,}
#' where \eqn{m} is the degree of the polynomial (usually \eqn{m \le 5}),
#' and \eqn{\epsilon_i}
#' is the error term.
#' The dependent variable \eqn{y} is the vector of
#' scores \eqn{s}
#' and \eqn{x}
#' to \eqn{x^m} are the orthogonal polynomials of
#' the vector of bin
#' numbers \eqn{b}.
#' Orthogonal polynomials are used in order to reduce the correlation
#' between the different
#' powers of \eqn{b}
#' and therefore avoid multicollinearity in the model. This is important,
#' because correlated
#' predictors lead to unstable coefficients, i.e., the coefficients of a
#' polynomial
#' regression model of
#' degree \eqn{m} can be greatly different from a model of degree \eqn{m + 1}.
#'
#' The orthogonal polynomials of vector \eqn{b} are obtained by
#' centering (subtracting the mean),
#' QR decomposition, and subsequent normalization.
#' Given the dependent variable \eqn{y} and the orthogonal
#' polynomials
#' of \eqn{b} \eqn{x} to
#' \eqn{x^m}, the model coefficients \eqn{\beta}
#' are chosen
#' in a way to minimize
#' the deviance between the actual and the predicted values characterized by
#' \deqn{M(x) = \beta_0 + \beta_1 x + \beta_2 x^2 + \cdots + \beta_m x^m}
#' \deqn{M = argmin_{M}{(\sum_{i = 1}^n{L(y_i, M(x_i))})},}
#' where L(actual value, predicted value) denotes
#' the loss function.
#'
#' Ordinary least squares is used as estimation method for the model
#' coefficients \eqn{\beta}. The loss
#' function of ordinary least squares is the sum of squared residuals
#' (SSR)
#' and is defined as follows
#' SSR\eqn{(y,
#' \hat{y}) = \sum_{i = 1}^n{(y_i - \hat{y}_i)^2},}
#' where \eqn{y} are the observed data
#' and \eqn{\hat{y}} the
#' model predictions.
#'
#' Thus the ordinary least squares estimate of the
#' coefficients \eqn{\hat{\beta}}
#' (including the
#' intercept \eqn{\hat{\beta}_0}) of the model \eqn{M} is defined by
#' \deqn{\hat{\beta} = argmin_{\beta}{(\sum_{i = 1}^n{(y_i - \beta_0 - \sum_{j = 1}^m{\beta_j x_i^j})^2})}.}
#'
#' After polynomial models of various degrees have been fitted to the
#' data, the F-test is used to select
#' the model that best fits the data. Since the SSR
#' monotonically decreases with
#' increasing model degree (model complexity), the relative decrease
#' of the SSR between the
#' simpler model and the more complex model must outweigh the increase
#' in model complexity between the two
#' models. The F-test gives the probability that a relative decrease of
#' the SSR between the simpler
#' and the more complex model given their respective degrees of freedom is due
#' to chance. A low p-value
#' indicates that the additional degrees of freedom of the more complex model
#' lead to a better fit of
#' the data than would be expected after a mere increase of degrees of freedom.
#'
#' The F-statistic is calculated as follows
#' \deqn{F = \frac{(SSR_1 - SSR_2) / (p_2 - p_1)}{SSR_2 / (n - p_2)},}
#' where \eqn{SSR_i} is the sum of squared residuals
#' and \eqn{p_i} is the number
#' of parameters of
#' model \eqn{i}. The number of data points, i.e., bins, is denoted as \eqn{n}.
#' \eqn{F} is distributed according to the F-distribution
#' with \eqn{df_1 = p_2 - p_1} and \eqn{df_2 = n - p_2}.
#'
#' @param x vector of values (e.g., enrichment values, normalized RBP
#' scores) per bin
#' @param p_values vector of p-values (e.g., significance of enrichment
#' values) per bin
#' @param x_label label of values (e.g., \code{"enrichment value"})
#' @param sorted_transcript_values vector of sorted transcript values, i.e.,
#' the fold change or signal-to-noise ratio or any other quantity that was used
#' to sort the transcripts that were passed to \code{run_matrix_spma} or
#' \code{run_kmer_spma} (default value is \code{NULL}). These values are
#' displayed as a semi-transparent area over the enrichment value heatmaps
#' of spectrum plots.
#' @param transcript_values_label label of transcript sorting criterion
#' (e.g., \code{"log fold change"}, default value is \code{"transcript value"}),
#' only shown if \code{!is.null(sorted_transcript_values)}
#' @param midpoint for enrichment values the midpoint should be \code{1},
#' for log enrichment values \code{0} (defaults to \code{0})
#' @param x_value_limits sets limits of the x-value color scale (used to
#' harmonize color scales of different spectrum plots), see \code{limits}
#' argument of \code{\link[ggplot2]{continuous_scale}} (defaults to
#' \code{NULL}, i.e., the data-dependent default scale range)
#' @param max_model_degree maximum degree of polynomial
#' @param max_cs_permutations maximum number of permutations performed in
#' Monte Carlo test for consistency score
#' @param min_cs_permutations minimum number of permutations performed in
#' Monte Carlo test for consistency score
#' @param e integer-valued stop criterion for consistency score Monte
#' Carlo test: aborting permutation process after
#' observing \code{e} random consistency values with more extreme values than
#' the actual consistency value
#' @return A list object of class \code{SpectrumScore} with the following
#' components:
#' \tabular{rl}{
#'   \code{adj_r_squared} \tab adjusted \eqn{R^2} of polynomial model\cr
#'   \code{degree} \tab maximum degree of polynomial\cr
#'   \code{residuals} \tab residuals  of polynomial model\cr
#'   \code{slope} \tab coefficient of the linear term of the polynomial model
#'    (spectrum "direction")\cr
#'   \code{f_statistic} \tab statistic of the F-test \cr
#'   \code{f_statistic_p_value} \tab p-value of F-test\cr
#'   \code{consistency_score} \tab normalized sum of deviance between the linear
#'   interpolation of the scores
#'    of two adjoining bins and the score of the middle bin, for each position
#'    in the spectrum \cr
#'   \code{consistency_score_p_value} \tab obtained by Monte Carlo sampling
#'   (randomly permuting
#'   the coordinates of the scores vector) \cr
#'   \code{consistency_score_n} \tab number of permutations\cr
#'   \code{plot} \tab
#' }
#' @examples
#' # random spectrum
#' score_spectrum(runif(n = 40, min = -1, max = 1), max_model_degree = 1)
#'
#' # two random spectrums with harmonized color scales
#' plot(score_spectrum(runif(n = 40, min = -1, max = 1), max_model_degree = 1,
#'      x_value_limits = c(-2.0, 2.0)))
#' plot(score_spectrum(runif(n = 40, min = -2, max = 2), max_model_degree = 1,
#'      x_value_limits = c(-2.0, 2.0)))
#'
#' # random spectrum with p-values
#' score_spectrum(runif(n = 40, min = -1, max = 1),
#'                p_values = runif(n = 40, min = 0, max = 1),
#'                max_model_degree = 1)
#'
#' # random spectrum with sorted transcript values
#' log_fold_change <- log(runif(n = 1000, min = 0, max = 1) /
#'                            runif(n = 1000, min = 0, max = 1))
#' score_spectrum(runif(n = 40, min = -1, max = 1),
#'                sorted_transcript_values = sort(log_fold_change),
#'                max_model_degree = 1)
#'
#' # non-random linear spectrum
#' signal <- seq(-1, 0.99, 2 / 40)
#' noise <- rnorm(n = 40, mean = 0, sd = 0.5)
#' score_spectrum(signal + noise, max_model_degree = 1,
#'                max_cs_permutations = 100000)
#'
#' # non-random quadratic spectrum
#' signal <- seq(-1, 0.99, 2 / 40)^2 - 0.5
#' noise <- rnorm(n = 40, mean = 0, sd = 0.2)
#' score_spectrum(signal + noise, max_model_degree = 2,
#'                max_cs_permutations = 100000)
#' @family SPMA functions
#' @importFrom dplyr arrange
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_rect
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_area
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom scales pretty_breaks
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 labs
#' @importFrom stats lm
#' @importFrom stats poly
#' @importFrom stats sd
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 geom_smooth
#' @importFrom stats predict
#' @importFrom stats pf
#' @importFrom grDevices pdf
#' @importFrom ggplot2 ggplot_gtable
#' @importFrom ggplot2 ggplot_build
#' @importFrom grDevices dev.off
#' @importFrom gridExtra arrangeGrob
#' @export
score_spectrum <- function(x, p_values = array(1, length(x)),
                           x_label = "log enrichment",
                           sorted_transcript_values = NULL,
                           transcript_values_label = "transcript value",
                           midpoint = 0,
                           x_value_limits = NULL,
                           max_model_degree = 3,
                           max_cs_permutations = 10000000,
                           min_cs_permutations = 5000, e = 5) {
    transcript_value <- position <- bin <- bin_xmin <- bin_xmax <-
        value <- p_value <- NULL

    num_bins <- length(x)

    if (num_bins < 7) {
        stop("too few bins")
    }
    if (max_model_degree > 5 || max_model_degree < 1) {
        stop("supported values for max_model_degree: 1, 2, 3, 4, 5")
    }

    bins <- seq_len(num_bins)

    if (is.null(sorted_transcript_values)) {
        sorted_transcript_values <- 1
    }
    enrichment_values_df <- data.frame(bin = bins,
                                       value = x,
                                       p_value = p_values,
                                       bin_xmin = bins - 0.5,
                                       bin_xmax = bins + 0.5,
                                       stringsAsFactors = FALSE)
    transcript_values_df <- data.frame(transcript_value = sorted_transcript_values,
                                       stringsAsFactors = FALSE)
    num_transcripts <- nrow(transcript_values_df)
    transcript_values_df$position <- seq_len(num_transcripts) / (num_transcripts / num_bins) + 0.5

    if (num_transcripts > 1) {
        value_range <- range(sorted_transcript_values)
        asterisk_y_pos <- (value_range[2] - value_range[1]) * 0.9 - abs(value_range[1])
    } else {
        asterisk_y_pos <- 0.9
    }

    enrichment_plot <- ggplot2::ggplot() +
        ggplot2::geom_rect(data = enrichment_values_df,
                           ggplot2::aes(xmin = bin_xmin, xmax = bin_xmax,
                                        ymin = -Inf, ymax = Inf,
                                        fill = value)) +
        ggplot2::geom_area(data = transcript_values_df,
                           ggplot2::aes(x = position,
                                        y = transcript_value),
                           fill = "black", alpha = 0.2) +
        ggplot2::scale_x_continuous(limits = c(0.5, num_bins + 0.5),
                                    breaks = bins,
                                    expand = c(0, 0)) +
        ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
        ggplot2::scale_fill_gradient2(
            low = "#004F91FF", mid = "white", high = "#832424FF",
            midpoint = midpoint, limits = x_value_limits
        )  +
        ggplot2::geom_text(data = enrichment_values_df,
                           ggplot2::aes(x = bin,
                                        label = ifelse(p_value <= 0.001, "***",
                                                       ifelse(p_value <= 0.01, "**",
                                                              ifelse(p_value <= 0.05, "*", "")))),
                           size = 4, y = asterisk_y_pos
        ) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(size = 8),
            axis.title.x = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_text(face = "bold"),
            panel.background = ggplot2::element_blank(),
            legend.position = "top",
            legend.title = ggplot2::element_text(face = "bold")
        ) +
        ggplot2::labs(fill = x_label,
                      y = ifelse(num_transcripts == 1,
                                 "", transcript_values_label))

    if (num_transcripts == 1) {
        enrichment_plot <- enrichment_plot + ggplot2::theme(
            axis.ticks.y = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank()
        )
    }

    # determine polynomial degree of gradient
    model <- stats::lm(x ~ stats::poly(bins, degree = max_model_degree))
    coeffs <- summary(model)$coefficients

    if (stats::sd(x) <= 1e-10) {
        degree <- 0
    } else {
        if (max_model_degree == 5 && coeffs[6, 4] <= 0.001) {
            degree <- 5
        } else if (max_model_degree >= 4 && coeffs[5, 4] <= 0.001) {
            degree <- 4
        } else if (max_model_degree >= 3 && coeffs[4, 4] <= 0.001) {
            degree <- 3
        } else if (max_model_degree >= 2 && coeffs[3, 4] <= 0.001) {
            degree <- 2
        } else if (coeffs[2, 4] <= 0.001) {
            degree <- 1
        } else {
            degree <- 0
        }
    }

    df <- data.frame(y = x, bin = bins)
    scatterplot <- ggplot2::ggplot(df,
                                   ggplot2::aes_string(x = "bin", y = "y")) +
        ggplot2::geom_point() +
        ggplot2::ylab(x_label) +
        ggplot2::scale_x_continuous(limits = c(0.5, num_bins + 0.5),
                                    breaks = bins,
                                    expand = c(0, 0)) +
        ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(size = 8),
            axis.title.x = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_text(face = "bold"),
            panel.background = ggplot2::element_blank(),
            strip.background = ggplot2::element_blank()
        )
    if (degree > 0) {
        selected_model <- stats::lm(x ~ stats::poly(bins, degree = degree))
        scatterplot <- scatterplot +
            ggplot2::geom_smooth(method = "lm",
                                 formula = y ~ stats::poly(x, degree = degree))
    } else {
        selected_model <- stats::lm(x ~ 1)
        scatterplot <- scatterplot +
            ggplot2::geom_smooth(method = "lm", formula = y ~ 1)
    }

    # calculate residuals of selected model
    predicted <- stats::predict(selected_model, data.frame(bins))
    residuals <- sum((predicted - x)^2)
    meta_info <- summary(selected_model)

    # calculate local consistency score
    local_consistency <- calculate_local_consistency(x, max_cs_permutations,
                                                     min_cs_permutations, e)

    if (degree > 0) {
        slope <- meta_info$coefficients[2, 1]
        f_statistic <- meta_info$fstatistic
        f_statistic_p_value <- stats::pf(f_statistic[[1]],
                                         df1 = f_statistic[[2]],
                                         df2 = f_statistic[[3]],
                                         lower.tail = FALSE
        )
    } else {
        slope <- 0
        f_statistic <- list(0, NA, NA)
        f_statistic_p_value <- 1
    }

    grDevices::pdf(NULL)
    gp1 <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(enrichment_plot))
    gp2 <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(scatterplot))
    invisible(grDevices::dev.off())

    gp1$widths <- gp2$widths

    return(.SpectrumScore(
        adj_r_squared = meta_info$adj.r.squared,
        degree = as.integer(degree),
        residuals = residuals,
        slope = slope,
        f_statistic = f_statistic[[1]],
        f_statistic_p_value = f_statistic_p_value,
        consistency_score = local_consistency$score,
        consistency_score_p_value = local_consistency$p_value,
        consistency_score_n = as.integer(local_consistency$n),
        plot = gridExtra::arrangeGrob(gp1, gp2, ncol = 1, heights = c(1, 1.4))
    ))
}

#' @title Simple spectrum classifier based on empirical thresholds
#'
#' @description
#' Spectra can be classified based on the aggregate spectrum classifier score.
#' If \code{sum(score) == 3} spectrum considered non-random, random otherwise.
#'
#' @param adj_r_squared adjusted \eqn{R^2} of polynomial model, returned by
#' \link{score_spectrum}
#' @param degree degree of polynomial, returned by \link{score_spectrum}
#' @param slope coefficient of the linear term of the polynomial model
#' (spectrum "direction"),
#' returned by \link{score_spectrum}
#' @param consistency_score_n number of performed permutations before
#' early stopping,
#' returned by \link{score_spectrum}
#' @param n_significant number of bins with statistically significant
#' enrichment
#' @param n_bins number of bins
#' @return a three-dimensional binary vector with the following components:
#' \tabular{rl}{
#'   \code{coordinate 1} \tab \code{adj_r_squared >= 0.4}\cr
#'   \code{coordinate 2} \tab \code{consistency_score_n > 1000000}\cr
#'   \code{coordinate 3} \tab \code{n_significant >= floor(n_bins / 10)}
#' }
#' @examples
#' n_bins <- 40
#'
#' # random spectrum
#' random_sp <- score_spectrum(runif(n = n_bins, min = -1, max = 1),
#'   max_model_degree = 1)
#' score <- classify_spectrum(
#'   get_adj_r_squared(random_sp), get_model_degree(random_sp),
#'   get_model_slope(random_sp), get_consistency_score_n(random_sp), 0, n_bins
#' )
#' sum(score)
#'
#' # non-random linear spectrum with strong noise component
#' signal <- seq(-1, 0.99, 2 / 40)
#' noise <- rnorm(n = 40, mean = 0, sd = 0.5)
#' linear_sp <- score_spectrum(signal + noise, max_model_degree = 1,
#'   max_cs_permutations = 100000)
#' score <- classify_spectrum(
#'   get_adj_r_squared(linear_sp), get_model_degree(linear_sp),
#'   get_model_slope(linear_sp), get_consistency_score_n(linear_sp), 10, n_bins
#' )
#' sum(score)
#' \dontrun{
#' # non-random linear spectrum with weak noise component
#' signal <- seq(-1, 0.99, 2 / 40)
#' noise <- rnorm(n = 40, mean = 0, sd = 0.2)
#' linear_sp <- score_spectrum(signal + noise, max_model_degree = 1,
#'   max_cs_permutations = 100000)
#' score <- classify_spectrum(
#'   get_adj_r_squared(linear_sp), get_model_degree(linear_sp),
#'   get_model_slope(linear_sp), get_consistency_score_n(linear_sp), 10, n_bins
#' )
#' sum(score)
#' }
#'
#' # non-random quadratic spectrum with strong noise component
#' signal <- seq(-1, 0.99, 2 / 40)^2 - 0.5
#' noise <- rnorm(n = 40, mean = 0, sd = 0.2)
#' quadratic_sp <- score_spectrum(signal + noise, max_model_degree = 2,
#'   max_cs_permutations = 100000)
#' score <- classify_spectrum(
#'   get_adj_r_squared(quadratic_sp), get_model_degree(quadratic_sp),
#'   get_model_slope(quadratic_sp),
#'   get_consistency_score_n(quadratic_sp), 10, n_bins
#' )
#' sum(score)
#' \dontrun{
#' # non-random quadratic spectrum with weak noise component
#' signal <- seq(-1, 0.99, 2 / 40)^2 - 0.5
#' noise <- rnorm(n = 40, mean = 0, sd = 0.1)
#' quadratic_sp <- score_spectrum(signal + noise, max_model_degree = 2)
#' score <- classify_spectrum(
#'   get_adj_r_squared(quadratic_sp), get_model_degree(quadratic_sp),
#'   get_model_slope(quadratic_sp),
#'   get_consistency_score_n(quadratic_sp), 10, n_bins
#' )
#' sum(score)
#' }
#' @family SPMA functions
#' @export
classify_spectrum <- function(adj_r_squared, degree, slope,
                              consistency_score_n,
                              n_significant, n_bins) {
    return(c(
        as.integer(adj_r_squared >= 0.4),
        as.integer(consistency_score_n > 1000000),
        as.integer(n_significant >= ceiling(n_bins / 10))
    ))
}
