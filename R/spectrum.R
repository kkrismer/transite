#' An S4 class to represent a scored spectrum
#'
#' @slot adj_r_squared adjusted \eqn{R^2} of polynomial model
#' @slot degree degree of polynomial (integer between 0 and 5)
#' @slot residuals residuals of the polynomial model
#' @slot slope coefficient of the linear term of the polynomial model
#' (spectrum "direction")
#' @slot f_statistic F statistic from the F test used to determine the
#' degree of the polynomial model
#' @slot f_statistic_p_value p-value associated with the F statistic
#' @slot consistency_score raw local consistency score of the spectrum
#' @slot consistency_score_p_value p-value associated with the local consistency
#' score
#' @slot consistency_score_n number of permutations performed to calculate
#' p-value of local consistency score (permutations performed before early
#' stopping criterion reached)
#' @slot plot spectrum plot
#' @examples
#' new("SpectrumScore",
#'     adj_r_squared = 0,
#'     degree = 0L,
#'     residuals = 0,
#'     slope = 0,
#'     f_statistic = 0,
#'     f_statistic_p_value = 1,
#'     consistency_score = 1,
#'     consistency_score_p_value = 1,
#'     consistency_score_n = 1000L,
#'     plot = NULL
#' )
#' @return Object of type SpectrumScore
#' @importFrom methods setClass
#' @name SpectrumScore-class
#' @rdname SpectrumScore-class
#' @exportClass SpectrumScore
.SpectrumScore <- setClass("SpectrumScore", slots = c(adj_r_squared = "numeric",
                                                      degree = "integer",
                                                      residuals = "numeric",
                                                      slope = "numeric",
                                                      f_statistic = "numeric",
                                                      f_statistic_p_value = "numeric",
                                                      consistency_score = "numeric",
                                                      consistency_score_p_value = "numeric",
                                                      consistency_score_n = "integer",
                                                      plot = "ANY"))

#' Getter Method spectrumAdjRSquared
#' @name SpectrumScore-class
#' @rdname SpectrumScore-class
#' @param object SpectrumScore object
#' @importFrom methods setGeneric
#' @exportMethod spectrumAdjRSquared
setGeneric("spectrumAdjRSquared", function(object) standardGeneric("spectrumAdjRSquared"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname SpectrumScore-class
#' @aliases spectrumAdjRSquared
#' @aliases spectrumAdjRSquared,SpectrumScore-method
setMethod("spectrumAdjRSquared", signature(object = "SpectrumScore"),
          function(object) object@adj_r_squared)

#' Getter Method spectrumDegree
#' @name SpectrumScore-class
#' @rdname SpectrumScore-class
#' @importFrom methods setGeneric
#' @exportMethod spectrumDegree
setGeneric("spectrumDegree", function(object) standardGeneric("spectrumDegree"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname SpectrumScore-class
#' @aliases spectrumDegree
#' @aliases spectrumDegree,SpectrumScore-method
setMethod("spectrumDegree", signature(object = "SpectrumScore"),
          function(object) object@degree)

#' Getter Method spectrumResiduals
#' @name SpectrumScore-class
#' @rdname SpectrumScore-class
#' @importFrom methods setGeneric
#' @exportMethod spectrumResiduals
setGeneric("spectrumResiduals",
           function(object) standardGeneric("spectrumResiduals"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname SpectrumScore-class
#' @aliases spectrumResiduals
#' @aliases spectrumResiduals,SpectrumScore-method
setMethod("spectrumResiduals", signature(object = "SpectrumScore"),
          function(object) object@residuals)

#' Getter Method spectrumSlope
#' @name SpectrumScore-class
#' @rdname SpectrumScore-class
#' @importFrom methods setGeneric
#' @exportMethod spectrumSlope
setGeneric("spectrumSlope", function(object) standardGeneric("spectrumSlope"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname SpectrumScore-class
#' @aliases spectrumSlope
#' @aliases spectrumSlope,SpectrumScore-method
setMethod("spectrumSlope", signature(object = "SpectrumScore"), function(object) object@slope)

#' Getter Method spectrumFStatistic
#' @name SpectrumScore-class
#' @rdname SpectrumScore-class
#' @importFrom methods setGeneric
#' @exportMethod spectrumFStatistic
setGeneric("spectrumFStatistic", function(object) standardGeneric("spectrumFStatistic"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname SpectrumScore-class
#' @aliases spectrumFStatistic
#' @aliases spectrumFStatistic,SpectrumScore-method
setMethod("spectrumFStatistic", signature(object = "SpectrumScore"),
          function(object) object@f_statistic)

#' Getter Method spectrumFStatisticPValue
#' @name SpectrumScore-class
#' @rdname SpectrumScore-class
#' @importFrom methods setGeneric
#' @exportMethod spectrumFStatisticPValue
setGeneric("spectrumFStatisticPValue",
           function(object) standardGeneric("spectrumFStatisticPValue"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname SpectrumScore-class
#' @aliases spectrumFStatisticPValue
#' @aliases spectrumFStatisticPValue,SpectrumScore-method
setMethod("spectrumFStatisticPValue", signature(object = "SpectrumScore"),
          function(object) object@f_statistic_p_value)

#' Getter Method spectrumConsistencyScore
#' @name SpectrumScore-class
#' @rdname SpectrumScore-class
#' @importFrom methods setGeneric
#' @exportMethod spectrumConsistencyScore
setGeneric("spectrumConsistencyScore",
           function(object) standardGeneric("spectrumConsistencyScore"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname SpectrumScore-class
#' @aliases spectrumConsistencyScore
#' @aliases spectrumConsistencyScore,SpectrumScore-method
setMethod("spectrumConsistencyScore", signature(object = "SpectrumScore"),
          function(object) object@consistency_score)

#' Getter Method spectrumConsistencyScorePValue
#' @name SpectrumScore-class
#' @rdname SpectrumScore-class
#' @importFrom methods setGeneric
#' @exportMethod spectrumConsistencyScorePValue
setGeneric("spectrumConsistencyScorePValue",
           function(object) standardGeneric("spectrumConsistencyScorePValue"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname SpectrumScore-class
#' @aliases spectrumConsistencyScorePValue
#' @aliases spectrumConsistencyScorePValue,SpectrumScore-method
setMethod("spectrumConsistencyScorePValue", signature(object = "SpectrumScore"),
          function(object) object@consistency_score_p_value)

#' Getter Method spectrumConsistencyScoreN
#' @name SpectrumScore-class
#' @rdname SpectrumScore-class
#' @importFrom methods setGeneric
#' @exportMethod spectrumConsistencyScoreN
setGeneric("spectrumConsistencyScoreN",
           function(object) standardGeneric("spectrumConsistencyScoreN"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname SpectrumScore-class
#' @aliases spectrumConsistencyScoreN
#' @aliases spectrumConsistencyScoreN,SpectrumScore-method
setMethod("spectrumConsistencyScoreN", signature(object = "SpectrumScore"),
          function(object) object@consistency_score_n)

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @importFrom methods is
#' @rdname SpectrumScore-class
#' @aliases show,SpectrumScore-method
setMethod("show", signature(object = "SpectrumScore"), function(object) {
    cat(is(object)[[1]], "\n",
        "  adjusted R squared: ", object@adj_r_squared, "\n",
        "  degree:  ", object@degree, "\n",
        "  residuals:  ", object@residuals, "\n",
        "  slope (linear coefficient):  ", object@slope, "\n",
        "  F statistic:  ", object@f_statistic, "\n",
        "  F statistic (p-value):  ", object@f_statistic_p_value, "\n",
        "  local consistency score:  ", object@consistency_score, "\n",
        "  local consistency score (p-value):  ", object@consistency_score_p_value, "\n",
        "  local consistency score (number of permutations):  ", object@consistency_score_n, "\n",
        sep = ""
    )
})

#' @param x SpectrumScore object
#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname SpectrumScore-class
#' @aliases plot,SpectrumScore-method
#' @exportMethod plot
setMethod("plot", signature(x = "SpectrumScore"), function(x) {
    grid::grid.draw(x@plot)
})

#' @title Subdivides Sequences into \emph{n} Bins
#'
#' @description
#' Preprocessing function for SPMA, divides transcript sequences
#' into \emph{n} bins.
#'
#' @param background_set character vector of named sequences (names are
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
#' toy_background_set <- c(
#'   "CAACAGCCUUAAUU", "CAGUCAAGACUCC", "CUUUGGGGAAU", "UCAUUUUAUUAAA",
#'   "AAUUGGUGUCUGGAUACUUCCCUGUACAU", "AUCAAAUUA", "AGAU", "GACACUUAAAGAUCCU",
#'   "UAGCAUUAACUUAAUG", "AUGGA", "GAAGAGUGCUCA", "AUAGAC", "AGUUC", "CCAGUAA"
#' )
#' # names are used as keys in the hash table (cached version only)
#' # ideally sequence identifiers (e.g., RefSeq ids) and
#' # sequence region labels (e.g., 3UTR for 3'-UTR)
#' names(toy_background_set) <- c(
#'   "NM_1_DUMMY|3UTR", "NM_2_DUMMY|3UTR", "NM_3_DUMMY|3UTR",
#'   "NM_4_DUMMY|3UTR", "NM_5_DUMMY|3UTR", "NM_6_DUMMY|3UTR",
#'   "NM_7_DUMMY|3UTR",
#'   "NM_8_DUMMY|3UTR", "NM_9_DUMMY|3UTR", "NM_10_DUMMY|3UTR",
#'   "NM_11_DUMMY|3UTR",
#'   "NM_12_DUMMY|3UTR", "NM_13_DUMMY|3UTR", "NM_14_DUMMY|3UTR"
#' )
#'
#' foreground_sets <- subdivideData(toy_background_set, n_bins = 7)
#'
#' # example data set
#' background_df <- transite:::ge$background_df
#' # sort sequences by signal-to-noise ratio
#' background_df <- dplyr::arrange(background_df, value)
#' # character vector of named sequences
#' background_set <- background_df$seq
#' names(background_set) <- paste0(background_df$refseq, "|",
#'   background_df$seq_type)
#'
#' foreground_sets <- subdivideData(background_set)
#' @family SPMA functions
#' @export
subdivideData <- function(background_set, n_bins = 40) {
    if (n_bins < 7 || n_bins > 100) {
        stop("value of bin_num is invalid, valid values are integers between 7 and 100")
    }
    if (n_bins > length(background_set)) {
        stop("more bins than transcripts")
    }

    bin_size <- floor(length(background_set) / n_bins)
    cuts <- seq_len(n_bins - 1) * bin_size
    cuts <- c(0, cuts, length(background_set))
    foreground_sets <- tapply(background_set,
                              cut(seq_len(length(background_set)),
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
#' data, the F-test
#' is used to select
#' the model that best fits the data. Since the SSR
#' monotonically decreases with
#' increasing model degree (model complexity), the relative decrease
#' of the SSR
#' between the
#' simpler model and the more complex model must outweigh the increase
#' in model complexity
#' between the two
#' models. The F-test gives the probability that a relative decrease of
#' the SSR
#' between the simpler
#' and the more complex model given their respective degrees of freedom is due
#' to chance.
#' A low p-value
#' indicates that the additional degrees of freedom of the more complex model
#' lead to a
#' better fit of
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
#' @param p_value vector of p-values (e.g., significance of enrichment
#' values) per bin
#' @param x_label label of values (e.g., \code{"enrichment value"})
#' @param midpoint for enrichment values the midpoint should be \code{1},
#' for log enrichment
#' values \code{0})
#' @param max_model_degree maximum degree of polynomial
#' @param max_cs_permutations maximum number of permutations performed in
#' Monte Carlo test
#' for consistency score
#' @param min_cs_permutations minimum number of permutations performed in
#' Monte Carlo test
#' for consistency score
#' @param e integer-valued stop criterion for consistency score Monte Carlo test: aborting
#' permutation
#' process after
#' observing \code{e} random consistency values with more extreme values than
#' the actual
#' consistency value
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
#' scoreSpectrum(runif(n = 40, min = -1, max = 1), max_model_degree = 1)
#'
#' # non-random linear spectrum
#' signal <- seq(-1, 0.99, 2 / 40)
#' noise <- rnorm(n = 40, mean = 0, sd = 0.5)
#' scoreSpectrum(signal + noise, max_model_degree = 1,
#'   max_cs_permutations = 100000)
#'
#' # non-random quadratic spectrum
#' signal <- seq(-1, 0.99, 2 / 40)^2 - 0.5
#' noise <- rnorm(n = 40, mean = 0, sd = 0.2)
#' scoreSpectrum(signal + noise, max_model_degree = 2,
#'   max_cs_permutations = 100000)
#' @family SPMA functions
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 scale_fill_gradient
#' @importFrom scales squish
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ggplot_gtable
#' @importFrom ggplot2 ggplot_build
#' @importFrom gridExtra arrangeGrob
#' @importFrom stats lm
#' @importFrom stats sd
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_smooth
#' @importFrom stats poly
#' @importFrom ggplot2 ylab
#' @importFrom stats predict
#' @importFrom stats pf
#' @importFrom methods new
#' @export
scoreSpectrum <- function(x, p_value = array(1, length(x)),
                          x_label = "log enrichment",
                          midpoint = 0,
                          max_model_degree = 3,
                          max_cs_permutations = 10000000,
                          min_cs_permutations = 5000, e = 5) {
    if (length(x) < 7) {
        print(x)
        message(x)
        stop("too few bins")
    }
    if (max_model_degree > 5 || max_model_degree < 1) {
        stop("supported values for max_model_degree: 1, 2, 3, 4, 5")
    }
    bin <- seq_len(length(x))

    df <- data.frame(
        bin = bin, value = x, p_value = p_value, RBP = rep("", length(bin)),
        stringsAsFactors = FALSE
    )

    enrichment_plot <- ggplot2::ggplot(df, ggplot2::aes_string(
        x = as.character("bin"),
        y = "RBP", fill = "value"
    )) +
        ggplot2::geom_tile(colour = "white") +
        ggplot2::scale_fill_gradient2(
            low = "#004F91FF", mid = "white", high = "#832424FF",
            midpoint = midpoint
        ) +
        ggplot2::scale_x_discrete(limits = as.character(df$bin)) +
        ggplot2::geom_text(ggplot2::aes(label = ifelse(p_value <= 0.001, "***",
                                                       ifelse(p_value <= 0.01, "**",
                                                              ifelse(p_value <= 0.05, "*", "")
                                                       )
        )),
        size = 4
        ) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(size = 8),
            axis.title.x = ggplot2::element_blank(),
            panel.background = ggplot2::element_blank(),
            legend.position = "top"
        ) +
        ggplot2::labs(fill = x_label, y = x_label)


    # determine polynomial degree of gradient
    model <- stats::lm(x ~ stats::poly(bin, degree = max_model_degree))
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

    df <- data.frame(y = x, bin = bin)
    # calculate residuals of selected model
    if (degree > 0) {
        selected_model <- stats::lm(x ~ stats::poly(bin, degree = degree))
        scatterplot <- ggplot2::ggplot(df, ggplot2::aes_string(x = "bin",
                                                               y = "y")) +
            ggplot2::geom_tile(colour = "transparent", width = 0, height = 0) +
            ggplot2::geom_point() +
            ggplot2::geom_smooth(method = "lm",
                                 formula = y ~ stats::poly(x, degree = degree)) +
            ggplot2::ylab(x_label) +
            ggplot2::scale_x_discrete(limits = as.character(df$bin)) +
            ggplot2::theme(
                axis.text.x = ggplot2::element_text(size = 8),
                axis.title.x = ggplot2::element_blank(),
                axis.title.y = ggplot2::element_blank(),
                panel.background = ggplot2::element_blank(),
                strip.background = ggplot2::element_blank()
            )
    } else {
        selected_model <- stats::lm(x ~ 1)
        scatterplot <- ggplot2::ggplot(df,
                                       ggplot2::aes_string(x = "bin",
                                                           y = "y")) +
            ggplot2::geom_tile(colour = "transparent", width = 0, height = 0) +
            ggplot2::geom_point() +
            ggplot2::geom_smooth(method = "lm", formula = y ~ 1) +
            ggplot2::ylab(x_label) +
            ggplot2::scale_x_discrete(limits = as.character(df$bin)) +
            ggplot2::theme(
                axis.text.x = ggplot2::element_text(size = 8),
                axis.title.x = ggplot2::element_blank(),
                axis.title.y = ggplot2::element_blank(),
                panel.background = ggplot2::element_blank(),
                strip.background = ggplot2::element_blank()
            )
    }

    predicted <- stats::predict(selected_model, data.frame(bin))
    residuals <- sum((predicted - x)^2)
    meta_info <- summary(selected_model)

    # calculate local consistency score
    local_consistency <- calculateLocalConsistency(x, max_cs_permutations,
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

    # combine figure
    gp1 <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(enrichment_plot))
    gp2 <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(scatterplot))

    gp2$widths <- gp1$widths

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
#' \link{scoreSpectrum}
#' @param degree degree of polynomial, returned by \link{scoreSpectrum}
#' @param slope coefficient of the linear term of the polynomial model
#' (spectrum "direction"),
#' returned by \link{scoreSpectrum}
#' @param consistency_score_n number of performed permutations before
#' early stopping,
#' returned by \link{scoreSpectrum}
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
#' random_sp <- scoreSpectrum(runif(n = n_bins, min = -1, max = 1),
#'   max_model_degree = 1)
#' score <- spectrumClassifier(
#'   spectrumAdjRSquared(random_sp), spectrumDegree(random_sp),
#'   spectrumSlope(random_sp), spectrumConsistencyScoreN(random_sp), 0, n_bins
#' )
#' sum(score)
#'
#' # non-random linear spectrum with strong noise component
#' signal <- seq(-1, 0.99, 2 / 40)
#' noise <- rnorm(n = 40, mean = 0, sd = 0.5)
#' linear_sp <- scoreSpectrum(signal + noise, max_model_degree = 1,
#'   max_cs_permutations = 100000)
#' score <- spectrumClassifier(
#'   spectrumAdjRSquared(linear_sp), spectrumDegree(linear_sp),
#'   spectrumSlope(linear_sp), spectrumConsistencyScoreN(linear_sp), 10, n_bins
#' )
#' sum(score)
#' \dontrun{
#' # non-random linear spectrum with weak noise component
#' signal <- seq(-1, 0.99, 2 / 40)
#' noise <- rnorm(n = 40, mean = 0, sd = 0.2)
#' linear_sp <- scoreSpectrum(signal + noise, max_model_degree = 1,
#'   max_cs_permutations = 100000)
#' score <- spectrumClassifier(
#'   spectrumAdjRSquared(linear_sp), spectrumDegree(linear_sp),
#'   spectrumSlope(linear_sp), spectrumConsistencyScoreN(linear_sp), 10, n_bins
#' )
#' sum(score)
#' }
#'
#' # non-random quadratic spectrum with strong noise component
#' signal <- seq(-1, 0.99, 2 / 40)^2 - 0.5
#' noise <- rnorm(n = 40, mean = 0, sd = 0.2)
#' quadratic_sp <- scoreSpectrum(signal + noise, max_model_degree = 2,
#'   max_cs_permutations = 100000)
#' score <- spectrumClassifier(
#'   spectrumAdjRSquared(quadratic_sp), spectrumDegree(quadratic_sp),
#'   spectrumSlope(quadratic_sp), spectrumConsistencyScoreN(quadratic_sp), 10, n_bins
#' )
#' sum(score)
#' \dontrun{
#' # non-random quadratic spectrum with weak noise component
#' signal <- seq(-1, 0.99, 2 / 40)^2 - 0.5
#' noise <- rnorm(n = 40, mean = 0, sd = 0.1)
#' quadratic_sp <- scoreSpectrum(signal + noise, max_model_degree = 2)
#' score <- spectrumClassifier(
#'   spectrumAdjRSquared(quadratic_sp), spectrumDegree(quadratic_sp),
#'   spectrumSlope(quadratic_sp), spectrumConsistencyScoreN(quadratic_sp), 10, n_bins
#' )
#' sum(score)
#' }
#' @family SPMA functions
#' @export
spectrumClassifier <- function(adj_r_squared, degree, slope,
                               consistency_score_n,
                               n_significant, n_bins) {
    return(c(
        as.integer(adj_r_squared >= 0.4),
        as.integer(consistency_score_n > 1000000),
        as.integer(n_significant >= ceiling(n_bins / 10))
    ))
}
