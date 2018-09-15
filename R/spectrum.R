#' @title Subdivides Sequences into \emph{n} Bins
#'
#' @description
#' Preprocessing function for SPMA, divides transcript sequences
#' into \emph{n} bins.
#'
#' @param background.set character vector of named sequences (names are
#' usually RefSeq
#' identifiers and sequence region labels,
#'  e.g., "NM_1_DUMMY|3UTR"). It is important that the sequences are
#'  already sorted by fold change,
#'  signal-to-noise ratio or any other meaningful measure.
#' @param n.bins specifies the number of bins in which the sequences
#' will be divided,
#' valid values are between 7 and 100
#' @return An array of \code{n.bins} length, containing the binned sequences
#' @examples
#' # toy example
#' toy.background.set <- c(
#'   "CAACAGCCTTAATT", "CAGTCAAGACTCC", "CTTTGGGGAAT", "TCATTTTATTAAA",
#'   "AATTGGTGTCTGGATACTTCCCTGTACAT", "ATCAAATTA", "AGAT", "GACACTTAAAGATCCT",
#'   "TAGCATTAACTTAATG", "ATGGA", "GAAGAGTGCTCA", "ATAGAC", "AGTTC", "CCAGTAA"
#' )
#' # names are used as keys in the hash table (cached version only)
#' # ideally sequence identifiers (e.g., RefSeq ids) and
#' # sequence region labels (e.g., 3UTR for 3'-UTR)
#' names(toy.background.set) <- c(
#'   "NM_1_DUMMY|3UTR", "NM_2_DUMMY|3UTR", "NM_3_DUMMY|3UTR",
#'   "NM_4_DUMMY|3UTR", "NM_5_DUMMY|3UTR", "NM_6_DUMMY|3UTR",
#'   "NM_7_DUMMY|3UTR",
#'   "NM_8_DUMMY|3UTR", "NM_9_DUMMY|3UTR", "NM_10_DUMMY|3UTR",
#'   "NM_11_DUMMY|3UTR",
#'   "NM_12_DUMMY|3UTR", "NM_13_DUMMY|3UTR", "NM_14_DUMMY|3UTR"
#' )
#'
#' foreground.sets <- subdivideData(toy.background.set, n.bins = 7)
#'
#' # exemplary data set
#' background.df <- ge$background
#' # sort sequences by signal-to-noise ratio
#' background.df <- dplyr::arrange(background.df, value)
#' # character vector of named sequences
#' background.set <- background.df$seq
#' names(background.set) <- paste0(background.df$refseq, "|",
#'   background.df$seq.type)
#'
#' foreground.sets <- subdivideData(background.set)
#' @family SPMA functions
#' @export
subdivideData <- function(background.set, n.bins = 40) {
    if (n.bins < 7 || n.bins > 100) {
        stop("value of bin.num is invalid, valid values are integers between 7 and 100")
    }
    if (n.bins > length(background.set)) {
        stop("more bins than transcripts")
    }

    bin.size <- floor(length(background.set) / n.bins)
    cuts <- seq_len(n.bins - 1) * bin.size
    cuts <- c(0, cuts, length(background.set))
    foreground.sets <- tapply(background.set,
                              cut(seq_len(length(background.set)),
                                  breaks = cuts), identity)
    return(foreground.sets)
}

SpectrumScore <- function(adj.r.squared, degree, residuals, slope,
                          f.statistic, f.statistic.p.value,
                          consistency.score, consistency.score.p.value,
                          consistency.score.n,
                          plot) {
    object <- list(
        adj.r.squared = adj.r.squared, degree = degree, residuals = residuals,
        slope = slope,
        f.statistic = f.statistic, f.statistic.p.value = f.statistic.p.value,
        consistency.score = consistency.score,
        consistency.score.p.value = consistency.score.p.value,
        consistency.score.n = consistency.score.n,
        plot = plot
    )
    class(object) <- append(class(object), "SpectrumScore")
    return(object)
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
#' \deqn{x_c = \frac{1}{n} \sum_{i = 1}^{n - 2}{\left|\frac{s_i + s_{i + 2}}{2} - s_{i + 1}\right|}.}
#' In order to obtain an estimate of the significance of a particular
#' score \eqn{x_c'},
#' Monte Carlo sampling
#' is performed by randomly permuting the coordinates of the scores
#' vector \eqn{s}
#' and recomputing \eqn{x_c}.
#' The probability estimate \eqn{\hat{p}} is given by the lower tail
#' version of the cumulative
#' distribution function
#' \deqn{\hat{Pr}(T(x)) &= \frac{\sum\limits_{i = 1}^n \mathbbm{1}\left(T(y_i) \le T(x)\right) + 1}{n + 1},}
#' where \eqn{\mathbbm{1}} is the indicator function, \eqn{n} is the
#' sample size, i.e.,
#' the number of performed permutations, and \eqn{T} equals \eqn{x_c}
#' in the above equation.
#'
#' An alternative approach to assess the consistency of a spectrum plot is
#' via polynomial regression.
#' In a first step, polynomial regression models of various degrees are
#' fitted to the data, i.e.,
#' the dependent variable \eqn{\boldsymbol{s}} (vector of scores), and
#' orthogonal
#' polynomials of the independent variable \eqn{\boldsymbol{b}} (vector of
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
#' The dependent variable \eqn{\boldsymbol{y}} is the vector of
#' scores \eqn{\boldsymbol{s}}
#' and \eqn{\boldsymbol{x}}
#' to \eqn{\boldsymbol{x}^m} are the orthogonal polynomials of
#' the vector of bin
#' numbers \eqn{\boldsymbol{b}}.
#' Orthogonal polynomials are used in order to reduce the correlation
#' between the different
#' powers of \eqn{\boldsymbol{b}}
#' and therefore avoid multicollinearity in the model. This is important,
#' because correlated
#' predictors lead to unstable coefficients, i.e., the coefficients of a
#' polynomial
#' regression model of
#' degree \eqn{m} can be greatly different from a model of degree \eqn{m + 1}.
#'
#' The orthogonal polynomials of vector \eqn{\boldsymbol{b}} are obtained by
#' centering (subtracting the mean),
#' QR decomposition, and subsequent normalization.
#' Given the dependent variable \eqn{\boldsymbol{y}} and the orthogonal
#' polynomials
#' of \eqn{\boldsymbol{b}} \eqn{\boldsymbol{x}} to
#' \eqn{\boldsymbol{x}^m}, the model coefficients \eqn{\boldsymbol{\beta}}
#' are chosen
#' in a way to minimize
#' the deviance between the actual and the predicted values characterized by
#' \deqn{\mathcal{M}(x) = \beta_0 + \beta_1 x + \beta_2 x^2 + \cdots + \beta_m x^m}
#' \deqn{\mathcal{M} = \argmin_{\mathcal{M}}{\left(\sum_{i = 1}^n{L(y_i, \mathcal{M}(x_i))}\right)},}
#' where \eqn{L(\text{actual value}, \text{predicted value})} denotes
#' the loss function.
#'
#' Ordinary least squares is used as estimation method for the model
#' coefficients \eqn{\boldsymbol{\beta}}. The loss
#' function of ordinary least squares is the sum of squared residuals
#' (\eqn{\mathit{SSR}})
#' and is defined as follows
#' \deqn{\mathit{SSR}(\boldsymbol{y},
#' \boldsymbol{\hat{y}}) = \sum_{i = 1}^n{(y_i - \hat{y}_i)^2},}
#' where \eqn{\boldsymbol{y}} are the observed data
#' and \eqn{\boldsymbol{\hat{y}}} the
#' model predictions.
#'
#' Thus the ordinary least squares estimate of the
#' coefficients \eqn{\boldsymbol{\hat{\beta}}}
#' (including the
#' intercept \eqn{\hat{\beta}_0}) of the model \eqn{\mathcal{M}} is defined by
#' \deqn{\boldsymbol{\hat{\beta}} = \argmin_{\boldsymbol{\beta}}{\left(\sum_{i = 1}^n{\left(y_i - \beta_0 - \sum_{j = 1}^m{\beta_j x_i^j}\right)^2}\right)}.}
#'
#' After polynomial models of various degrees have been fitted to the
#' data, the F-test
#' is used to select
#' the model that best fits the data. Since the \eqn{\mathit{SSR}}
#' monotonically decreases with
#' increasing model degree (model complexity), the relative decrease
#' of the \eqn{\mathit{SSR}}
#' between the
#' simpler model and the more complex model must outweigh the increase
#' in model complexity
#' between the two
#' models. The F-test gives the probability that a relative decrease of
#' the \eqn{\mathit{SSR}}
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
#' \deqn{F = \frac{(\mathit{SSR}_1 - \mathit{SSR}_2) / (p_2 - p_1)}{\mathit{SSR}_2 / (n - p_2)},}
#' where \eqn{\mathit{SSR}_i} is the sum of squared residuals
#' and \eqn{p_i} is the number
#' of parameters of
#' model \eqn{i}. The number of data points, i.e., bins, is denoted as \eqn{n}.
#' \eqn{F} is distributed according to the F-distribution
#' with \eqn{df_1 = p_2 - p_1} and \eqn{df_2 = n - p_2}.
#'
#'
#' @param x vector of values (e.g., enrichment values, normalized RBP
#' scores) per bin
#' @param p.value vector of p-values (e.g., significance of enrichment
#' values) per bin
#' @param x.label label of values (e.g., \code{"enrichment value"})
#' @param midpoint for enrichment values the midpoint should be \code{1},
#' for log enrichment
#' values \code{0})
#' @param max.model.degree maximum degree of polynomial
#' @param max.cs.permutations maximum number of permutations performed in
#' Monte Carlo test
#' for consistency score
#' @param min.cs.permutations minimum number of permutations performed in
#' Monte Carlo test
#' for consistency score
#' @param e stop criterion for consistency score Monte Carlo test: aborting
#' permutation
#' process after
#' observing \code{e} random consistency values with more extreme values than
#' the actual
#' consistency value
#' @return A list object of class \code{SpectrumScore} with the following
#' components:
#' \tabular{rl}{
#'   \code{adj.r.squared} \tab adjusted \eqn{R^2} of polynomial model\cr
#'   \code{degree} \tab maximum degree of polynomial\cr
#'   \code{residuals} \tab residuals  of polynomial model\cr
#'   \code{slope} \tab coefficient of the linear term of the polynomial model
#'    (spectrum "direction")\cr
#'   \code{f.statistic} \tab statistic of the F-test \cr
#'   \code{f.statistic.p.value} \tab p-value of F-test\cr
#'   \code{consistency.score} \tab normalized sum of deviance between the linear
#'   interpolation of the scores
#'    of two adjoining bins and the score of the middle bin, for each position
#'    in the spectrum \cr
#'   \code{consistency.score.p.value} \tab obtained by Monte Carlo sampling
#'   (randomly permuting
#'   the coordinates of the scores vector) \cr
#'   \code{consistency.score.n} \tab number of permutations\cr
#'   \code{plot} \tab
#' }
#' @examples
#' # random spectrum
#' scoreSpectrum(runif(n = 40, min = -1, max = 1), max.model.degree = 1)
#'
#' # non-random linear spectrum
#' signal <- seq(-1, 0.99, 2 / 40)
#' noise <- rnorm(n = 40, mean = 0, sd = 0.5)
#' scoreSpectrum(signal + noise, max.model.degree = 1,
#'   max.cs.permutations = 100000)
#'
#' # non-random quadratic spectrum
#' signal <- seq(-1, 0.99, 2 / 40)^2 - 0.5
#' noise <- rnorm(n = 40, mean = 0, sd = 0.2)
#' scoreSpectrum(signal + noise, max.model.degree = 2,
#'   max.cs.permutations = 100000)
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
#' @export
scoreSpectrum <- function(x, p.value = array(1, length(x)),
                          x.label = "log enrichment",
                          midpoint = 0,
                          max.model.degree = 3,
                          max.cs.permutations = 10000000,
                          min.cs.permutations = 5000, e = 5) {
    if (length(x) < 7) {
        stop("too few bins")
    }
    if (max.model.degree > 5 || max.model.degree < 1) {
        stop("supported values for max.model.degree: 1, 2, 3, 4, 5")
    }
    bin <- seq_len(length(x))

    df <- data.frame(
        bin = bin, value = x, p.value = p.value, RBP = rep("", length(bin)),
        stringsAsFactors = FALSE
    )

    enrichment.plot <- ggplot2::ggplot(df, ggplot2::aes_string(
        x = as.character("bin"),
        y = "RBP", fill = "value"
    )) +
        ggplot2::geom_tile(colour = "white") +
        ggplot2::scale_fill_gradient2(
            low = "#004F91FF", mid = "white", high = "#832424FF",
            midpoint = midpoint
        ) +
        ggplot2::scale_x_discrete(limits = as.character(df$bin)) +
        ggplot2::geom_text(ggplot2::aes(label = ifelse(p.value <= 0.001, "***",
                                                       ifelse(p.value <= 0.01, "**",
                                                              ifelse(p.value <= 0.05, "*", "")
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
        ggplot2::labs(fill = x.label, y = x.label)


    # determine polynomial degree of gradient
    model <- stats::lm(x ~ stats::poly(bin, degree = max.model.degree))
    coeffs <- summary(model)$coefficients

    if (stats::sd(x) <= 1e-10) {
        degree <- 0
    } else {
        if (max.model.degree == 5 && coeffs[6, 4] <= 0.001) {
            degree <- 5
        } else if (max.model.degree >= 4 && coeffs[5, 4] <= 0.001) {
            degree <- 4
        } else if (max.model.degree >= 3 && coeffs[4, 4] <= 0.001) {
            degree <- 3
        } else if (max.model.degree >= 2 && coeffs[3, 4] <= 0.001) {
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
        selected.model <- stats::lm(x ~ stats::poly(bin, degree = degree))
        scatterplot <- ggplot2::ggplot(df, ggplot2::aes_string(x = "bin",
                                                               y = "y")) +
            ggplot2::geom_tile(colour = "transparent", width = 0, height = 0) +
            ggplot2::geom_point() +
            ggplot2::geom_smooth(method = "lm",
                                 formula = y ~ stats::poly(x, degree = degree)) +
            ggplot2::ylab(x.label) +
            ggplot2::scale_x_discrete(limits = as.character(df$bin)) +
            ggplot2::theme(
                axis.text.x = ggplot2::element_text(size = 8),
                axis.title.x = ggplot2::element_blank(),
                axis.title.y = ggplot2::element_blank(),
                panel.background = ggplot2::element_blank(),
                strip.background = ggplot2::element_blank()
            )
    } else {
        selected.model <- stats::lm(x ~ 1)
        scatterplot <- ggplot2::ggplot(df,
                                       ggplot2::aes_string(x = "bin",
                                                           y = "y")) +
            ggplot2::geom_tile(colour = "transparent", width = 0, height = 0) +
            ggplot2::geom_point() +
            ggplot2::geom_smooth(method = "lm", formula = y ~ 1) +
            ggplot2::ylab(x.label) +
            ggplot2::scale_x_discrete(limits = as.character(df$bin)) +
            ggplot2::theme(
                axis.text.x = ggplot2::element_text(size = 8),
                axis.title.x = ggplot2::element_blank(),
                axis.title.y = ggplot2::element_blank(),
                panel.background = ggplot2::element_blank(),
                strip.background = ggplot2::element_blank()
            )
    }

    predicted <- stats::predict(selected.model, data.frame(bin))
    residuals <- sum((predicted - x)^2)
    meta.info <- summary(selected.model)

    # calculate local consistency score
    local.consistency <- calculateLocalConsistency(x, max.cs.permutations,
                                                   min.cs.permutations, e)

    if (degree > 0) {
        slope <- meta.info$coefficients[2, 1]
        f.statistic <- meta.info$fstatistic
        f.statistic.p.value <- stats::pf(f.statistic[[1]],
                                         df1 = f.statistic[[2]],
                                         df2 = f.statistic[[3]],
                                         lower.tail = FALSE
        )
    } else {
        slope <- 0
        f.statistic <- list(0, NA, NA)
        f.statistic.p.value <- 1
    }

    # combine figure
    gp1 <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(enrichment.plot))
    gp2 <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(scatterplot))

    gp2$widths <- gp1$widths

    return(SpectrumScore(
        meta.info$adj.r.squared, degree, residuals, slope,
        f.statistic[[1]], f.statistic.p.value,
        local.consistency$score, local.consistency$p.value, local.consistency$n,
        gridExtra::arrangeGrob(gp1, gp2, ncol = 1, heights = c(1, 1.4))
    ))
}

#' @title Simple spectrum classifier based on empirical thresholds
#'
#' @description
#' Spectra can be classified based on the aggregate spectrum classifier score.
#' If \code{sum(score) == 3} spectrum considered non-random, random otherwise.
#'
#' @param adj.r.squared adjusted \eqn{R^2} of polynomial model, returned by
#' \link{scoreSpectrum}
#' @param degree degree of polynomial, returned by \link{scoreSpectrum}
#' @param slope coefficient of the linear term of the polynomial model
#' (spectrum "direction"),
#' returned by \link{scoreSpectrum}
#' @param consistency.score.n number of performed permutations before
#' early stopping,
#' returned by \link{scoreSpectrum}
#' @param n.significant number of bins with statistically significant
#' enrichment
#' @param n.bins number of bins
#' @return a three-dimensional binary vector with the following components:
#' \tabular{rl}{
#'   \code{coordinate 1} \tab \code{adj.r.squared >= 0.4}\cr
#'   \code{coordinate 2} \tab \code{consistency.score.n > 1000000}\cr
#'   \code{coordinate 3} \tab \code{n.significant >= floor(n.bins / 10)}
#' }
#' @examples
#' n.bins <- 40
#'
#' # random spectrum
#' random.sp <- scoreSpectrum(runif(n = n.bins, min = -1, max = 1),
#'   max.model.degree = 1)
#' score <- spectrumClassifier(
#'   random.sp$adj.r.squared, random.sp$degree, random.sp$slope,
#'   random.sp$consistency.score.n, 0, n.bins
#' )
#' sum(score)
#'
#' # non-random linear spectrum with strong noise component
#' signal <- seq(-1, 0.99, 2 / 40)
#' noise <- rnorm(n = 40, mean = 0, sd = 0.5)
#' linear.sp <- scoreSpectrum(signal + noise, max.model.degree = 1,
#'   max.cs.permutations = 100000)
#' score <- spectrumClassifier(
#'   linear.sp$adj.r.squared, linear.sp$degree, linear.sp$slope,
#'   linear.sp$consistency.score.n, 10, n.bins
#' )
#' sum(score)
#' \dontrun{
#' # non-random linear spectrum with weak noise component
#' signal <- seq(-1, 0.99, 2 / 40)
#' noise <- rnorm(n = 40, mean = 0, sd = 0.2)
#' linear.sp <- scoreSpectrum(signal + noise, max.model.degree = 1,
#'   max.cs.permutations = 100000)
#' score <- spectrumClassifier(
#'   linear.sp$adj.r.squared, linear.sp$degree, linear.sp$slope,
#'   linear.sp$consistency.score.n, 10, n.bins
#' )
#' sum(score)
#' }
#'
#' # non-random quadratic spectrum with strong noise component
#' signal <- seq(-1, 0.99, 2 / 40)^2 - 0.5
#' noise <- rnorm(n = 40, mean = 0, sd = 0.2)
#' quadratic.sp <- scoreSpectrum(signal + noise, max.model.degree = 2,
#'   max.cs.permutations = 100000)
#' score <- spectrumClassifier(
#'   quadratic.sp$adj.r.squared, quadratic.sp$degree,
#'   quadratic.sp$slope, quadratic.sp$consistency.score.n, 10, n.bins
#' )
#' sum(score)
#' \dontrun{
#' # non-random quadratic spectrum with weak noise component
#' signal <- seq(-1, 0.99, 2 / 40)^2 - 0.5
#' noise <- rnorm(n = 40, mean = 0, sd = 0.1)
#' quadratic.sp <- scoreSpectrum(signal + noise, max.model.degree = 2)
#' score <- spectrumClassifier(
#'   quadratic.sp$adj.r.squared, quadratic.sp$degree,
#'   quadratic.sp$slope, quadratic.sp$consistency.score.n, 10, n.bins
#' )
#' sum(score)
#' }
#' @family SPMA functions
#' @export
spectrumClassifier <- function(adj.r.squared, degree, slope,
                               consistency.score.n,
                               n.significant, n.bins) {
    return(c(
        as.integer(adj.r.squared >= 0.4),
        as.integer(consistency.score.n > 1000000),
        as.integer(n.significant >= ceiling(n.bins / 10))
    ))
}
