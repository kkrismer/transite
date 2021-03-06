#' @title P-value aggregation
#'
#' @description
#' \code{p_combine} is used to combine the p-values of independent
#' significance tests.
#'
#' @details
#' The problem can be specified as follows: Given a vector of \eqn{n}
#' p-values \eqn{p_1, ..., p_n}, find \eqn{p_c}, the combined p-value of the
#' \eqn{n} significance tests.
#' Most of the methods introduced here combine the p-values in order to obtain
#' a test statistic, which follows a known probability distribution.
#' The general procedure can be stated as:
#' \deqn{T(h, C) = \sum^n_{i = 1}{h(p_i)} * C}
#' The function \eqn{T}, which returns the test statistic \eqn{t}, takes
#' two arguments.
#' \eqn{h} is a function defined on the interval \eqn{[0, 1]} that transforms
#' the individual
#' p-values, and \eqn{C} is a correction term.
#'
#' Fisher's method (1932), also known as the inverse chi-square method is
#' probably the most widely
#' used method
#' for combining p-values. Fisher used the fact that if \eqn{p_i} is
#' uniformly distributed
#' (which p-values are under the null hypothesis), then \eqn{-2 \log{p_i}}
#' follows a chi-square
#' distribution with two degrees of freedom. Therefore, if p-values are
#' transformed as follows,
#' \deqn{h(p) = -2 \log{p},}
#' and the correction term \eqn{C} is neutral, i.e., equals \eqn{1}, the
#' following statement can be
#' made about the sampling distribution of the test statistic \eqn{T_f}
#' under the null hypothesis:
#' \eqn{t_f} is distributed as chi-square with \eqn{2n} degrees of freedom,
#' where \eqn{n} is the number of p-values.
#'
#' Stouffer's method, or the inverse normal method, uses a p-value
#' transformation
#' function \eqn{h} that
#' leads to a test statistic that follows the standard normal
#' distribution by transforming
#' each p-value to its corresponding normal score. The correction
#' term scales the sum of the normal
#' scores by the root of the number of p-values.
#' \deqn{h(p) = \Phi^{-1}(1 - p)}
#' \deqn{C = \frac{1}{\sqrt{n}}}
#' Under the null hypothesis, \eqn{t_s} is distributed as standard normal.
#' \eqn{\Phi^{-1}} is the inverse of the cumulative standard
#' normal distribution function.
#'
#' An extension of Stouffer's method with weighted p-values is
#' called Liptak's method.
#'
#' The logit method by Mudholkar and George uses the following transformation:
#' \deqn{h(p) = -\ln(p / (1 - p))}
#' When the sum of the transformed p-values is corrected in the following way:
#' \deqn{C = \sqrt{\frac{3(5n + 4)}{\pi^2 n (5n + 2)}},}
#' the test statistic \eqn{t_m} is approximately t-distributed with
#' \eqn{5n + 4} degrees of freedom.
#'
#' In Tippett's method the smallest p-value is used as the test
#' statistic \eqn{t_t} and the
#' combined significance is calculated as follows:
#' \deqn{Pr(t_t) = 1 - (1 - t_t)^n}
#'
#' @param p vector of p-values
#' @param method one of the following: Fisher (1932)
#' (\code{'fisher'}), Stouffer (1949),
#' Liptak (1958) (\code{'SL'}), Mudholkar and George (1979)
#' (\code{'MG'}), and Tippett (1931)
#' (\code{'tippett'})
#' @param w weights, only used in combination with Stouffer-Liptak.
#' If \code{is.null(w)} then
#' weights are set in an unbiased way
#'
#' @return A list with the following components:
#' \tabular{rl}{
#'   \code{statistic} \tab the test statistic\cr
#'   \code{p_value} \tab the corresponding p-value\cr
#'   \code{method} \tab the method used\cr
#'   \code{statistic_name} \tab the name of the test statistic
#' }
#'
#' @examples
#' p_combine(c(0.01, 0.05, 0.5))
#'
#' p_combine(c(0.01, 0.05, 0.5), method = "tippett")
#' @importFrom stats pchisq
#' @importFrom stats qnorm
#' @importFrom stats pnorm
#' @importFrom stats pt
#' @export
p_combine <- function(p, method = c("fisher", "SL", "MG", "tippett"), w = NULL) {
    method <- match.arg(method, choices = c("fisher", "SL", "MG", "tippett"))
    p <- p[!is.na(p)]
    n <- length(p)
    if (max(p) - .Machine$double.ep > 1) {
        stop("invalid input, p > 1")
    }

    if (n == 0) {
        warning("vector of p-values is empty")
        if (method == "fisher") {
            return(list(statistic = NA, p_value = NA, method = "Fisher (1932)",
                        statistic_name = "Xsq"))
        } else if (method == "SL") {
            return(list(
                statistic = NA, p_value = NA,
                method = "Stouffer (1949), Liptak (1958)",
                statistic_name = "Z"
            ))
        } else if (method == "MG") {
            return(list(
                statistic = NA, p_value = NA,
                method = "Mudholkar and George (1979)",
                statistic_name = "L"
            ))
        } else if (method == "tippett") {
            return(list(
                statistic = NA, p_value = NA,
                method = "Tippett (1931)",
                statistic_name = "p_min"
            ))
        }
    } else {
        if (method == "fisher") {
            Xsq <- -2 * sum(log(p))
            p_val <- stats::pchisq(Xsq, df = 2 * n, lower.tail = FALSE)
            return(list(
                statistic = Xsq, p_value = p_val, method = "Fisher (1932)",
                statistic_name = "Xsq"
            ))
        } else if (method == "SL") {
            if (is.null(w)) {
                w <- rep(1, n) / n
            } else {
                if (length(w) != n) {
                    stop("length of p and w must be equal")
                }
            }
            Zi <- stats::qnorm(1 - p)
            Z <- sum(w * Zi) / sqrt(sum(w^2))
            p_value <- 1 - stats::pnorm(Z)
            return(list(
                statistic = Z, p_value = p_value,
                method = "Stouffer (1949), Liptak (1958)",
                statistic_name = "Z"
            ))
        } else if (method == "MG") {
            L <- sum((-1) * log(p / (1 - p)))
            p_value <- stats::pt(L * sqrt((15 * n + 12) / (pi^2 * n * (5 * n + 2))),
                                 df = 5 * n + 4, lower.tail = FALSE
            )
            return(list(
                statistic = L, p_value = p_value,
                method = "Mudholkar and George (1979)",
                statistic_name = "L"
            ))
        } else if (method == "tippett") {
            return(list(
                statistic = min(p), p_value = 1 - (1 - min(p))^n,
                method = "Tippett (1931)",
                statistic_name = "p_min"
            ))
        }
    }
}
