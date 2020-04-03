#' Getter Method get_adj_r_squared
#' @name SpectrumScore-class
#' @rdname SpectrumScore-class
#' @param object SpectrumScore object
#' @importFrom methods setGeneric
#' @exportMethod get_adj_r_squared
setGeneric("get_adj_r_squared", function(object) standardGeneric("get_adj_r_squared"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname SpectrumScore-class
#' @aliases get_adj_r_squared
#' @aliases get_adj_r_squared,SpectrumScore-method
setMethod("get_adj_r_squared", signature(object = "SpectrumScore"),
          function(object) object@adj_r_squared)

#' Getter Method get_model_degree
#' @name SpectrumScore-class
#' @rdname SpectrumScore-class
#' @importFrom methods setGeneric
#' @exportMethod get_model_degree
setGeneric("get_model_degree", function(object) standardGeneric("get_model_degree"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname SpectrumScore-class
#' @aliases get_model_degree
#' @aliases get_model_degree,SpectrumScore-method
setMethod("get_model_degree", signature(object = "SpectrumScore"),
          function(object) object@degree)

#' Getter Method get_model_residuals
#' @name SpectrumScore-class
#' @rdname SpectrumScore-class
#' @importFrom methods setGeneric
#' @exportMethod get_model_residuals
setGeneric("get_model_residuals",
           function(object) standardGeneric("get_model_residuals"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname SpectrumScore-class
#' @aliases get_model_residuals
#' @aliases get_model_residuals,SpectrumScore-method
setMethod("get_model_residuals", signature(object = "SpectrumScore"),
          function(object) object@residuals)

#' Getter Method get_model_slope
#' @name SpectrumScore-class
#' @rdname SpectrumScore-class
#' @importFrom methods setGeneric
#' @exportMethod get_model_slope
setGeneric("get_model_slope", function(object) standardGeneric("get_model_slope"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname SpectrumScore-class
#' @aliases get_model_slope
#' @aliases get_model_slope,SpectrumScore-method
setMethod("get_model_slope", signature(object = "SpectrumScore"), function(object) object@slope)

#' Getter Method get_model_f_statistic
#' @name SpectrumScore-class
#' @rdname SpectrumScore-class
#' @importFrom methods setGeneric
#' @exportMethod get_model_f_statistic
setGeneric("get_model_f_statistic", function(object) standardGeneric("get_model_f_statistic"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname SpectrumScore-class
#' @aliases get_model_f_statistic
#' @aliases get_model_f_statistic,SpectrumScore-method
setMethod("get_model_f_statistic", signature(object = "SpectrumScore"),
          function(object) object@f_statistic)

#' Getter Method get_model_f_statistic_p_value
#' @name SpectrumScore-class
#' @rdname SpectrumScore-class
#' @importFrom methods setGeneric
#' @exportMethod get_model_f_statistic_p_value
setGeneric("get_model_f_statistic_p_value",
           function(object) standardGeneric("get_model_f_statistic_p_value"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname SpectrumScore-class
#' @aliases get_model_f_statistic_p_value
#' @aliases get_model_f_statistic_p_value,SpectrumScore-method
setMethod("get_model_f_statistic_p_value", signature(object = "SpectrumScore"),
          function(object) object@f_statistic_p_value)

#' Getter Method get_consistency_score
#' @name SpectrumScore-class
#' @rdname SpectrumScore-class
#' @importFrom methods setGeneric
#' @exportMethod get_consistency_score
setGeneric("get_consistency_score",
           function(object) standardGeneric("get_consistency_score"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname SpectrumScore-class
#' @aliases get_consistency_score
#' @aliases get_consistency_score,SpectrumScore-method
setMethod("get_consistency_score", signature(object = "SpectrumScore"),
          function(object) object@consistency_score)

#' Getter Method get_consistency_score_p_value
#' @name SpectrumScore-class
#' @rdname SpectrumScore-class
#' @importFrom methods setGeneric
#' @exportMethod get_consistency_score_p_value
setGeneric("get_consistency_score_p_value",
           function(object) standardGeneric("get_consistency_score_p_value"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname SpectrumScore-class
#' @aliases get_consistency_score_p_value
#' @aliases get_consistency_score_p_value,SpectrumScore-method
setMethod("get_consistency_score_p_value", signature(object = "SpectrumScore"),
          function(object) object@consistency_score_p_value)

#' Getter Method get_consistency_score_n
#' @name SpectrumScore-class
#' @rdname SpectrumScore-class
#' @importFrom methods setGeneric
#' @exportMethod get_consistency_score_n
setGeneric("get_consistency_score_n",
           function(object) standardGeneric("get_consistency_score_n"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname SpectrumScore-class
#' @aliases get_consistency_score_n
#' @aliases get_consistency_score_n,SpectrumScore-method
setMethod("get_consistency_score_n", signature(object = "SpectrumScore"),
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
