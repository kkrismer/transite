#' Getter Method get_id
#' @name RBPMotif-class
#' @rdname RBPMotif-class
#' @param object RBPMotif object
#' @importFrom methods setGeneric
#' @exportMethod get_id
setGeneric("get_id", function(object) standardGeneric("get_id"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname RBPMotif-class
#' @aliases get_id
#' @aliases get_id,RBPMotif-method
setMethod("get_id", signature(object = "RBPMotif"),
          function(object) object@id)

#' Getter Method get_rbps
#' @name RBPMotif-class
#' @rdname RBPMotif-class
#' @importFrom methods setGeneric
#' @exportMethod get_rbps
setGeneric("get_rbps", function(object) standardGeneric("get_rbps"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname RBPMotif-class
#' @aliases get_rbps
#' @aliases get_rbps,RBPMotif-method
setMethod("get_rbps", signature(object = "RBPMotif"),
          function(object) object@rbps)

#' Getter Method get_motif_matrix
#' @name RBPMotif-class
#' @rdname RBPMotif-class
#' @importFrom methods setGeneric
#' @exportMethod get_motif_matrix
setGeneric("get_motif_matrix", function(object) standardGeneric("get_motif_matrix"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname RBPMotif-class
#' @aliases get_motif_matrix
#' @aliases get_motif_matrix,RBPMotif-method
setMethod("get_motif_matrix", signature(object = "RBPMotif"),
          function(object) object@matrix)

#' Getter Method get_hexamers
#' @name RBPMotif-class
#' @rdname RBPMotif-class
#' @importFrom methods setGeneric
#' @exportMethod get_hexamers
setGeneric("get_hexamers", function(object) standardGeneric("get_hexamers"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname RBPMotif-class
#' @aliases get_hexamers
#' @aliases get_hexamers,RBPMotif-method
setMethod("get_hexamers", signature(object = "RBPMotif"),
          function(object) object@hexamers)

#' Getter Method get_heptamers
#' @name RBPMotif-class
#' @rdname RBPMotif-class
#' @importFrom methods setGeneric
#' @exportMethod get_heptamers
setGeneric("get_heptamers",
           function(object) standardGeneric("get_heptamers"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname RBPMotif-class
#' @aliases get_heptamers
#' @aliases get_heptamers,RBPMotif-method
setMethod("get_heptamers", signature(object = "RBPMotif"),
          function(object) object@heptamers)

#' Getter Method get_width
#' @name RBPMotif-class
#' @rdname RBPMotif-class
#' @importFrom methods setGeneric
#' @exportMethod get_width
setGeneric("get_width",
           function(object) standardGeneric("get_width"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname RBPMotif-class
#' @aliases get_width
#' @aliases get_width,RBPMotif-method
setMethod("get_width", signature(object = "RBPMotif"),
          function(object) object@length)

#' Getter Method get_iupac
#' @name RBPMotif-class
#' @rdname RBPMotif-class
#' @importFrom methods setGeneric
#' @exportMethod get_iupac
setGeneric("get_iupac", function(object) standardGeneric("get_iupac"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname RBPMotif-class
#' @aliases get_iupac
#' @aliases get_iupac,RBPMotif-method
setMethod("get_iupac", signature(object = "RBPMotif"),
          function(object) object@iupac)

#' Getter Method get_type
#' @name RBPMotif-class
#' @rdname RBPMotif-class
#' @importFrom methods setGeneric
#' @exportMethod get_type
setGeneric("get_type", function(object) standardGeneric("get_type"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname RBPMotif-class
#' @aliases get_type
#' @aliases get_type,RBPMotif-method
setMethod("get_type", signature(object = "RBPMotif"),
          function(object) object@type)

#' Getter Method get_species
#' @name RBPMotif-class
#' @rdname RBPMotif-class
#' @importFrom methods setGeneric
#' @exportMethod get_species
setGeneric("get_species", function(object) standardGeneric("get_species"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname RBPMotif-class
#' @aliases get_species
#' @aliases get_species,RBPMotif-method
setMethod("get_species", signature(object = "RBPMotif"),
          function(object) object@species)

#' Getter Method get_source
#' @name RBPMotif-class
#' @rdname RBPMotif-class
#' @importFrom methods setGeneric
#' @exportMethod get_source
setGeneric("get_source", function(object) standardGeneric("get_source"))

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @rdname RBPMotif-class
#' @aliases get_source
#' @aliases get_source,RBPMotif-method
setMethod("get_source", signature(object = "RBPMotif"),
          function(object) object@src)

#' @importFrom methods setMethod
#' @importFrom methods signature
#' @importFrom methods is
#' @rdname RBPMotif-class
#' @aliases show,RBPMotif-method
setMethod("show", signature(object = "RBPMotif"), function(object) {
    cat(is(object)[[1]], "\n",
        "  id: ", object@id, "\n",
        "  RBPs:  ", paste(object@rbps, collapse = ", "), "\n",
        "  length (in nt):  ", object@length, "\n",
        "  IUPAC:  ", object@iupac, "\n",
        "  type:  ", object@type, "\n",
        "  species:  ", object@species, "\n",
        "  source:  ", object@src, "\n",
        sep = ""
    )
})

#' @param x RBPMotif object
#' @importFrom methods setMethod
#' @importFrom methods signature
#' @importFrom ggseqlogo ggseqlogo
#' @rdname RBPMotif-class
#' @aliases plot,RBPMotif-method
#' @exportMethod plot
setMethod("plot", signature(x = "RBPMotif"), function(x) {
    ppm <- t(get_ppm(x))
    return(ggseqlogo::ggseqlogo(ppm))
})
