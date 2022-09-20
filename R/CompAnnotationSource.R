#' @include AllGenerics.R

#' @title Compound Annotation Sources
#'
#' @aliases CompAnnotationSource-class
#'
#' @description
#'
#' `CompAnnotationSource`s (i.e. classes extending the base virtual
#' `CompAnnotationSource` class) define and provide access to a (potentially
#' remote) compound annotation resource. This aims to simplify the integration
#' of external annotation resources by automating the actual connection
#' (or data resource download) process from the user. In addition, since the
#' reference resource is not directly exposed to the user it allows integration
#' of annotation resources that do not allow access to the full data.
#'
#' Objects extending `CompAnnotationSource` available in this package are:
#'
#' - [CompDbSource()]: annotation source referencing an annotation source in the
#'   `[CompDb()]` format ( from the `CompoundDb` Bioconductor package).
#'
#' Classes extending `CompAnnotationSource` need to implement the `matchSpectra`
#' method with parameters `query`, `target` and `param` where `query` is
#' the `Spectra` object with the (experimental) query spectra, `target` the
#' object extending the `CompAnnotationSource` and `param` the parameter object
#' defining the similarity calculation (e.g. [CompareSpectraParam()]. The method
#' is expected to return a [MatchedSpectra] object.
#'
#' `CompAnnotationSource` objects are not expected to contain any annotation
#' data. Access to the annotation data (in form of a `Spectra` object) is
#' suggested to be only established within the object's `matchSpectra` method.
#' This would also enable parallel processing of annotations as no e.g. database
#' connection would have to be shared across processes.
#'
#' @param query for `matchSpectra`: [Spectra] object with the query spectra.
#'
#' @param target for `matchSpectra`: object extending [CompAnnotationSource]
#'   (such as [CompDbSource]) with the target (reference) spectra to compare
#'   `query` against.
#'
#' @param object A `CompAnnotationSource` object.
#'
#' @param param for `matchSpectra`: parameter object (such as
#'   [CompareSpectraParam]) defining the settings for the matching.
#'
#' @param x A `CompAnnotationSource` object.
#'
#' @param ... additional parameters passed to `matchSpectra`.
#'
#' @section Methods that need to be implemented:
#'
#' For an example implementation see [CompDbSource()].
#'
#' - `matchSpectra`: function to match experimental MS2 spectra against the
#'   annotation source. See [matchSpectra()] for parameters.
#'
#' - `metadata`: function to provide metadata on the annotation resource (host,
#'   source, version etc).
#'
#' - `show` (optional): method to provide general information on the data
#'   source.
#'
#' @author Johannes Rainer, Nir Shachaf
#'
#' @name CompAnnotationSource
NULL

#' @exportClass CompAnnotationSource
setClass(
    "CompAnnotationSource",
    contains = "VIRTUAL",
    slots = c(
        version = "character"),
    prototype = prototype(version = "0.1"))

#' @rdname CompAnnotationSource
setMethod("matchSpectra",
          signature(query = "Spectra", target = "CompAnnotationSource",
                    param = "Param"),
          function(query, target, param, ...)
              stop("Not implemented for ", class(target), "."))

#' @rdname CompAnnotationSource
setMethod("show", "CompAnnotationSource", function(object) {
    cat("Object of class", class(object)[1L], "\n")
})

#' @importMethodsFrom S4Vectors metadata
#'
#' @rdname CompAnnotationSource
setMethod("metadata", "CompAnnotationSource", function(x, ...) {
    stop("Not implemented for ", class(x), ".")
})
