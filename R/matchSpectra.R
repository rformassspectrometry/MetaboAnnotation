#' @title Spectral matching
#'
#' @description
#' 
#' The `matchSpectra` method matches (compares) spectra from `query` with those
#' from `target` based on settings specified with `param` and returns the result
#' from this as a [MatchedSpectra] object.
#'
#' @param query [Spectra] object with the (experimental) spectra.
#'
#' @param target spectral data to compare against. Can be another [Spectra].
#'
#' @param param parameter object containing the settings for the matching (e.g.
#'     eventual prefiltering settings, cut-off value for similarity above which
#'     spectra are considered matching etc).
#'
#' @param ... optional parameters.
#' 
#' @return a [MatchedSpectra] object with the spectra matching results.
#' 
#' @author Johannes Rainer
#'
#' @seealso [CompareSpectraParam()] for the comparison between [Spectra]
#'     objects.
#' 
#' @export
setGeneric("matchSpectra", function(query, target, param, ...)
           standardGeneric("matchSpectra"))

#' @title Matching MS Spectra against a reference
#'
#' @aliases CompareSpectraParam-class
#' 
#' @name CompareSpectraParam
#' 
#' @description
#' 
#' `matchSpectra` with both `query` and `target` being a [Spectra] object
#' matches each spectra in `query` against all spectra in `target` and reports
#' matches with a similarity that passes the `THRESHFUN` condition. The
#' parameters for the matching can be specified with one of the `param` objects
#' listed below:
#'
#' - `CompareSpectraParam`: the *generic* parameter object allowing to set all
#'   settings for the [compareSpectra()] call that is used to perform the
#'   similarity calculation. This includes `MAPFUN` and `FUN` defining the
#'   peak-mapping and similarity calculation functions and `ppm` and `tolerance`
#'   to define an acceptable difference between m/z values of the compared
#'   peaks. Additional parameters to the `compareSpectra` call
#'   can be passed along with `...`. See the help of [Spectra()] for more
#'   information on these parameters. Parameters `requirePrecursor` (default
#'   `TRUE`) and `requirePrecursorPeak` (default `FALSE`) allow to pre-filter
#'   the target spectra prior to the actual similarity calculation for each
#'   individual query spectrum. This can considerably improve performance.
#'   Finally, parameter `THRESHFUN` allows to define a function to be applied to
#'   the similarity scores to define which matches to report. See below for more
#'   details. 
#'
#' @param BPPARAM for `matchSpectra`: parallel processing setup (see the
#'   `BiocParallel` package for more information). By default parallel
#'   processing is disabled.
#' 
#' @param FUN `function` used to calculate similarity between spectra. Defaults
#'   for `CompareSpectraParam` to [MsCoreUtils::ndotproduct()]. See
#'   [MsCoreUtils::ndotproduct()] for details.
#'
#' @param MAPFUN `function` used to map peaks between the compared spectra.
#'    Defaults for `CompareSpectraParam` to [joinPeaks()]. See
#'    [compareSpectra()] for details.
#'
#' @param param for `matchSpectra`: parameter object (such as
#'   `CompareSpectraParam`) defining the settings for the matching.
#' 
#' @param ppm `numeric(1)` for a relative, m/z-dependent, maximal accepted
#'   difference between m/z values. This will be used in `compareSpectra` as
#'   well as for eventual precursor m/z matching.
#'
#' @param query for `matchSpectra`: [Spectra] object with the query spectra.
#' 
#' @param requirePrecursor `logical(1)` whether only target spectra are
#'   considered in the similarity calculation with a precursor m/z that matches
#'   the precursor m/z of the query spectrum (considering also `ppm` and
#'   `tolerance`). With `requirePrecursor = TRUE` (the default) the function
#'   will complete much faster, but will not find any hits for target (or query
#'   spectra) with missing precursor m/z.
#'
#' @param requirePrecursorPeak `logical(1)` whether only target spectra will be
#'   considered in the spectra similarity calculation that have a peak with an
#'   m/z matching the precursor m/z of the query spectrum. Defaults to
#'   `requirePrecursorPeak = FALSE`.
#'
#' @param target for `matchSpectra`: [Spectra] object with the target
#'   (reference) spectra to compare `query` against.
#'
#' @param tolerance `numeric(1)` for an absolute maximal accepted difference
#'   between m/z values. This will be used in `compareSpectra` as well as for
#'   eventual precursor m/z matching.
#' 
#' @param THRESHFUN `function` applied to the similarity score to define which
#'   target spectra are considered *matching*. Defaults to
#'   `THRESHFUN = function(x) which(x >= 0.7)` hence selects
#'   all target spectra matching a query spectrum with a similarity higher or
#'   equal than `0.7`. Any function that takes a numeric vector with similarity
#'   scores (as returned by [compareSpectra()]) as input and returns a
#'   `logical` or `integer` vector with the matches is supported. 
#'
#' @param ... for `CompareSpectraParam`: additional parameters passed along
#'   to the [compareSpectra()] call.
#' 
#' @author Johannes Rainer
#'
#' @importClassesFrom ProtGenerics Param
#'
#' @importClassesFrom Spectra Spectra
#' 
#' @importFrom MsCoreUtils ndotproduct
#'
#' @importFrom Spectra joinPeaks
#' 
#' @rdname CompareSpectraParam
#'
#' @exportClass CompareSpectraParam
NULL

setClass("CompareSpectraParam",
         slots = c(
             MAPFUN = "function",
             tolerance = "numeric",
             ppm = "numeric",
             FUN = "function",
             dots = "list",
             requirePrecursor = "logical",
             requirePrecursorPeak = "logical",
             THRESHFUN = "function"),
         contains = "Param",
         prototype = prototype(
             MAPFUN = joinPeaks,
             tolerance = 0,
             ppm = 5,
             FUN = MsCoreUtils::ndotproduct,
             dots = list(),
             requirePrecursor = TRUE,
             requirePrecursorPeak = FALSE,
             THRESHFUN = function(x) which(x >= 0.7)
         ),
         validity = function(object) {
             msg <- NULL
             if (length(object@tolerance) != 1 || object@tolerance < 0)
                 msg <- c("'tolerance' has to be a positive number of length 1")
             if (length(object@ppm) != 1 || object@ppm < 0)
                 msg <- c("'ppm' has to be a positive number of length 1")
             msg
         })

#' @rdname CompareSpectraParam
#'
#' @importMethodsFrom Spectra compareSpectra spectraNames<- precursorMz
#'
#' @importFrom MsCoreUtils ppm
#'
#' @importFrom BiocParallel bplapply
#' 
#' @importFrom methods new
#'
#' @export
CompareSpectraParam <- function(MAPFUN = joinPeaks, tolerance = 0, ppm = 5,
                                FUN = MsCoreUtils::ndotproduct,
                                requirePrecursor = TRUE,
                                requirePrecursorPeak = FALSE,
                                THRESHFUN = function(x) which(x >= 0.7),
                                ...) {
    new("CompareSpectraParam", MAPFUN = MAPFUN, tolerance = tolerance,
        ppm = ppm, FUN = FUN, requirePrecursor = requirePrecursor[1L],
        requirePrecursorPeak = requirePrecursorPeak[1L],
        THRESHFUN = THRESHFUN, dots = list(...))
}

.compare_spectra_parms_list <- function(x) {
    c(list(MAPFUN = x@MAPFUN, tolerance = x@tolerance, ppm = x@ppm,
           FUN = x@FUN), x@dots)
}

#' @importMethodsFrom Spectra spectraNames containsMz filterPrecursorMz
#'
#' @rdname CompareSpectraParam
#' 
#' @export
setMethod(
    "matchSpectra",
    signature(query = "Spectra", target = "Spectra",
              param = "CompareSpectraParam"),
    function(query, target, param, BPPARAM = BiocParallel::SerialParam()) {
        parms <- .compare_spectra_parms_list(param)
        snames <- spectraNames(target)
        spectraNames(target) <- seq_along(target)
        maps <- bplapply(seq_along(query), function(i, qry, trgt, parms, tf,
                                                    precMz, precMzPeak) {
            qi <- qry[i]
            if (precMz) {
                pmz <- precursorMz(qi)
                pmz <- pmz + c(-1, 1) * ppm(pmz, parms$ppm) + parms$tolerance
                trgt <- filterPrecursorMz(trgt, mz = pmz)
            }
            if (precMzPeak)
                trgt <- trgt[containsMz(trgt, mz = precursorMz(qi),
                                        ppm = parms$ppm,
                                        tolerance = parms$tolerance)]
            cor <- base::do.call(compareSpectra,
                                 c(list(x = qi, y = trgt), parms))
            keep <- tf(cor)
            data.frame(query_idx = rep(i, length(keep)),
                       target_idx = as.integer(spectraNames(trgt)[keep]),
                       score = cor[keep])
        }, qry = query, trgt = target, parms = parms,
        precMz = param@requirePrecursor,
        precMzPeak = param@requirePrecursorPeak,
        tf = param@THRESHFUN, BPPARAM = BPPARAM)
        maps <- do.call(rbind, maps)
        res <- MatchedSpectra(query, target, maps)
        res@metadata <- list(param = param)
        res
    })

