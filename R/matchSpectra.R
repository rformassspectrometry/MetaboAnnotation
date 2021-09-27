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
#' @aliases CompareSpectraParam-class MatchForwardReverseParam-class
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
#'   individual query spectrum. Target spectra can also be pre-filtered based on
#'   retention time if parameter `toleranceRt` is set to a value different than
#'   the default `toleranceRt = Inf`. Only target spectra with a retention time
#'   within the query's retention time +/- `toleranceRt` are considered. While
#'   these pre-filters can considerably improve performance, it should be noted
#'   that no matches will be found between query and target spectra with missing
#'   values in the considered variable (precursor m/z or retention time). For
#'   target spectra without retention times (such as for `Spectra` from a
#'   public reference database such as MassBank) the default `toleranceRt = Inf`
#'   should thus be used.
#'   Finally, parameter `THRESHFUN` allows to define a function to be applied to
#'   the similarity scores to define which matches to report. See below for more
#'   details.
#'
#' - `MatchForwardReverseParam`: performs spectra matching as with
#'   `CompareSpectraParam` but reports, similar to MS-DIAL, also the *reverse*
#'   similarity score and the *presence ratio*. In detail, the matching of query
#'   spectra to target spectra is performed by considering all peaks from the
#'   query and all peaks from the target (reference) spectrum (i.e. *forward*
#'   matching using an *outer join*-based peak matching strategy). For matching
#'   spectra also the *reverse* similarity is calculated considering only peaks
#'   present in the target (reference) spectrum (i.e. using a *right join*-based
#'   peak matching). This is reported as spectra variable `"reverse_score"`.
#'   In addition, the ratio between the number of matched peaks and the total
#'   number of peaks in the target (reference) spectra is reported as the
#'   *presence ratio* (spectra variable `"presence_ratio"`) and the total
#'   number of matched peaks as `"matched_peaks_count"`. See examples below
#'   for details. Parameter `THRESHFUN_REVERSE` allows to define an additional
#'   *threshold function* to filter matches. If `THRESHFUN_REVERSE` is defined
#'   only matches with a spectra similarity fulfilling both `THRESHFUN` **and**
#'   `THRESHFUN_REVERSE` are returned. With the default
#'   `THRESHFUN_REVERSE = NULL` all matches passing `THRESHFUN` are reported.
#'
#' @param BPPARAM for `matchSpectra`: parallel processing setup (see the
#'   `BiocParallel` package for more information). Parallel processing is
#'   disabled by default.
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
#'   spectra) with missing precursor m/z. It is suggested to check first the
#'   availability of the precursor m/z in `target` and `query`.
#'
#' @param requirePrecursorPeak `logical(1)` whether only target spectra will be
#'   considered in the spectra similarity calculation that have a peak with an
#'   m/z matching the precursor m/z of the query spectrum. Defaults to
#'   `requirePrecursorPeak = FALSE`. It is suggested to check first the
#'   availability of the precursor m/z in `query`, as no match will be reported
#'   for query spectra with missing precursor m/z.
#'
#' @param target for `matchSpectra`: [Spectra] object with the target
#'   (reference) spectra to compare `query` against.
#'
#' @param tolerance `numeric(1)` for an absolute maximal accepted difference
#'   between m/z values. This will be used in `compareSpectra` as well as for
#'   eventual precursor m/z matching.
#'
#' @param toleranceRt `numeric(1)` defining the maximal accepted (absolute)
#'   difference in retention time between query and target spectra. By default
#'   (with `toleranceRt = Inf`) this filter is not considered. See help of
#'   `CompareSpectraParam` above for more information.
#'
#' @param THRESHFUN `function` applied to the similarity score to define which
#'   target spectra are considered *matching*. Defaults to
#'   `THRESHFUN = function(x) which(x >= 0.7)` hence selects
#'   all target spectra matching a query spectrum with a similarity higher or
#'   equal than `0.7`. Any function that takes a numeric vector with similarity
#'   scores from the comparison of a query spectrum with all target spectra (as
#'   returned by [compareSpectra()]) as input and returns a
#'   `logical` vector (same dimensions as the similarity scores) or an integer
#'   with the matches is supported.
#'
#' @param THRESHFUN_REVERSE for `MatchForwardReverseParam`: optional additional
#'   *thresholding function* to filter the results on the reverse score. If
#'   specified the same format than `THRESHFUN` is expected.
#'
#' @param ... for `CompareSpectraParam`: additional parameters passed along
#'   to the [compareSpectra()] call.
#'
#' @author Johannes Rainer, Michael Witting
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
#'
#' @examples
#'
#' library(Spectra)
#' library(msdata)
#' fl <- system.file("TripleTOF-SWATH", "PestMix1_DDA.mzML", package = "msdata")
#' pest_ms2 <- filterMsLevel(Spectra(fl), 2L)
#'
#' ## subset to selected spectra.
#' pest_ms2 <- pest_ms2[c(808, 809, 945:955)]
#'
#' ## Load a small example MassBank data set
#' load(system.file("extdata", "minimb.RData", package = "MetaboAnnotation"))
#'
#' ## Match spectra with the default similarity score (normalized dot product)
#' csp <- CompareSpectraParam(requirePrecursor = TRUE, ppm = 10)
#' mtches <- matchSpectra(pest_ms2, minimb, csp)
#'
#' mtches
#'
#' ## Are there any matching spectra for the first query spectrum?
#' mtches[1]
#' ## No
#'
#' ## And for the second query spectrum?
#' mtches[2]
#' ## The second query spectrum matches 4 target spectra. The scores for these
#' ## matches are:
#' mtches[2]$score
#'
#' ## To access the score for the full data set
#' mtches$score
#'
#' ## Below we use a THRESHFUN that returns for each query spectrum the (first)
#' ## best matching target spectrum.
#' csp <- CompareSpectraParam(requirePrecursor = FALSE, ppm = 10,
#'     THRESHFUN = function(x) which.max(x))
#' mtches <- matchSpectra(pest_ms2, minimb, csp)
#' mtches
#'
#' ## Each of the query spectra is matched to one target spectrum
#' length(mtches)
#' matches(mtches)
#'
#' ## See the package vignette for details, descriptions and more examples.
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
             THRESHFUN = "function",
             toleranceRt = "numeric"
             ),
         contains = "Param",
         prototype = prototype(
             MAPFUN = joinPeaks,
             tolerance = 0,
             ppm = 5,
             FUN = MsCoreUtils::ndotproduct,
             dots = list(),
             requirePrecursor = TRUE,
             requirePrecursorPeak = FALSE,
             THRESHFUN = function(x) which(x >= 0.7),
             toleranceRt = Inf
         ),
         validity = function(object) {
             msg <- NULL
             if (length(object@tolerance) != 1 || object@tolerance < 0)
                 msg <- c("'tolerance' has to be a positive number of length 1")
             if (length(object@ppm) != 1 || object@ppm < 0)
                 msg <- c("'ppm' has to be a positive number of length 1")
             msg <- c(msg, .valid_threshfun(object@THRESHFUN))
             if (length(object@toleranceRt) != 1 || object@toleranceRt < 0)
                 msg <- c("'toleranceRt' has to be a positive number of",
                          " length 1")
             msg
         })

setClass("MatchForwardReverseParam",
         slots = c(THRESHFUN_REVERSE = "ANY"),
         contains = "CompareSpectraParam",
         prototype = prototype(
             THRESHFUN_REVERSE = NULL
         ),
         validity = function(object) {
             .valid_threshfun(object@THRESHFUN_REVERSE, "THRESHFUN_REVERSE")
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
                                toleranceRt = Inf, ...) {
    new("CompareSpectraParam", MAPFUN = MAPFUN, tolerance = tolerance,
        ppm = ppm, FUN = FUN, requirePrecursor = requirePrecursor[1L],
        requirePrecursorPeak = requirePrecursorPeak[1L],
        THRESHFUN = THRESHFUN, toleranceRt = toleranceRt, dots = list(...))
}

#' @rdname CompareSpectraParam
#'
#' @export
MatchForwardReverseParam <- function(MAPFUN = joinPeaks, tolerance = 0, ppm = 5,
                                     FUN = MsCoreUtils::ndotproduct,
                                     requirePrecursor = TRUE,
                                     requirePrecursorPeak = FALSE,
                                     THRESHFUN = function(x) which(x >= 0.7),
                                     THRESHFUN_REVERSE = NULL,
                                     toleranceRt = Inf, ...) {
    dots <- list(...)
    if (any(names(dots) == "type")) {
        warning("Specifying a join type with parameter 'type' is not ",
                "supported. Will ignore parameter 'type'")
        dots$type <- NULL
    }
    new("MatchForwardReverseParam", MAPFUN = MAPFUN, tolerance = tolerance,
        ppm = ppm, FUN = FUN, requirePrecursor = requirePrecursor[1L],
        requirePrecursorPeak = requirePrecursorPeak[1L],
        THRESHFUN = THRESHFUN, THRESHFUN_REVERSE = THRESHFUN_REVERSE,
        toleranceRt = toleranceRt, dots = dots)
}

.valid_threshfun <- function(x, variable = "THRESHFUN") {
    msg <- NULL
    if (length(x)) {
        if (!is.function(x))
            return(paste0(variable, " needs to be a function."))
        res <- x(1:10)
        if (!is.logical(res) & !is.integer(res))
            return(paste0(variable,
                          " needs to return integer or logical vector."))
        if (is.logical(res) && length(res) != 10)
            return(paste0(variable, " needs to return a logical vector with",
                          " the same length than target spectra."))
        if (is.integer(res) && any(res < 1 | res > 10))
            return(paste0(variable, " needs to return an integer vector with ",
                          "values between 1 and the length of target spectra."))
    }
    msg
}

.compare_spectra_parms_list <- function(x) {
    c(list(MAPFUN = x@MAPFUN, tolerance = x@tolerance, ppm = x@ppm,
           FUN = x@FUN, toleranceRt = x@toleranceRt), x@dots)
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
        if (length(query) == 1 || param@requirePrecursor ||
            param@requirePrecursorPeak || is.finite(param@toleranceRt)) {
            if (is(BPPARAM, "SerialParam"))
                .match_spectra(query, target, param)
            else .match_spectra_parallel(query, target, param, BPPARAM)
        }
        else .match_spectra_without_precursor(query, target, param)
    })

.match_spectra <- function(query, target, param) {
    parms <- .compare_spectra_parms_list(param)
    if (is.null(spectraNames(target)))
        spectraNames(target) <- seq_along(target)
    snames <- spectraNames(target)
    queryl <- length(query)
    res <- vector("list", queryl)
    for (i in seq_len(queryl)) {
        res[[i]] <- .get_matches_spectra(i, query, target, parms,
                                         param@THRESHFUN,
                                         param@requirePrecursor,
                                         param@requirePrecursorPeak,
                                         param@toleranceRt,
                                         sn = snames)
    }
    maps <- do.call(rbind, res)
    res <- MatchedSpectra(query, target, maps)
    res@metadata <- list(param = param)
    res
}

.match_spectra_parallel <- function(query, target, param, BPPARAM) {
    parms <- .compare_spectra_parms_list(param)
    if (is.null(spectraNames(target)))
        spectraNames(target) <- seq_along(target)
    snames <- spectraNames(target)
    maps <- bplapply(seq_along(query), .get_matches_spectra,
                     query = query,
                     target = target,
                     parlist = parms,
                     THRESHFUN = param@THRESHFUN,
                     precMz = param@requirePrecursor,
                     precMzPeak = param@requirePrecursorPeak,
                     toleranceRt = param@toleranceRt,
                     sn = snames, BPPARAM = BPPARAM)
    maps <- do.call(rbind, maps)
    res <- MatchedSpectra(query, target, maps)
    res@metadata <- list(param = param)
    res
}

#' This version does not use a loop but performs the comparison all in one. It
#' can thus only be performed if no precursor check/filter is used.
#'
#' @noRd
.match_spectra_without_precursor <- function(query, target, param) {
    parlist <- .compare_spectra_parms_list(param)
    cor <- do.call(compareSpectra, c(list(x = query, y = target), parlist))
    res <- lapply(seq_len(nrow(cor)), function(i) param@THRESHFUN(cor[i, ]))
    if (is.logical(res[[1L]]))
        res <- lapply(res, which)
    tidx <- unlist(res, use.names = FALSE)
    if (length(tidx)) {
        qidx <- rep(seq_along(res), lengths(res))
        maps <- data.frame(query_idx = qidx, target_idx = tidx,
                           score = cor[cbind(qidx, tidx)])
    }
    else maps <- data.frame(query_idx = integer(), target_idx = integer(),
                            score = numeric())
    res <- MatchedSpectra(query, target, maps)
    res@metadata <- list(param = param)
    res
}

#' @importMethodsFrom Spectra filterRt rtime
#'
#' @noRd
.get_matches_spectra <- function(i, query, target, parlist, THRESHFUN, precMz,
                                 precMzPeak, toleranceRt = Inf, sn) {
    qi <- query[i]
    if (is.finite(toleranceRt))
        target <- filterRt(
            target, rt = rtime(qi) + c(-toleranceRt, toleranceRt))
    if (precMz) {
        pmz <- precursorMz(qi)
        pmz <- pmz + c(-1, 1) * (ppm(pmz, parlist$ppm) + parlist$tolerance)
        target <- filterPrecursorMz(target, mz = pmz)
    }
    if (precMzPeak)
        target <- target[containsMz(target, mz = precursorMz(qi),
                                ppm = parlist$ppm,
                                tolerance = parlist$tolerance)]
    if (!length(target))
        return(data.frame(query_idx = integer(),
                          target_idx = integer(),
                          score = numeric()))
    cor <- base::do.call(compareSpectra,
                         c(list(x = qi, y = target), parlist))
    keep <- THRESHFUN(cor)
    if (is.logical(keep))
        keep <- which(keep)
    kl <- length(keep)
    data.frame(query_idx = rep(i, kl),
               target_idx = match(spectraNames(target)[keep], sn),
               score = cor[keep])
}

#' @rdname CompareSpectraParam
#'
#' @importFrom methods as
#'
#' @importMethodsFrom Spectra peaksData
#'
#' @export
setMethod(
    "matchSpectra",
    signature(query = "Spectra", target = "Spectra",
              param = "MatchForwardReverseParam"),
    function(query, target, param, BPPARAM = BiocParallel::SerialParam()) {
        res <- matchSpectra(query, target, as(param, "CompareSpectraParam"),
                            BPPARAM = BPPARAM)
        ## Loop over the matches and assign additional stuff...
        nm <- nrow(res@matches)
        res@matches$reverse_score <- rep(NA_real_, nm)
        res@matches$presence_ratio <- rep(NA_real_, nm)
        res@matches$matched_peaks_count <- rep(NA_real_, nm)
        parms_rv <- .compare_spectra_parms_list(param)
        parms_rv$type <- "right"
        for (i in seq_len(nm)) {
            spl <- c(
                list(x = peaksData(query[res@matches$query_idx[i]])[[1L]],
                     y = peaksData(target[res@matches$target_idx[i]])[[1L]]),
                parms_rv)
            map <- do.call(param@MAPFUN, spl)
            spl$x <- map$x
            spl$y <- map$y
            cor <- do.call(param@FUN, spl)
            res@matches$reverse_score[i] <- cor
            nmatched <- sum(!is.na(map$x[, 1L]))
            res@matches$presence_ratio[i] <- nmatched /
                nrow(map$y)
            res@matches$matched_peaks_count[i] <- nmatched
        }
        if (length(param@THRESHFUN_REVERSE))
            res@matches <- res@matches[param@THRESHFUN_REVERSE(
                                                 res@matches$reverse_score), ,
                                       drop = FALSE]
        res@metadata$param <- param
        res
    })
