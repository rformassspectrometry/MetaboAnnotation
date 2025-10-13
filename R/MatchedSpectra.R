#' @title Representation of Spectra matches
#'
#' @name MatchedSpectra
#'
#' @aliases MatchedSpectra MatchedSpectra-class
#'
#' @description
#'
#' Matches between query and target spectra can be represented by the
#' `MatchedSpectra` object. Functions like the [matchSpectra()] function will
#' return this type of object. By default, all data accessors work as
#' *left joins* between the *query* and the *target* spectra, i.e. values are
#' returned for each *query* spectrum with eventual duplicated entries (values)
#' if the query spectrum matches more than one target spectrum.
#'
#' @section Creation, subset and filtering:
#'
#' `MatchedSpectra` objects are the result object from the [matchSpectra()].
#' While generally not needed, `MatchedSpectra` objects can also be created
#' with the `MatchedSpectra` function providing the `query` and `target`
#' `Spectra` objects as well as a `data.frame` with the *matches* between
#' query and target elements. This data frame is expected to have columns
#' `"query_idx"`, `"target_idx"` with the `integer` indices of query and
#' target objects that are *matched* and a column `"score"` with a `numeric`
#' score for the match.
#'
#' `MatchedSpectra` objects can be subset using:
#'
#' - `[` subset the `MatchedSpectra` selecting `query` spectra to keep with
#'   parameter `i`. The `target` spectra will by default be returned as-is.
#'
#' - `pruneTarget` *cleans* the `MatchedSpectra` object by removing non-matched
#'   target spectra.
#'
#' In addition, `MatchedSpectra` can be filtered with any of the filtering
#' approaches defined for [Matched()] objects: [SelectMatchesParam()],
#' [TopRankedMatchesParam()] or [ScoreThresholdParam()].
#'
#'
#' @section Extracting data:
#'
#' - `$` extracts a single spectra variable from the `MatchedSpectra` `x`. Use
#'   `spectraVariables` to get all available spectra variables. Prefix
#'   `"target_"` is used for spectra variables from the *target* `Spectra`. The
#'   matching scores are available as *spectra variable* `"score"`.
#'   Similar to a left join between the query and target spectra, this function
#'   returns a value for each query spectrum with eventual duplicated values for
#'   query spectra matching more than one target spectrum. If spectra variables
#'   from the target spectra are extracted, an `NA` is reported for *query*
#'   spectra that don't match any target spectra. See examples below for more
#'   details.
#'
#' - `length` returns the number of **query** spectra.
#'
#' - `matchedData` same as `spectraData` below.
#'
#' - `query` returns the *query* `Spectra`.
#'
#' - `queryVariables` returns the `spectraVariables` of *query*.
#'
#' - `spectraData` returns spectra variables from the query and/or target
#'   `Spectra` as a `DataFrame`. Parameter `columns` allows to define which
#'   variables should be returned (defaults to
#'   `columns = spectraVariables(object)`), spectra variable names of the target
#'   spectra need to be prefixed with `target_` (e.g. `target_msLevel` to get
#'   the MS level from target spectra). The score from the matching function is
#'   returned as spectra variable `"score"`. Similar to `$`, this function
#'   performs a *left join* of spectra variables from the *query* and *target*
#'   spectra returning all values for all query spectra (eventually returning
#'   duplicated elements for query spectra matching multiple target spectra)
#'   and the values for the target spectra matched to the respective query
#'   spectra. See help on `$` above or examples below for details.
#'
#' - `spectraVariables` returns all available spectra variables in the *query*
#'   and *target* spectra. The prefix `"target_"` is used to label spectra
#'   variables of target spectra (e.g. the name of the spectra variable for the
#'   MS level of target spectra is called `"target_msLevel"`).
#'
#' - `target` returns the *target* `Spectra`.
#'
#' - `targetVariables` returns the `spectraVariables` of *target* (prefixed
#'    with `"target_"`).
#'
#' - `whichTarget` returns an `integer` with the indices of the spectra in
#'   *target* that match at least on spectrum in *query*.
#'
#' - `whichQuery` returns an `integer` with the indices of the spectra in
#'   *query* that match at least on spectrum in *target*.
#'
#' @section Data manipulation and plotting:
#'
#' - `addProcessing`: add a processing step to both the *query* and *target*
#'   `Spectra` in `object`. Additional parameters for `FUN` can be passed *via*
#'   `...`. See `addProcessing` documentation in [Spectra::Spectra()] for more
#'   information.
#'
#' - `plotSpectraMirror`: creates a mirror plot between the query and each
#'   matching target spectrum. Can only be applied to a `MatchedSpectra` with a
#'   single query spectrum. Setting parameter `scalePeaks = TRUE` will scale
#'   the peak intensities per spectrum to a total sum of one for a better
#'   graphical visualization. Additional plotting parameters can be passed
#'   through `...`. The parameters `tolerance` and `ppm` used for matching
#'   will be passed down and reflected in the plot.
#'
#' - `setBackend`: allows to change the *backend* of both the query and target
#'   [Spectra::Spectra()] object. The function will return a `MatchedSpectra`
#'   object with the query and target `Spectra` changed to the specified
#'   `backend`, which can be any backend extending [Spectra::MsBackend].
#'
#' @param backend for `setBackend`: instance of an object extending
#'   [Spectra::MsBackend]. See help for [Spectra::setBackend()]
#'   for more details.
#'
#' @param columns for `spectraData`: `character` vector with spectra variable
#'   names that should be extracted.
#'
#' @param FUN for `addProcessing`: function to be applied to the peak matrix
#'   of each spectrum in `object`. See [Spectra::Spectra()] for more details.
#'
#' @param main for `plotSpectraMirror`: an optional title for each plot.
#'
#' @param matches `data.frame` with columns `"query_idx"` (`integer`),
#'   `"target_idx"` (`integer`) and `"score"` (`numeric`) representing the
#'   *n:m* mapping of elements between the `query` and the `target` `Spectra`.
#'
#' @param name for `$`: the name of the spectra variable to extract.
#'
#' @param object `MatchedSpectra` object.
#'
#' @param ppm For `plotSpectraMirror()`: relative m/z tolerance (in parts per
#' million, ppm) for matching peaks (see [MsCoreUtils::common()] for details).
#' If not provided by the user, the value from the `param` object is used;
#' if that is missing, the default is 20.
#'
#' @param tolerance For `plotSpectraMirror()`: absolute m/z tolerance for
#' matching peaks (see [MsCoreUtils::common()] for details).
#' If not provided by the user, the value from the `param` object is used;
#' if that is missing, the default is 0.
#'
#' @param scalePeaks for `plotSpectraMirror`: `logical(1)` if peak intensities
#'   (per spectrum) should be scaled to a total sum of one (per spectrum) prior
#'   to plotting.
#'
#' @param spectraVariables for `addProcessing`: `character` with additional
#'   spectra variables that should be passed along to the function defined
#'   with `FUN`. See [Spectra::Spectra()] for details.
#'
#' @param target `Spectra` with the spectra against which `query` has been
#'   matched.
#'
#' @param query `Spectra` with the query spectra.
#'
#' @param x `MatchedSpectra` object.
#'
#' @param xlab for `plotSpectraMirror`: the label for the x-axis.
#'
#' @param ylab for `plotSpectraMirror`: the label for the y-axis.
#'
#' @param ... for `addProcessing`: additional parameters for the function `FUN`.
#'   For `plotSpectraMirror`: additional parameters passed to the plotting
#'   functions.
#'
#' @return See individual method desciption above for details.
#'
#' @importClassesFrom Spectra Spectra
#'
#' @importFrom Spectra Spectra
#'
#' @exportClass MatchedSpectra
#'
#' @author Johannes Rainer
#'
#' @rdname MatchedSpectra
#'
#' @seealso [Matched()] for additional functions available for `MatchedSpectra`.
#'
#' @examples
#'
#' ## Creating a dummy MatchedSpectra object.
#' library(Spectra)
#' df1 <- DataFrame(
#'     msLevel = 2L, rtime = 1:10,
#'     spectrum_id = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j"))
#' df2 <- DataFrame(
#'     msLevel = 2L, rtime = rep(1:10, 20),
#'     spectrum_id = rep(c("A", "B", "C", "D", "E"), 20))
#' sp1 <- Spectra(df1)
#' sp2 <- Spectra(df2)
#' ## Define matches between query spectrum 1 with target spectra 2 and 5,
#' ## query spectrum 2 with target spectrum 2 and query spectrum 4 with target
#' ## spectra 8, 12 and 15.
#' ms <- MatchedSpectra(
#'     sp1, sp2, matches = data.frame(query_idx = c(1L, 1L, 2L, 4L, 4L, 4L),
#'                                    target_idx = c(2L, 5L, 2L, 8L, 12L, 15L),
#'                                    score = 1:6))
#'
#' ## Which of the query spectra match at least one target spectrum?
#' whichQuery(ms)
#'
#' ## Extracting spectra variables: accessor methods for spectra variables act
#' ## as "left joins", i.e. they return a value for each query spectrum, with
#' ## eventually duplicated elements if one query spectrum matches more than
#' ## one target spectrum.
#'
#' ## Which target spectrum matches at least one query spectrum?
#' whichTarget(ms)
#'
#' ## Extracting the retention times of the query spectra.
#' ms$rtime
#'
#' ## We have duplicated retention times for query spectrum 1 (matches 2 target
#' ## spectra) and 4 (matches 3 target spectra). The retention time is returned
#' ## for each query spectrum.
#'
#' ## Extracting retention times of the target spectra. Note that only retention
#' ## times for target spectra matching at least one query spectrum are returned
#' ## and an NA is reported for query spectra without matching target spectrum.
#' ms$target_rtime
#'
#' ## The first query spectrum matches target spectra 2 and 5, thus their
#' ## retention times are returned as well as the retention time of the second
#' ## target spectrum that matches also query spectrum 2. The 3rd query spectrum
#' ## does match any target spectrum, thus `NA` is returned. Query spectrum 4
#' ## matches target spectra 8, 12, and 15, thus the next reported retention
#' ## times are those from these 3 target spectra. None of the remaining 6 query
#' ## spectra matches any target spectra and thus `NA` is reported for each of
#' ## them.
#'
#' ## With `queryIndex` and `targetIndex` it is possible to extract the indices
#' ## of the matched query-index pairs
#' queryIndex(ms)
#' targetIndex(ms)
#'
#' ## The first match is between query index 1 and target index 2, the second
#' ## match between query index 1 and target index 5 and so on.
#' ## We could use these indices to extract a `Spectra` object containing only
#' ## matched target spectra and assign a spectra variable with the indices of
#' ## the query spectra
#' matched_target <- target(ms)[targetIndex(ms)]
#' matched_target$query_index <- queryIndex(ms)
#'
#' ## This `Spectra` object thus contains information from the matching, but
#' ## is a *conventional* `Spectra` object that could be used for further
#' ## analyses.
#'
#' ## `spectraData` can be used to extract all (or selected) spectra variables
#' ## from the object. Same as with `$`, a left join between the specta
#' ## variables from the query spectra and the target spectra is performed. The
#' ## prefix `"target_"` is used to label the spectra variables from the target
#' ## spectra. Below we extract selected spectra variables from the object.
#' res <- spectraData(ms, columns = c("rtime", "spectrum_id",
#'     "target_rtime", "target_spectrum_id"))
#' res
#' res$spectrum_id
#' res$target_spectrum_id
#'
#' ## Again, all values for query spectra are returned and for query spectra not
#' ## matching any target spectrum NA is reported as value for the respecive
#' ## variable.
#'
#' ## The example matched spectra object contains all query and all target
#' ## spectra. Below we subset the object keeping only query spectra that are
#' ## matched to at least one target spectrum.
#' ms_sub <- ms[whichQuery(ms)]
#'
#' ## ms_sub contains now only 3 query spectra:
#' length(query(ms_sub))
#'
#' ## while the original object contains all 10 query spectra:
#' length(query(ms))
#'
#' ## Both object contain however still the full target `Spectra`:
#' length(target(ms))
#' length(target(ms_sub))
#'
#' ## With the `pruneTarget` we can however reduce also the target spectra to
#' ## only those that match at least one query spectrum
#' ms_sub <- pruneTarget(ms_sub)
#' length(target(ms_sub))
NULL

#' @importFrom Spectra Spectra
#'
#' @noRd
setClass(
    "MatchedSpectra",
    contains = "Matched",
    prototype = prototype(
        query = Spectra(),
        target = Spectra(),
        matches = data.frame(query_idx = integer(),
                             target_idx = integer(),
                             score = numeric()),
        metadata = list(),
        version = "0.1")
)

#' @importFrom Spectra Spectra
#'
#' @export
#'
#' @rdname MatchedSpectra
MatchedSpectra <- function(query = Spectra(), target = Spectra(),
                           matches = data.frame(query_idx = integer(),
                                                target_idx = integer(),
                                                score = numeric())) {
    .matched_spectra(query = query, target = target, matches = matches)
}

.matched_spectra <- function(query = Spectra(), target = Spectra(),
                             matches = data.frame(query_idx = integer(),
                                                  target_idx = integer(),
                                                  score = numeric()),
                             metadata = list(), validate = TRUE) {
    res <- new("MatchedSpectra")
    slot(res, "query", check = FALSE) <- query
    slot(res, "target", check = FALSE) <- target
    slot(res, "matches", check = FALSE) <- matches
    slot(res, "metadata", check = FALSE) <- metadata
    if (validate)
        validObject(res)
    res
}

setValidity("MatchedSpectra", function(object) {
    msg <- .validate_matches_format(object@matches)
    if (length(msg)) return(msg)
    msg <- .validate_matches_content(object@matches, length(object@query),
                                     length(object@target))
    if (length(msg)) return(msg)
    TRUE
})

#' @importMethodsFrom ProtGenerics spectraVariables
#'
#' @importMethodsFrom BiocGenerics colnames
#'
#' @rdname MatchedSpectra
#'
#' @export
setMethod("spectraVariables", "MatchedSpectra", function(object) {
    svq <- spectraVariables(query(object))
    svt <- spectraVariables(target(object))
    cns <- colnames(object@matches)
    c(svq, paste0("target_", svt), cns[!cns %in% c("query_idx", "target_idx")])
})

#' @rdname MatchedSpectra
setMethod("queryVariables", "MatchedSpectra", function(object) {
    spectraVariables(query(object))
})

#' @rdname MatchedSpectra
setMethod("targetVariables", "MatchedSpectra", function(object) {
    paste0("target_", spectraVariables(target(object)))
})

#' @exportMethod colnames
#'
#' @rdname MatchedSpectra
setMethod("colnames", "MatchedSpectra", function(x) {
    spectraVariables(x)
})

#' @rdname MatchedSpectra
setMethod("$", "MatchedSpectra", function(x, name) {
    if (name %in% spectraVariables(x@query))
        qry <- spectraData(x@query, name)
    else qry <- spectraData(x@query, "msLevel")
    if (name %in% paste0("target_", spectraVariables(x@target)))
        trg <- spectraData(x@target, sub("^target_", "", name))
    else trg <- data.frame()
    .dollar(qry, trg, x@matches, name)
})

#' @importMethodsFrom Spectra spectraData
#'
#' @rdname MatchedSpectra
setMethod("spectraData", "MatchedSpectra",
          function(object, columns = spectraVariables(object)) {
              mtches <- object@matches
              ## Better to subset target Spectra before passing it on.
              cols_trg <- spectraVariables(object@target)
              cols_trg <- cols_trg[paste0("target_", cols_trg) %in% columns]
              if (length(cols_trg)) {
                  trg <- as.data.frame(
                      spectraData(object@target[mtches$target_idx], cols_trg))
                  mtches$target_idx <- seq_len(nrow(mtches))
              } else
                  trg <- data.frame()
              cols_qry <- columns[columns %in% spectraVariables(object@query)]
              .matchedData(
                  spectraData(object@query, unique(c("msLevel", cols_qry))),
                  trg, mtches, columns)
          })

#' @rdname MatchedSpectra
setMethod("matchedData", "MatchedSpectra",
          function(object, columns = spectraVariables(object), ...) {
              spectraData(object, columns)
          })

#' @importMethodsFrom Spectra addProcessing
#'
#' @rdname MatchedSpectra
#'
#' @export
setMethod(
    "addProcessing", "MatchedSpectra",
    function(object, FUN, ..., spectraVariables = character()) {
        if (missing(FUN))
            return(object)
        object@query <- addProcessing(object@query, FUN = FUN, ...,
                                      spectraVariables = spectraVariables)
        object@target <- addProcessing(object@target, FUN = FUN, ...,
                                       spectraVariables = spectraVariables)
        object
    })

#' @importMethodsFrom Spectra plotSpectraMirror
#'
#' @rdname MatchedSpectra
#'
#' @importFrom graphics par
#'
#' @importFrom grDevices n2mfrow
#'
#' @importFrom Spectra scalePeaks
#'
#' @importMethodsFrom ProtGenerics as.list
#'
#' @export
setMethod(
    "plotSpectraMirror", "MatchedSpectra",
    function(x, xlab = "m/z", ylab = "intensity", main = "",
             scalePeaks = FALSE, ...) {
        if (length(query(x)) != 1)
            stop("Length of 'query(x)' has to be 1.")
        dots <- list(...)
        pl <- as.list(x@metadata[["param"]])
        ppm_res <- .res_setting(dots = dots, param_list = pl,
                                 name = "ppm", default = 20)
        tol_res <- .res_setting(dots = dots, param_list = pl,
                                 name = "tolerance",default = 10)
        y <- x@target[x@matches$target_idx]
        if (!length(y))
            y <- Spectra(DataFrame(msLevel = 2L))
        if (scalePeaks) {
            x <- scalePeaks(query(x)[1L])
            y <- scalePeaks(y)
        } else x <- query(x)[1L]
        par(mfrow = n2mfrow(length(y)))
        for (i in seq_along(y))
            plotSpectraMirror(x = x, y = y[i],
                              xlab = xlab, ylab = ylab, main = main,
                              ppm = ppm_res, tolerance = tol_res,...)
    })

#' Helper for ppm and tolerance setting.
#' @noRd
.res_setting <- function(dots = list(), param_list, name, default) {
    if (!is.null(dots[[name]]))
        return(dots[[name]])
    if (!is.null(param_list[[name]]))
        return(param_list[[name]])
    return(default)
}

#' @importMethodsFrom ProtGenerics setBackend
#'
#' @rdname MatchedSpectra
#'
#' @export
setMethod("setBackend", c("MatchedSpectra", "MsBackend"),
          function(object, backend, ...) {
              object@query <- setBackend(object@query, backend = backend, ...)
              object@target <- setBackend(object@target, backend = backend, ...)
              object
          })
