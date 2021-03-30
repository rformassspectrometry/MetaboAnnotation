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
#' @section Creation and subsetting:
#'
#' `MatchedSpectra` objects can be created with the `MatchedSpectra` function
#' providing the `query` and `target` `Spectra` as well as a `data.frame` with
#' the
#'
#' - `[` subset the `MatchedSpectra` selecting `query` spectra to keep with
#'   parameter `i`. The `target` spectra will by default be returned as-is.
#' 
#' @section Extracting data:
#'
#' - `$` extracts a single spectra variable from the `MatchedSpectra` `x`. Use
#'   `spectraVariables` to get all available spectra variables. Prefix
#'   `"target_"` is used for spectra variables from the *target* `Spectra`.
#'   Similar to a left join between the query and target spectra, this function
#'   returns a value for each query spectrum with eventual duplicated values for
#'   query spectra matching more than one target spectrum. If spectra variables
#'   from the target spectra are extracted, an `NA` is reported for *query*
#'   spectra that don't match any target spectra. See examples below for more
#'   details.
#' 
#' - `length` returns the number of matches.
#'
#' - `spectraVariables` returns all available spectra variables in the *query*
#'   and *target* spectra. The prefix `"target_"` is used to label spectra
#'   variables of target spectra (e.g. the name of the spectra variable for the
#'   MS level of target spectra is called `"target_msLevel"`).  
#' 
#' - `target` returns the *target* `Spectra`.
#' 
#' - `query` returns the *query* `Spectra`.
#'
#' - `whichTarget` returns an `integer` with the indices of the spectra in
#'   *target* that match at least on spectrum in *query*.
#'
#' - `whichQuery` returns an `integer` with the indices of the spectra in
#'   *query* that match at least on spectrum in *target*.
#'
#' @param drop for `[`: ignored.
#' 
#' @param i `integer` or `logical` defining the `query` spectra to keep.
#'
#' @param j for `[`: ignored.
#'
#' @param matches `data.frame` with columns `"query_idx"` (`integer`),
#'     `"target_idx"` (`integer`) and `"score"` (`numeric`) representing the
#'     *n:m* mapping of elements between the `query` and the `target` `Spectra`.
#'
#' @param object `MatchedSpectra` object.
#' 
#' @param target `Spectra` with the spectra against which `query` has been
#'     matched.
#'
#' @param query `Spectra` with the query spectra.
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
NULL

setClass(
    "MatchedSpectra",
    slots = c(
        query = "Spectra",
        target = "Spectra",
        matches = "data.frame",
        ## metadata
        metadata = "list",
        version = "character"
    ),
    prototype = prototype(
        query = Spectra(),
        target = Spectra(),
        matches = data.frame(query_idx = integer(),
                             target_idx = integer(),
                             score = numeric()),
        metadata = list(),
        version = "0.1")
)

#' @export
#'
#' @rdname MatchedSpectra
MatchedSpectra <- function(query = Spectra(), target = Spectra(),
                           matches = data.frame(query_idx = integer(),
                                                target_idx = integer(),
                                                score = numeric())) {
    new("MatchedSpectra", query = query, target = target, matches = matches)
}

setValidity("MatchedSpectra", function(object) {
    msg <- .validate_matches_format(object@matches)
    if (length(msg)) return(msg)
    msg <- .validate_matches_content(object@matches, length(object@query),
                                     length(object@target))
    if (length(msg)) return(msg)
    TRUE
})

#' @export
#'
#' @rdname MatchedSpectra
setMethod("length", "MatchedSpectra", function(x) nrow(x@matches))

#' @exportMethod show
#'
#' @importMethodsFrom methods show
#'
#' @rdname MatchedSpectra
setMethod("show", "MatchedSpectra", function(object) {
    cat("Object of class", class(object)[1L], "\n")
    cat("Number of matches:", length(object@matches), "\n")
    cat("Number of query spectra: ", length(object@query), " (",
        length(unique(object@matches$query_idx)), " matched)\n", sep = "")
    cat("Number of target spectra: ", length(object@target), " (",
        length(unique(object@matches$target_idx)), " matched)\n", sep = "")
})

#' @exportMethod [
#'
#' @rdname MatchedSpectra
setMethod("[", "MatchedSpectra", function(x, i, j, ..., drop = FALSE) {
    if (missing(i))
        return(x)
    if (is.logical(i))
        i <- which(i)
    .subset_matches_nodim(x, i)
})

#' @rdname MatchedSpectra
#' 
#' @export
target <- function(object) {
    object@target
}

#' @rdname MatchedSpectra
#' 
#' @export
query <- function(object) {
    object@query
}

#' @rdname MatchedSpectra
#' 
#' @export
whichTarget <- function(object) {
    unique(object@matches$target_idx)
}

#' @rdname MatchedSpectra
#' 
#' @export
whichQuery <- function(object) {
    unique(object@matches$query_idx)
}

#' @importMethodsFrom ProtGenerics spectraVariables
#'
#' @rdname MatchedSpectra
#' 
#' @export
setMethod("spectraVariables", "MatchedSpectra", function(object) {
    svq <- spectraVariables(query(object))
    svt <- spectraVariables(target(object))
    c(svq, paste0("target_", svt))
})

#' @importMethodsFrom S4Vectors $
#'
#' @rdname MatchedSpectra
#' 
#' @export
setMethod("$", "MatchedSpectra", function(x, name) {
    idxs <- .fill_index(seq_along(x@query), x@matches$query_idx)
    if (length(grep("^target_", name))) {
        vals <- do.call("$", list(x@target[x@matches$target_idx],
                                  sub("target_", "", name)))
        keep <- idxs %in% x@matches$query_idx
        idxs[keep] <- seq_len(sum(keep))
        idxs[!keep] <- NA
        vals[idxs]
    } else
        do.call("$", list(x@query, name))[idxs]
})

.fill_index <- function(x, y) {
    sort(c(setdiff(x, y), y))
}

## Other methods
## - spectraData: combine both spectraData DataFrames. Add NAs for not matching
##   query spectra.
## - $
## - pruneTarget: keep only spectra with a match in query

.validate_matches_format <- function(x) {
    msg <- NULL
    if (!is.data.frame(x))
        msg <- c(msg, "'matches' should be a 'data.frame'")
    else {
        if (!all(c("query_idx", "target_idx", "score") %in% colnames(x)))
            return(c(msg, paste0("Not all required column names \"query_idx\",",
                                 " \"target_idx\" and \"score\" found in",
                                 " 'matches'")))
        if (!is.integer(x$target_idx))
            msg <- c(msg,
                     "column \"target_idx\" is expected to be of type integer")
        if (!is.integer(x$query_idx))
            msg <- c(msg,
                     "column \"query_idx\" is expected to be of type integer")
    }
    msg
}

.validate_matches_content <- function(x, nquery, ntarget) {
    msg <- NULL
    if (!all(x$query_idx %in% seq_len(nquery)))
        msg <- c(msg, "indices in \"query_idx\" are out of bounds")
    if (!all(x$target_idx %in% seq_len(ntarget)))
        msg <- c(msg, "indices in \"target_idx\" are out of bounds")
    msg
}

#' Subsetting of a matched object with slots query, target and matches. Should
#' work on all such objects for which `target` and `query` don't have dimensions
#'
#' @param i has to be an `integer` vector with the indices of the query elements
#'     to keep.
#' 
#' @importFrom methods slot<-
#'
#' @noRd
.subset_matches_nodim <- function(x, i) {
    slot(x, "query", check = FALSE) <- x@query[i]
    mtches <- x@matches[x@matches$query_idx %in% i, , drop = FALSE]
    ## Support handling duplicated indices.
    mtches <- split.data.frame(
        x@matches, f = as.factor(x@matches$query_idx))[as.character(i)]
    lns <- vapply(mtches, function(z)
        if (length(z)) nrow(z) else 0L, integer(1))
    mtches <- do.call(rbind, mtches[lengths(mtches) > 0])
    rownames(mtches) <- NULL
    mtches$query_idx <- rep(seq_along(i), lns)
    slot(x, "matches", check = FALSE) <- mtches
    x
}
