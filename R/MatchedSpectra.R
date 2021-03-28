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
#' return this type of object.
#'
#' @section Extracting data:
#'
#' - `length` returns the number of matches.
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
    cat("Number of query spectra: ", length(object@query), ", ",
        length(unique(object@matches$query_idx)), "matched\n", sep = "")
    cat("Number of target spectra: ", length(object@target), ", ",
        length(unique(object@matches$target_idx)), "matched\n", sep = "")
})

#' @export
target <- function(object) {
    object@target
}

#' @export
query <- function(object) {
    object@query
}

#' @export
whichTarget <- function(object) {
    unique(object@matches$query_idx)
}

#' @export
whichQuery <- function(object) {
    unique(object@matches$target_idx)
}

## Other methods
## - spectraVariables
## - spectraData
## - $
## - [
## - prune

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

## .subset_matches_nodim <- function(x, i) {
##     keepq <- match(x@matches$query_idx, i)
##     keepq <- keep[!is.na(keep)]
##     slot(x, "query", check = FALSE) <- x@query[i]
##     slot(x, "matches", check = FALSE) <- x@matches[keep, , drop = FALSE]
##     slot(x, "target", check = FALSE) <- x@target[x@matches$target_idx]
##     x
## }
