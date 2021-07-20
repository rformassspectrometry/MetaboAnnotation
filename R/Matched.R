#' @title Representation of generic objects matches
#'
#' @name Matched
#'
#' @aliases Matched Matched-class [,Matched-method
#'
#' @description
#'
#' Matches between query and target generic objects can be represented by the
#' `Matched` object. By default, all data accessors work as
#' *left joins* between the *query* and the *target* object, i.e. values are
#' returned for each *query* object with eventual duplicated entries (values)
#' if the query object matches more than one target object. See also
#' *Creation and subsetting* as well as *Extracting data* sections below for
#' details and more information
#'
#' @section Creation and subsetting:
#'
#' `Matched` objects can be created with the `Matched` function
#' providing the `query` and `target` object as well as the `matches`
#' `data.frame` with two columns of integer indices defining which elements
#' from *query* match which element from *target*.
#'
#' A `Matched` object can be subsetted in the following ways:
#'
#' - `[` subset the `Matched` selecting `query` object elements to keep with
#'   parameter `i`. The resulting `Matched` object will contain all the matches
#'   for the selected query elements. The `target` object will by default be
#'   returned as-is.
#'
#' - `pruneTarget` *cleans* the `Matched` object by removing non-matched
#'   target elements.
#'
#' @section Extracting data:
#'
#' - `$` extracts a single variable from the `Matched` `x`. The variables that
#'   can be extracted can be listed using `colnames(x)`. These variables can
#'   belong to *query*, *target* or be related to the matches (e.g. the
#'   score of each match). If the *query* (*target*) object is two dimensional,
#'   its columns can be extracted (prefix` "target_"` is used for columns in the
#'   *target* object) otherwise (if *query* respectively *target* does have only
#'   a single dimension (e.g. is a `list` or a `character`) the whole object can
#'   is extracted with `x$query` (`x$target`).
#'   The matching scores are available as *variable*
#'   `"score"`. Similar to a left join between the query and target elements,
#'   this function returns a value for each query element, with eventual
#'   duplicated values for query elements matching more than one target element.
#'   If variables from the target `data.frame` are extracted, an `NA` is
#'   reported for the entries corresponding to *query* elements that don't match
#'   any target element. See examples below for more details.
#'
#' - `length` returns the number of **query** elements.
#'
#' - `matchedData` allows to extract multiple variables contained in the
#'   `Matched` object as a `DataFrame`. Parameter `columns` allows to
#'   define which columns (or variables) should be returned (defaults to
#'   `columns = colnames(object)`). Each single column in the returned
#'   `DataFrame` is constructed in the same way as in `$`. That is, like `$`,
#'   this function performs a *left join* of variables from the *query* and
#'   *target* objects returning all values for all query elements (eventually
#'   returning duplicated elements for query elements matching multiple target
#'   elements) and the values for the target elements matched to the respective
#'   query elements (or `NA` if the target element is not matched to any query
#'   element).
#'
#' - `target` returns the *target* object.
#'
#' - `query` returns the *query* object.
#'
#' - `whichTarget` returns an `integer` with the indices of the elements in
#'   *target* that match at least one element in *query*.
#'
#' - `whichQuery` returns an `integer` with the indices of the elements in
#'   *query* that match at least one element in *target*.
#'
#' @param columns for `matchedData`: `character` vector with column names of
#'   variables that should be extracted.
#'
#' @param drop for `[`: ignored.
#'
#' @param i `integer` or `logical` defining the `query` elements to keep.
#'
#' @param j for `[`: ignored.
#'
#' @param matches `data.frame` with columns `"query_idx"` (`integer`),
#'   `"target_idx"` (`integer`) and `"score"` (`numeric`) representing the n:m
#'   mapping of elements between the `query` and the `target` objects.
#'
#' @param name for `$`: the name of the column (or variable) to extract.
#'
#' @param object a `Matched` object.
#'
#' @param target object with the elements against which `query` has been
#'   matched.
#'
#' @param query object with the query elements.
#'
#' @param x `Matched` object.
#'
#' @param ... additional parameters.
#'
#' @return See individual method description above for details.
#'
#' @exportClass Matched
#'
#' @author Andrea Vicini, Johannes Rainer
#'
#' @rdname Matched
#'
#' @examples
#'
#' ## Creating a Matched object.
#' q1 <- data.frame(col1 = 1:5, col2 = 6:10)
#' t1 <- data.frame(col1 = 11:16, col2 = 17:22)
#' ## Define matches between query row 1 with target row 2 and, query row 2
#' ## with target rows 2,3,4 and query row 5 with target row 5.
#' mo <- Matched(
#'     q1, t1, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
#'                                 target_idx = c(2L, 2L, 3L, 4L, 5L),
#'                                  score = seq(0.5, 0.9, by = 0.1)))
#' mo
#'
#' ## Which of the query elements (rows) match at least one target
#' ## element (row)?
#' whichQuery(mo)
#'
#' ## Which target elements (rows) match at least one query element (row)?
#' whichTarget(mo)
#'
#' ## Extracting variable "col1" from query object .
#' mo$col1
#'
#' ## We have duplicated values for the entries of `col1` related to query
#' ## elements (rows) matched to multiple rows of the target object). The
#' ## value of `col1` is returned for each element (row) in the query.
#'
#' ## Extracting variable "col1" from target object. To access columns from
#' ## target we have to prefix the name of the column by `"target_"`.
#' ## Note that only values of `col1` for rows matching at least one query
#' ## row are returned and an NA is reported for query rows without matching
#' ## target rows.
#' mo$target_col1
#'
#' ## The 3rd and 4th query rows do not match any target row, thus `NA` is
#' ## returned.
#'
#' ## `matchedData` can be used to extract all (or selected) columns
#' ## from the object. Same as with `$`, a left join between the columns
#' ## from the query and the target is performed. Below we extract selected
#' ## columns from the object as a DataFrame.
#' res <- matchedData(mo, columns = c("col1", "col2", "target_col1", "target_col2"))
#' res
#' res$col1
#' res$target_col1
#'
#' ## The example matched object contains all query and all target
#' ## elements (rows). Below we subset the object keeping only query rows that
#' ## are matched to at least one target row.
#' mo_sub <- mo[whichQuery(mo)]
#'
#' ## mo_sub contains now only 3 query rows:
#' nrow(query(mo_sub))
#'
#' ## while the original object contains all 5 query rows:
#' nrow(query(mo))
#'
#' ## Both objects contain however still the full target object:
#' nrow(target(mo))
#' nrow(target(mo_sub))
#'
#' ## With the `pruneTarget` we can however reduce also the target rows to
#' ## only those that match at least one query row
#' mo_sub <- pruneTarget(mo_sub)
#' nrow(target(mo_sub))
#'
#' ## Creating a Matched object with a `data.frame` for `query` and a `vector`
#' ## for `target`. The matches are specified in the same way as the example
#' ## before.
#'
#' q1 <- data.frame(col1 = 1:5, col2 = 6:10)
#' t2 <- 11:16
#' mo <- Matched(q1, t2, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
#'     target_idx = c(2L, 2L, 3L, 4L, 5L), score = seq(0.5, 0.9, by = 0.1)))
#'
#' ## *target* is a simple vector and has thus no columns. The matched values
#' ## from target, if it does not have dimensions and hence column names, can
#' ## be retrieved with `$target`
#' mo$target
#'
#' ## Note that in this case "target" is returned by the function `colnames`
#' colnames(mo)
#'
#' ## As before, we can extract all data as a `DataFrame`
#' res <- matchedData(mo)
#' res
#'
#' ## Note that the columns of the obtained `DataFrame` are the same as the
#' ## corresponding vectors obtained with `$`
#'
#' res$col1
#' res$target
#'
#' ## Also subsetting and pruning works in the same way as the example above.
#'
#' mo_sub <- mo[whichQuery(mo)]
#'
#' ## mo_sub contains now only 3 query rows:
#' nrow(query(mo_sub))
#'
#' ## while the original object contains all 5 query rows:
#' nrow(query(mo))
#'
#' ## Both object contain however still the full target object:
#' length(target(mo))
#' length(target(mo_sub))
#'
#' ## Reducing the target elements to only those that match at least one query row
#' mo_sub <- pruneTarget(mo_sub)
#' nrow(target(mo_sub))
NULL

setClass(
    "Matched",
    slots = c(
        query = "ANY",
        target = "ANY",
        matches = "data.frame",
        metadata = "list",
        version = "character"
    ),
    prototype = prototype(
        query = list(),
        target = list(),
        matches = data.frame(query_idx = integer(),
                             target_idx = integer(),
                             score = numeric()),
        metadata = list(),
        version = "0.1")
)

#' @export
#'
#' @rdname Matched
Matched <- function(query = list(), target = list(),
                    matches = data.frame(query_idx = integer(),
                                         target_idx = integer(),
                                         score = numeric())) {
    new("Matched", query = query, target = target, matches = matches)
}

setValidity("Matched", function(object) {
    msg <- .validate_matches_format(object@matches)
    if (length(msg)) return(msg)
    msg <- .validate_matches_content(object@matches, .nelements(object@query),
                                     .nelements(object@target))
    if(any(c("query", "target") %in% colnames(object@matches)))
        return("\"query\" and \"target\" can't be used as matches column names")
    if (length(msg)) return(msg)
    msg <- .validate_qt(object@query)
    if (length(msg)) return(msg)
    msg <- .validate_qt(object@target)
    if (length(msg)) return(msg)
    TRUE
})

#' @export
#'
#' @rdname Matched
setMethod("length", "Matched", function(x) .nelements(x@query))

#' @exportMethod show
#'
#' @importMethodsFrom methods show
#'
#' @rdname Matched
setMethod("show", "Matched", function(object) {
    cat("Object of class", class(object)[1L], "\n")
    cat("Total number of matches:", nrow(object@matches), "\n")
    cat("Number of query objects: ", .nelements(object@query),
        " (", length(unique(object@matches$query_idx)), " matched)\n", sep = "")
    cat("Number of target objects: ", .nelements(object@target),
        " (", length(unique(object@matches$target_idx)), " matched)\n", sep = "")
})

#' @exportMethod [
#'
#' @rdname Matched
setMethod("[", "Matched", function(x, i, j, ..., drop = FALSE) {
    if (missing(i))
        return(x)
    if (is.logical(i))
        i <- which(i)
    .subset_matches_nodim(x, i)
})

#' @rdname Matched
#'
#' @export
target <- function(object) {
    object@target
}

#' @rdname Matched
#'
#' @export
query <- function(object) {
    object@query
}

#' @rdname Matched
#'
#' @export
whichTarget <- function(object) {
    unique(object@matches$target_idx)
}

#' @rdname Matched
#'
#' @export
whichQuery <- function(object) {
    unique(object@matches$query_idx)
}

#' @importMethodsFrom S4Vectors $
#'
#' @rdname Matched
#'
#' @export
setMethod("$", "Matched", function(x, name) {
  .dollar2(x@query, x@target, x@matches, name)
})

#' @importFrom BiocGenerics colnames
#'
#' @exportMethod colnames
#'
#' @rdname Matched
setMethod("colnames", "Matched", function(x) {
  .colnames(x@query, x@target, x@matches)
})

#' @rdname Matched
#'
#' @exportMethod matchedData
setGeneric("matchedData", function(object, ...)
    standardGeneric("matchedData"))

#' @importMethodsFrom S4Vectors cbind
#'
#' @importFrom S4Vectors DataFrame
#'
#' @rdname Matched
#'
#' @export
setMethod("matchedData", "Matched", function(object,
                                             columns = colnames(object), ...) {
    .matchedData2(object@query, object@target, object@matches, columns, ... )
})


.nelements <- function(x) {
    ifelse(is.null(d <- dim(x)), length(x), d[1])
}

#' @importFrom methods is
#'
#' @noRd
.extract_elements <- function(x, i, j, drop = FALSE) {
    if (length(dim(x))) {
        if (missing(j)) j <- seq_len(dim(x)[2])
        res <- x[i, j, drop = drop]
    } else res <- x[i]
    if(is(x)[1] == "list" && any(na <- is.na(i)))
        res[na] <- NA
    res
}

#' Subsetting of a matched object with slots query, target and matches.
#'
#' @param i has to be an `integer` vector with the indices of the query elements
#'     to keep.
#'
#' @importFrom methods slot<-
#'
#' @noRd
.subset_matches_nodim <- function(x, i) {
    if(!all(i %in% seq_len(length(x))))
        stop("subscript contains out-of-bounds indices")
    slot(x, "query", check = FALSE) <- .extract_elements(x@query, i)
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

.fill_index <- function(x, y) {
    sort(c(setdiff(x, y), y))
}

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

.validate_qt <- function(x){
    msg <- NULL
    ndim <- length(dim(x))
    if(!(ndim %in% c(0,2)))
        msg <- c(msg, "unsupported dimensions in either \"query\" or \"target\"")
    if (ndim == 2 && !inherits(x, "SummarizedExperiment") && is.null(colnames(x)))
        msg <- c(msg, paste0("either \"query\" or \"target\" have 2",
                 " dimensions but no column names"))
    msg
}

.colnames <- function(query, target, matches) {
  cns <- colnames(matches)
  cnq <- NULL
  cnt <- NULL
  if (length(dim(target)) == 2){
    cnt <- colnames(target)
    if (length(cnt)) cnt <- paste0("target_", cnt)
  }
  if (is.null(dim(target)))
    cnt <- "target"
  if (length(dim(query)) == 2)
    cnq <- colnames(query)
  if (is.null(dim(query)))
    cnq <- "query"
  c(cnq, cnt, cns[!cns %in% c("query_idx", "target_idx")])
}

.dollar <- function(query, target, matches, name) {
  if(name %in% .colnames(query, target, matches))
  {
    idxs <- .fill_index(seq_len(.nelements(query)), matches$query_idx)
    if (name %in% colnames(matches)) {
      keep <- idxs %in% matches$query_idx
      idxs[keep] <- seq_len(sum(keep))
      idxs[!keep] <- NA
      return(matches[idxs, name])
    }
    if (name == "query" && is.null(dim(query)))
      return(query[idxs])
    if (name == "target" && is.null(dim(target))){
      vals <- target[matches$target_idx]
      keep <- idxs %in% matches$query_idx
      idxs[keep] <- seq_len(sum(keep))
      idxs[!keep] <- NA
      return(.extract_elements(vals, idxs))
    }
    if (length(grep("^target_", name))) {
      vals <- target[matches$target_idx, sub("target_", "", name)]
      keep <- idxs %in% matches$query_idx
      idxs[keep] <- seq_len(sum(keep))
      idxs[!keep] <- NA
      vals[idxs]
    } else
      query[idxs, name]
  }
}

.dollar2 <- function(query, target, matches, name) {
  if(name %in% .colnames(query, target, matches))
  {
    not_mtchd <- setdiff(seq_len(.nelements(query)), matches$query_idx)
    idxs_qry <- c(matches$query_idx, not_mtchd)
    ord <- order(idxs_qry) # Needed only if we want the entries of the returned vector 
    # in the same way as query is ordered but in principle we could avoid the ordering
    if (name %in% colnames(matches)) {
      idxs_mtch <- c(seq_len(nrow(matches)), rep(NA, length(not_mtchd)))[ord]
      return(matches[idxs_mtch, name])
    }
    if (name == "target" || length(grep("^target_", name))){
      idxs_trg <- c(matches$target_idx, rep(NA, length(not_mtchd)))[ord] 
      .extract_elements(target, idxs_trg, sub("target_", "", name), drop = TRUE)
    }else
      .extract_elements(query, idxs_qry[ord], name, drop = TRUE)
  }
}


.matchedData <- function(query, target, matches, columns, ...) {
  cnms <- .colnames(query, target, matches)
  if (any(!columns %in% cnms))
    stop("column(s) ", paste0(columns[!columns %in% cnms],
                              collapse = ", "), " not available")
  idxs <- .fill_index(seq_len(.nelements(query)), matches$query_idx)
  from_target <- grepl("^target_", columns) | columns == "target"
  from_matches <- columns %in% colnames(matches)
  from_query <- !(from_target | from_matches)
  res_q <- NULL
  res_t <- NULL
  res_m <- NULL
  if (any(from_query))
    res_q <- .extract_elements(query, idxs, columns[from_query])
  if (any(from_target)) {
    keep <- idxs %in% matches$query_idx
    target_idxs <- idxs
    target_idxs[keep] <- seq_len(sum(keep))
    target_idxs[!keep] <- NA
    res_t <- .extract_elements(target,
                               matches$target_idx[target_idxs],
                               sub("target_", "", columns[from_target]))
  }
  if (any(from_matches)) {
    keep <- idxs %in% matches$query_idx
    idxs[keep] <- seq_len(sum(keep))
    idxs[!keep] <- NA
    res_m <- matches[idxs, columns[from_matches], drop = FALSE]
  }
  if(!is.null(res_q) && is.null(dim(query))) res_q <- I(res_q)
  if(!is.null(res_t) && is.null(dim(target))) res_t <- I(res_t)
  any_qtm <- c(any(from_query), any(from_target), any(from_matches))
  res <- DataFrame(do.call(cbind, list(res_q, res_t, res_m)[any_qtm]))
  colnames(res) <- c(columns[from_query], columns[from_target],
                     columns[from_matches])
  res[, columns, drop = FALSE]
}


.matchedData2 <- function(query, target, matches, columns, ...) {
  cnms <- .colnames(query, target, matches)
  if (any(!columns %in% cnms))
    stop("column(s) ", paste0(columns[!columns %in% cnms],
                              collapse = ", "), " not available")
  not_mtchd <- setdiff(seq_len(.nelements(query)), matches$query_idx)
  idxs_qry <- c(matches$query_idx, not_mtchd)
  ord <- order(idxs_qry) # Needed only if we want the rows of the returned matrix ordered 
  # in the same way as query is ordered but in principle we could avoid the ordering 
  from_target <- grepl("^target_", columns) | columns == "target"
  from_matches <- columns %in% colnames(matches)
  from_query <- !(from_target | from_matches)
  res_q <- NULL
  res_t <- NULL
  res_m <- NULL
  if (any(from_query))
    res_q <- .extract_elements(query, idxs_qry[ord], columns[from_query])
  if (any(from_target)) {
    idxs_trg <- c(matches$target_idx, rep(NA, length(not_mtchd)))[ord] 
    res_t <- .extract_elements(target, idxs_trg, 
                               sub("target_", "", columns[from_target]))
  }
  if (any(from_matches)) {
    idxs_mtch <- c(seq_len(nrow(matches)), rep(NA, length(not_mtchd)))[ord]
    res_m <- matches[idxs_mtch, columns[from_matches], drop = FALSE]
  }
  if(!is.null(res_q) && is.null(dim(query))) res_q <- I(res_q)
  if(!is.null(res_t) && is.null(dim(target))) res_t <- I(res_t)
  any_qtm <- c(any(from_query), any(from_target), any(from_matches))
  res <- DataFrame(do.call(cbind, list(res_q, res_t, res_m)[any_qtm]))
  colnames(res) <- c(columns[from_query], columns[from_target],
                     columns[from_matches])
  res[, columns, drop = FALSE]
}

#' @rdname Matched
#'
#' @importFrom methods validObject
#'
#' @export
pruneTarget <- function(object) {
    keep <- whichTarget(object)
    object@target <- .extract_elements(object@target, keep)
    object@matches$target_idx <- match(object@matches$target_idx, keep)
    validObject(object)
    object
}


removeMatches <- function(object, idxs)
{
    object@matches <- object@matches[-idxs, ]
    validObject(object)
    object
}

keepMatches <- function(object, idxs)
{
  object@matches <- object@matches[idxs, ]
  validObject(object)
  object
}
