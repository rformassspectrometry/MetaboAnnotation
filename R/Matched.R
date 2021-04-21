#' @title Representation of generic objects matches
#'
#' @name Matched
#' 
#' @aliases Matched Matched-class
#' 
#' @description
#' 
#' Matches between query and target generic objects can be represented by the
#' `Matched` object. By default, all data accessors work as
#' *left joins* between the *query* and the *target* object, i.e. values are
#' returned for each *query* object with eventual duplicated entries (values)
#' if the query object matches more than one target object. 
#'
#' @section Creation and subsetting:
#'
#' `Matched` objects can be created with the `Matched` function
#' providing the `query` and `target` object as well as a `data.frame` with
#' the
#'
#' - `[` subset the `Matched` selecting `query` object elements to keep with
#'   parameter `i`. The `target` object will by default be returned as-is.
#' 
#' - `pruneTarget` *cleans* the `Matched` object by removing non-matched
#'   target elements.
#' 
#' @section Extracting data:
#'
#' - `$` extracts a single column from the `Matched` `x` if the objects of the 
#'   match are `data.frame`s. Prefix
#'   `"target_"` is used for columns in the *target* `data.frame`. The
#'   matching scores are available as *variable* `"score"`.
#'   Similar to a left join between the query and target spectra, this function
#'   returns a value for each query element with eventual duplicated values for
#'   query elements matching more than one target element. If columns
#'   from the target `data.frame` are extracted, an `NA` is reported for *query*
#'   rows that don't match any target rows. See examples below for more
#'   details.
#' 
#' - `length` returns the number of **query** elements.
#'
#' - `data` extracts the data contined in the `Matched` object as a `DataFrame`. 
#'   Parameter `columns` allows to define which
#'   columns should be returned (defaults to
#'   `columns = colnames(object)`), column variable names of the target
#'   object need to be prefixed with `target_`. The score from the matching function is
#'   returned as variable `"score"`. Similar to `$`, this function
#'   performs a *left join* of variables from the *query* and *target*
#'   objects returning all values for all query elements (eventually returning
#'   duplicated elements for query elements matching multiple target elements)
#'   and the values for the target elements matched to the respective query
#'   elements. See help on `$` above or examples below for details.
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
#' @param columns for `data`: `character` vector with column names of variables
#'   that should be extracted.
#' 
#' @param drop for `[`: ignored.
#'
#' @param i `integer` or `logical` defining the `query` elements to keep.
#'
#' @param j for `[`: ignored.
#' 
#' @param matches `data.frame` with columns `"query_idx"` (`integer`),
#'   `"target_idx"` (`integer`) and `"score"` (`numeric`) representing the
#'   *n:m* mapping of elements between the `query` and the `target` objects.
#'
#' @param name for `$`: the name of the columns to extract.
#' 
#' @param object `Matched` object.
#'
#' @param target object with the elements against which `query` has been
#'   matched.
#'
#' @param query object with the query elements.
#'
#' @param x `Matched` object.
#'
#' @return See individual method description above for details. 
#' 
#' @exportClass Matched
#'
#' @author Johannes Rainer
#' 
#' @rdname Matched
#'
#' @examples
#'
#' ## Creating a Matched object.
#' q1 <- data.frame(col1 = 1:5, col2 = 6:10)
#' t1 <- data.frame(col1 = 11:16, col2 = 17:22)
#' ## Define matches between query row 1 with target row 2 and, query row 2 
#' with target rows 2,3,4 and query row 5 with target row 5.
#' mo <- Matched(
#'     q1, t1, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
#'                                 target_idx = c(2L, 2L, 3L, 4L, 5L),
#'                                  score = seq(0.5, 0.9, by = 0.1)))
#' 
#' ## Which of the query elements (rows) match at least one target element (row)?
#' whichQuery(mo)
#'
#' 
#' ## Which target elements (rows) match at least one query element (row)?
#' whichTarget(mo)
#' 
#' ## Extracting variable "col1" from query object .
#' mo$col1
#'
#' ## We have duplicated values for the entires of `col1` related to query elements
#' (rows) matched to multiple rows of the target object). The value of `col1` 
#' is returned for each element (row) in the query spectrum.
#'
#' ## Extracting variable "col1" from target object. Note that only values of
#' ## `col1` for rows matching at least one query row are returned
#' ## and an NA is reported for query rows without matching target rows.
#' mo$target_col1
#'
#' ## The 3rd and 4th query rows do not match any target row, thus `NA` is returned.
#'
#' ## `data` can be used to extract all (or selected) columns
#' ## from the object. Same as with `$`, a left join between the columns
#' ## from the query and the target is performed. The
#' ## prefix `"target_"` is used to label the columns from the target. 
#' Below we extract selected columns from the object as a DataFrame.
#' res <- data(mo, columns = c("col1", "col2",
#'     "target_col1", "target_col2"))
#' res
#' res$col1
#' res$target_col1
#'
#' ## Again, all values for query columns are returned and for query elements not
#' ## matching any target element NA is reported as value for the respective
#' ## variable.
#'
#' ## The example matched object contains all query and all target
#' ## elements (rows). Below we subset the object keeping only query rows that are
#' ## matched to at least one target row.
#' mo_sub <- mo[whichQuery(mo)]
#'
#' ## ms_sub contains now only 3 query rows:
#' nrow(query(ms_sub))
#'
#' ## while the original object contains all 5 query rows:
#' nrow(query(mo))
#'
#' ## Both object contain however still the full target object:
#' nrow(target(mo))
#' nrow(target(ms_sub))
#'
#' ## With the `pruneTarget` we can however reduce also the target rows to
#' ## only those that match at least one query row
#' ms_sub <- pruneTarget(ms_sub)
#' nrow(target(ms_sub))
NULL


setClass(
    "Matched",
    slots = c(
        query = "ANY",
        target = "ANY",
        matches = "data.frame",
        ## metadata
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
    msg <- .validate_matches_content(object@matches, 
                                     ifelse(is.null(d <- dim(object@query)),
                                            length(object@query), d[1]),
                                     ifelse(is.null(d <- dim(object@target)),
                                            length(object@target), d[1]))
    if (length(msg)) return(msg)
    # if(is(object@query)[1] != is(object@target)[1])
    #     return("query and target are not of the same type")
    TRUE
})



#' @export
#'
#' @rdname Matched
setMethod("length", "Matched", function(x) ifelse(is.null(d <- dim(x@query)),
                                                  length(x@query), d[1])) 

#' @exportMethod show
#'
#' @importMethodsFrom methods show
#'
#' @rdname Matched
setMethod("show", "Matched", function(object) {
    cat("Object of class", class(object)[1L], "\n")
    cat("Total number of matches:", nrow(object@matches), "\n")
    cat("Number of query objects: ",
        ifelse(is.null(d <- dim(object@query)), length(object@query), d[1]), 
        " (", length(unique(object@matches$query_idx)), " matched)\n", sep = "")
    cat("Number of target objects: ",
        ifelse(is.null(d <- dim(object@target)), length(object@target), d[1]), 
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
    if(name %in% colnames(x))
    {
        idxs <- .fill_index(seq_len(length(x)), x@matches$query_idx)
        if (name %in% colnames(x@matches)) {
            keep <- idxs %in% x@matches$query_idx
            idxs[keep] <- seq_len(sum(keep))
            idxs[!keep] <- NA
            return(x@matches[idxs, name])
        }
        if (length(grep("^target_", name))) {
            vals <- do.call("$", list(x@target[x@matches$target_idx, ],
                                      sub("target_", "", name)))
            keep <- idxs %in% x@matches$query_idx
            idxs[keep] <- seq_len(sum(keep))
            idxs[!keep] <- NA
            vals[idxs]
        } else
            do.call("$", list(x@query, name))[idxs] 
    }
})

#' @exportMethod colnames
#' 
#' @rdname Matched
setMethod("colnames", "Matched", function(x) {
    cns <- colnames(x@matches)
    cnq <- NULL
    cnt <- NULL
    if(is.data.frame(target(x)))  cnt <- paste0("target_", colnames(target(x)))
    if(is.data.frame(query(x))) cnq <- colnames(query(x))
    c(cnq, cnt, cns[!cns %in% c("query_idx", "target_idx")])
})

# I'm not sure if I should implement this with SetMethod and if I should change 
# DataFrame to data.frame

#' @importMethodsFrom S4Vectors cbind
#'
#' @importFrom S4Vectors DataFrame
#' 
#' @rdname Matched
#'
#' @export
#' 
data <- function(object, columns = colnames(object)){
    if (any(!columns %in% colnames(object)))
        stop("column(s) ", paste0(columns[!columns %in% colnames(object)],
                                  collapse = ", "), " not available")
    if(length(d <- dim(object@query)) == 2) query_idxs <- seq_len(d[1])
    else query_idxs <- seq_along(object@query)
    idxs <- .fill_index(query_idxs, object@matches$query_idx)
    from_target <- grepl("^target_", columns)
    from_matches <- columns %in% colnames(object@matches)
    from_query <- !(from_target | from_matches)
    res <- DataFrame()
    if (any(from_query))
        res <- object@query[idxs, from_query, drop = FALSE]
    if (any(from_target)) {
        spt <- object@target[object@matches$target_idx,
                             sub("target_", "", columns[from_target]), drop = FALSE]
        colnames(spt) <- paste0("target_", colnames(spt))
        keep <- idxs %in% object@matches$query_idx
        target_idxs <- idxs
        target_idxs[keep] <- seq_len(sum(keep))
        target_idxs[!keep] <- NA
        if (nrow(res) == length(idxs))
            res <- cbind(res, spt[target_idxs, , drop = FALSE])
        else res <- spt[target_idxs, , drop = FALSE]
    }
    if (any(from_matches)) {
        keep <- idxs %in% object@matches$query_idx
        idxs[keep] <- seq_len(sum(keep))
        idxs[!keep] <- NA
        if (nrow(res) == length(idxs))
            res <- cbind(res,
                         DataFrame(object@matches[idxs,
                                                  columns[from_matches],
                                                  drop = FALSE]))
        else res <- DataFrame(object@matches[idxs,
                                             columns[from_matches],
                                             drop = FALSE])
    }
    res[, columns, drop = FALSE]  
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
    if(is.null(dim(x@query))) slot(x, "query", check = FALSE) <- x@query[i]
    else slot(x, "query", check = FALSE) <- x@query[i, ]
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

#' @rdname Matched
#'
#' @importFrom methods validObject
#' 
#' @export
pruneTarget <- function(object) {
    keep <- whichTarget(object)
    if (is.null(dim(object@target))) object@target <- object@target[keep]
    else object@target <- object@target[keep, ]
    object@matches$target_idx <- match(object@matches$target_idx, keep)
    validObject(object)
    object
}
