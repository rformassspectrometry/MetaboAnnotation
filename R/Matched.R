#' @title Representation of generic objects matches
#'
#' @name Matched
#'
#' @aliases Matched Matched-class [,Matched-method
#'
#' @description
#'
#' Matches between *query* and *target* generic objects can be represented by
#' the `Matched` object. By default, all data accessors work as
#' *left joins* between the *query* and the *target* object, i.e. values are
#' returned for each *query* object with eventual duplicated entries (values)
#' if the *query* object matches more than one *target* object. See also
#' *Creation and subsetting* as well as *Extracting data* sections below for
#' details and more information.
#'
#' @section Creation and subsetting:
#'
#' `Matched` object is returned as result from the [matchValues()] function.
#'
#' Alternatively, `Matched` objects can also be created with the `Matched`
#' function providing the `query` and `target` objects as well as the `matches`
#' `data.frame` with two columns of integer indices defining which elements
#' from *query* match which element from *target*.
#'
#' - `[` subset the object selecting `query` object elements to keep with
#'   parameter `i`. The resulting object will contain all the matches
#'   for the selected query elements. The `target` object will by default be
#'   returned as-is.
#'
#' - `addMatches`: add new matches to an existing object. Parameters
#'   `queryValue` and `targetValue` allow to define which element(s) in
#'   `query` and `target` should be considered matching. If `isIndex = TRUE`,
#'   both `queryValue` and `targetValue` are considered to be integer indices
#'   identifying the matching elements in `query` and `target`, respectively.
#'   Alternatively (with `isIndex = FALSE`) `queryValue` and `targetValue` can
#'   be elements in columns `queryColname` or `targetColname` which can be used
#'   to identify the matching elements. Note that in this case
#'   **only the first** matching pair is added. Parameter `score` allows to
#'   provide the score for the match. It can be a numeric with the score or a
#'   `data.frame` with additional information on the manually added matches. In
#'   both cases its length (or number of rows) has to match the length of
#'   `queryValue`. See examples below for more information.
#'
#' - `filterMatches`: keeps or removes matches corresponding to certain indexes
#'   or values of `query` and `target`. If `queryValue` and `targetValue` are
#'   provided, matches for these value pairs are kept or removed. Parameter
#'   `index` allows to filter matches providing their index in the [matches()]
#'   matrix. Note that `filterMatches` removes only matches from the [matches()]
#'   matrix from the `Matched` object but thus not alter the `query` or `target`
#'   in the object. See examples below for more information.
#'
#' - `pruneTarget` *cleans* the object by removing non-matched
#'   **target** elements.
#'
#'
#' @section Extracting data:
#'
#' - `$` extracts a single variable from the `Matched` `x`. The variables that
#'   can be extracted can be listed using `colnames(x)`. These variables can
#'   belong to *query*, *target* or be related to the matches (e.g. the
#'   score of each match). If the *query* (*target*) object is two dimensional,
#'   its columns can be extracted (prefix` "target_"` is used for columns in the
#'   *target* object) otherwise if *query* (*target*) has only a single
#'   dimension (e.g. is a `list` or a `character`) the whole object can be
#'   extracted with `x$query` (`x$target`). More precisely, when
#'   *query* (*target*) is a `SummarizedExperiment` the columns from
#'   `rowData(query)` (rowData(`target`)) are extracted; when *query* (*target*)
#'   is a `QFeatures` the columns from `rowData` of the assay specified in the
#'   `queryAssay` (`targetAssay`) slot are extracted. The matching scores
#'   are available as *variable* `"score"`. Similar to a left join between the
#'   query and target elements, this function returns a value for each query
#'   element, with eventual duplicated values for query elements matching more
#'   than one target element. If variables from the target `data.frame` are
#'   extracted, an `NA` is reported for the entries corresponding to *query*
#'   elements that don't match any target element. See examples below for
#'   more details.
#'
#' - `length` returns the number of **query** elements.
#'
#' - `matchedData` allows to extract multiple variables contained in the
#'   `Matched` object as a `DataFrame`. Parameter `columns` allows to
#'   define which columns (or variables) should be returned (defaults to
#'   `columns = colnames(object)`). Each single column in the returned
#'   `DataFrame` is constructed in the same way as in `$`. That is, like `$`,
#'   this function performs a *left join* of variables from the *query* and
#'   *target* objects returning all values for all *query* elements
#'   (eventually returning duplicated elements for query elements matching
#'   multiple target elements) and the values for the target elements matched
#'   to the respective query elements (or `NA` if the target element is not
#'   matched to any query element).
#'
#' - `matches` returns a `data.frame` with the actual matching information with
#'   columns `"query_idx"` (index of the element in `query`), `"target_idx"`
#'   (index of the element in `target`) `"score"` (the score of the match) and
#'   eventual additional columns.
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
#' @param index for `filterMatches`: indices of the matches to keep (if
#'  `keep = TRUE`) or to drop if (`keep = FALSE`).
#'
#' @param isIndex for `addMatches`: specifies if `queryValue` and
#' `targetValue` are expected to be vectors of indexes.
#'
#' @param j for `[`: ignored.
#'
#' @param keep for `filterMatches`: `logical`. If `keep = TRUE` the matches are
#' kept, if `keep = FALSE` they are removed.
#'
#' @param matches `data.frame` with columns `"query_idx"` (`integer`),
#'   `"target_idx"` (`integer`) and `"score"` (`numeric`) representing the n:m
#'   mapping of elements between the `query` and the `target` objects.
#'
#' @param metadata `list` with optional additional metadata.
#'
#' @param name for `$`: the name of the column (or variable) to extract.
#'
#' @param object a `Matched` object.
#'
#' @param score for `addMatches`: `numeric` (same length than `queryValue`) or
#'   `data.frame` (same number of rows than `queryValue`) specifying the scores
#'   for the matches to add. If not specified, a `NA` will be used as score.
#'
#' @param target object with the elements against which `query` has been
#'   matched.
#'
#' @param targetAssay `character` that needs to be specified when `target` is
#'   `QFeatures` and corresponds to the name of one of the assay in `target`.
#'   In this case, the `Matched` object represents the matches between data in
#'   `query` and the `rowData` of that assay.
#'   
#' @param targetColname if `query` is 2-dimensional: column of `target` against
#'   which elements of `targetValue` are compared.
#'
#' @param targetValue for `filterMatches`: vector of values to search for in
#'   `target` (if `target` is 1-dimensional) or in column `targetColname` of
#'   `target` (if `target` is 2-dimensional). For `addMatches`: either an
#'   index in `target` or value in column `targetColname` of `target` defining
#'   (together with `queryValue`) the pair of query and target elements for
#'   which a match should be manually added. Lengths of `queryValue` and
#'   `targetValue` have to match.
#'
#' @param query object with the query elements.
#' 
#' @param queryAssay `character` that needs to be specified when `query` is
#'   `QFeatures` and corresponds to the name of one of the assay in `query`.
#'   In this case, the `Matched` object represents the matches between
#'   the `rowData` of that assay and data in `target`.
#'
#' @param queryColname if `query` is 2-dimensional: column of `query` against
#'   which elements of `queryValue` are compared.
#'
#' @param queryValue for `filterMatches`: vector of values to search for in
#'   `query` (if `query` is 1-dimensional) or in column `queryColname` of
#'   `query` (if `query` is 2-dimensional). For `addMatches`: either an index
#'   in `query` or value in column `queryColname` of `query` defining (together
#'   with `targetValue`) the pair of query and target elements for which a
#'   match should be manually added. Lengths of `queryValue` and
#'   `targetValue` have to match.
#'
#' @param x `Matched` object.
#'
#' @param ... additional parameters.
#'
#' @return See individual method description above for details.
#'
#' @seealso [MatchedSpectra()] for matched [Spectra()] objects.
#'
#' @exportClass Matched
#'
#' @author Andrea Vicini, Johannes Rainer
#'
#' @rdname Matched
#'
#' @examples
#'
#' ## Creating a `Matched` object.
#' q1 <- data.frame(col1 = 1:5, col2 = 6:10)
#' t1 <- data.frame(col1 = 11:16, col2 = 17:22)
#' ## Define matches between query row 1 with target row 2 and, query row 2
#' ## with target rows 2,3,4 and query row 5 with target row 5.
#' mo <- Matched(
#'     q1, t1, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
#'                                  target_idx = c(2L, 2L, 3L, 4L, 5L),
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
#' ########
#' ## Creating a `Matched` object with a `data.frame` for `query` and a `vector`
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
#' ## Reducing the target elements to only those that match at least one query
#' ## row
#' mo_sub <- pruneTarget(mo_sub)
#' length(target(mo_sub))
#'
#' ########
#' ## Filtering `Matched` with `filterMatches`
#'
#' ## Inspecting the matches in `mo`:
#' mo$col1
#' mo$target
#'
#' ## We have thus target *12* matched to both query elements with values 1 and
#' ## 2, and query element 2 is matching 3 target elements. Let's assume we want
#' ## to resolve this multiple mappings to keep from them only the match between
#' ## query 1 (column `"col1"` containing value `1`) with target 1 (value `12`)
#' ## and query 2 (column `"col1"` containing value `2`) with target 2 (value
#' ## `13`). In addition we also want to keep query element 5 (value `5` in
#' ## column `"col1"`) with the target with value `15`:
#' mo_sub <- filterMatches(mo, queryValue = c(1, 2, 5), queryColname = "col1",
#'     targetValue = c(12, 13, 15))
#' matchedData(mo_sub)
#'
#' ## Alternatively to specifying the matches to filter with `queryValue` and
#' ## `targetValue` it is also possible to specify directly the index of the
#' ## match(es) in the `matches` `data.frame`:
#' matches(mo)
#'
#' ## To keep only matches like in the example above we could use:
#' mo_sub <- filterMatches(mo, index = c(1, 3, 5))
#' matchedData(mo_sub)
#'
#' ## Note also that, instead of keeping the specified matches, it would be
#' ## possible to remove them by setting `keep = FALSE`. Below we remove
#' ## selected matches from the object:
#' mo_sub <- filterMatches(mo, queryValue = c(2, 2), queryColname = "col1",
#'     targetValue = c(12, 14), keep = FALSE)
#' mo_sub$col1
#' mo_sub$target
#'
#' ########
#' ## Adding matches using `addMatches`
#'
#' ## `addMatches` allows to manually add matches. Below we add a new match
#' ## between the `query` element with a value of `1` in column `"col1"` and
#' ## the target element with a value of `15`. Parameter `score` allows to
#' ## assign a score value to the match.
#' mo_add <- addMatches(mo, queryValue = 1, queryColname = "col1",
#'     targetValue = 15, score = 1.40)
#' matchedData(mo_add)
#' ## Matches are always sorted by `query`, thus, the new match is listed as
#' ## second match.
#'
#' ## Alternatively, we can also provide a `data.frame` with parameter `score`
#' ## which enables us to add additional information to the added match. Below
#' ## we define the score and an additional column specifying that this match
#' ## was added manually. This information will then also be available in the
#' ## `matchedData`.
#' mo_add <- addMatches(mo, queryValue = 1, queryColname = "col1",
#'     targetValue = 15, score = data.frame(score = 5, manual = TRUE))
#' matchedData(mo_add)
#'
#' ## The match will get a score of NA if we're not providing any score.
#' mo_add <- addMatches(mo, queryValue = 1, queryColname = "col1",
#'     targetValue = 15)
#' matchedData(mo_add)
#'
#' ## Creating a `Matched` object with a `SummarizedExperiment` for `query` and
#' ## a `vector` for `target`. The matches are specified in the same way as
#' ## the example before.
#' library(SummarizedExperiment)
#' q1 <- SummarizedExperiment(
#'   assays = data.frame(matrix(NA, 5, 2)),
#'   rowData = data.frame(col1 = 1:5, col2 = 6:10),
#'   colData = data.frame(cD1 = c(NA, NA), cD2 = c(NA, NA)))
#' t1 <- data.frame(col1 = 11:16, col2 = 17:22)
#' ## Define matches between row 1 in rowData(q1) with target row 2 and,
#' ## rowData(q1) row 2 with target rows 2,3,4 and rowData(q1) row 5 with target
#' ## row 5.
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
#' ## Extracting variable "col1" from rowData(q1).
#' mo$col1
#'
#' ## We have duplicated values for the entries of `col1` related to rows of
#' ## rowData(q1) matched to multiple rows of the target data.frame t1. The
#' ## value of `col1` is returned for each row in the rowData of query.
#'
#' ## Extracting variable "col1" from target object. To access columns from
#' ## target we have to prefix the name of the column by `"target_"`.
#' ## Note that only values of `col1` for rows matching at least one row in
#' ## rowData of query are returned and an NA is reported for those without
#' ## matching target rows.
#' mo$target_col1
#'
#' ## The 3rd and 4th query rows do not match any target row, thus `NA` is
#' ## returned.
#'
#' ## `matchedData` can be used to extract all (or selected) columns
#' ## from the object. Same as with `$`, a left join between the columns
#' ## from the query and the target is performed. Below we extract selected
#' ## columns from the object as a DataFrame.
#' res <- matchedData(mo, columns = c("col1", "col2", "target_col1",
#'                                   "target_col2"))
#' res
#' res$col1
#' res$target_col1
#'
#' ## The example `Matched` object contains all rows in the
#' ## `rowData` of the `SummarizedExperiment` and all target rows. Below we
#' ## subset the object keeping only rows that are matched to at least one
#' ## target row.
#' mo_sub <- mo[whichQuery(mo)]
#'
#' ## mo_sub contains now a `SummarizedExperiment` with only 3 rows:
#' nrow(query(mo_sub))
#'
#' ## while the original object contains a `SummarizedExperiment` with all 5
#' ## rows:
#' nrow(query(mo))
#'
#' ## Both objects contain however still the full target object:
#' nrow(target(mo))
#' nrow(target(mo_sub))
#'
#' ## With the `pruneTarget` we can however reduce also the target rows to
#' ## only those that match at least one in the `rowData` of query
#' mo_sub <- pruneTarget(mo_sub)
#' nrow(target(mo_sub))
NULL

setClass(
    "Matched",
    slots = c(
        query = "ANY",
        target = "ANY",
        matches = "data.frame",
        queryAssay = "character",
        targetAssay = "character",
        metadata = "list",
        version = "character"
    ),
    prototype = prototype(
        query = list(),
        target = list(),
        matches = data.frame(query_idx = integer(),
                             target_idx = integer(),
                             score = numeric()),
        queryAssay = character(),
        targetAssay = character(),
        metadata = list(),
        version = "0.1")
)

#' @export
#'
#' @rdname Matched
Matched <- function(query = list(), target = list(),
                    matches = data.frame(query_idx = integer(),
                                         target_idx = integer(),
                                         score = numeric()),
                    queryAssay = character(), targetAssay = character(),
                    metadata = list()) {
    new("Matched", query = query, target = target, matches = matches,
        queryAssay = queryAssay, targetAssay = targetAssay, metadata = metadata)
}

setValidity("Matched", function(object) {
    msg <- .validate_matches_format(object@matches)
    if (length(msg)) return(msg)
    msg <- .validate_matches_content(object@matches, length(object),
                                     .nelements_t(object))
    if(any(c("query", "target") %in% colnames(object@matches)))
        return("\"query\" and \"target\" can't be used as matches column names")
    if (length(msg)) return(msg)
    msg <- .validate_qt(object@query)
    if (length(msg)) return(msg)
    msg <- .validate_qt(object@target)
    if (length(msg)) return(msg)
    msg <- .validate_assay(object@query, object@queryAssay)
    if (length(msg)) return(msg)
    msg <- .validate_assay(object@target, object@targetAssay, "Target")
    if (length(msg)) return(msg)
    TRUE
})

.validate_assay <- function(x, assay, what = "Query") {
    if (inherits(x, "QFeatures")) {
        if (length(assay) != 1)
            return(paste0("\"","assay", what,  "\" must be of length 1"))
        if (!assay %in% names(x)) 
            return(paste0("No assay \"", assay, "\" in \"",tolower(what), "\""))
    }
}

#' @export
#'
#' @rdname Matched
setMethod("length", "Matched",
          function(x) .nelements(.objectToMatch(x@query, x@queryAssay)))

.nelements_t <- function(x) .nelements(.objectToMatch(x@target, x@targetAssay))

#' @exportMethod show
#'
#' @importMethodsFrom methods show
#'
#' @rdname Matched
setMethod("show", "Matched", function(object) {
    cat("Object of class", class(object)[1L], "\n")
    cat("Total number of matches:", nrow(object@matches), "\n")
    cat("Number of query objects: ", length(object),
        " (", length(unique(object@matches$query_idx)), " matched)\n", sep = "")
    cat("Number of target objects: ", .nelements_t(object), " (",
        length(unique(object@matches$target_idx)), " matched)\n", sep = "")
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
matches <- function(object) {
    object@matches
}

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
  .dollar(.objectToMatch(x@query, x@queryAssay),
          .objectToMatch(x@target, x@targetAssay), x@matches, name)
})

#' @importFrom BiocGenerics colnames
#'
#' @exportMethod colnames
#'
#' @rdname Matched
setMethod("colnames", "Matched", function(x) {
  .colnames(.objectToMatch(x@query, x@queryAssay),
            .objectToMatch(x@target, x@targetAssay), x@matches)
})

#' @importMethodsFrom S4Vectors cbind
#'
#' @importFrom S4Vectors DataFrame
#'
#' @rdname Matched
#'
#' @export
setMethod("matchedData", "Matched", function(object,
                                             columns = colnames(object), ...) {
    .matchedData(.objectToMatch(object@query, object@queryAssay),
                 .objectToMatch(object@target, object@targetAssay),
                 object@matches, columns, ... )
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
    if (is(x, "list") && any(na <- is.na(i)))
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
    slot(x, "query", check = FALSE) <- .subset_qt(x@query, x@queryAssay, i)
    mtches <- x@matches[x@matches$query_idx %in% i, , drop = FALSE]
    if (nrow(mtches)) {
        ## Support handling duplicated indices.
        mtches <- split.data.frame(
            mtches, f = as.factor(mtches$query_idx))[as.character(i)]
        lns <- vapply(mtches, function(z)
            if (length(z)) nrow(z) else 0L, integer(1))
        mtches <- do.call(rbind, mtches[lengths(mtches) > 0])
        rownames(mtches) <- NULL
        mtches$query_idx <- rep(seq_along(i), lns)
    }
    slot(x, "matches", check = FALSE) <- mtches
    x
}

.subset_qt <- function(x, assay, i, j, drop = FALSE) {
    if (is(x, "QFeatures")) {
        if (!assay %in% names(x))
            stop("Invalid assay name.")
        x[[assay]] <- .extract_elements(x[[assay]], i, j, drop)
    } else {
        x <- .extract_elements(x, i, j, drop)
    }
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
        msg <- c(msg, "indices in \"query_idx\" are out-of-bounds")
    if (!all(x$target_idx %in% seq_len(ntarget)))
        msg <- c(msg, "indices in \"target_idx\" are out-of-bounds")
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

.cnt <- function(target) {
    ndim <- length(dim(target))
    if (ndim == 2) {
        cnt <- colnames(target)
        if (length(cnt)) cnt <- paste0("target_", cnt)
    } else if (ndim == 0) {
        cnt <- "target"
    } else {
        stop("unsupported dimensions in \"target\"")
    }
    cnt
}

.colnames <- function(query, target, matches) {
  cns <- colnames(matches)
  cnq <- NULL
  cnt <- .cnt(target)
  if (length(dim(query)) == 2)
    cnq <- colnames(query)
  if (is.null(dim(query)))
    cnq <- "query"
  c(cnq, cnt, cns[!cns %in% c("query_idx", "target_idx")])
}

.dollar <- function(query, target, matches, name) {
  if(name %in% .colnames(query, target, matches))
  {
    not_mtchd <- setdiff(seq_len(.nelements(query)), matches$query_idx)
    idxs_qry <- c(matches$query_idx, not_mtchd)
    ord <- order(idxs_qry)
    if (name %in% colnames(matches)) {
      idxs_mtch <- c(seq_len(nrow(matches)), rep(NA, length(not_mtchd)))[ord]
      return(matches[idxs_mtch, name])
    }
    if (name %in% .cnt(target)) {
      idxs_trg <- c(matches$target_idx, rep(NA, length(not_mtchd)))[ord]
      .extract_elements(target, idxs_trg, sub("target_", "", name), drop = TRUE)
    }else
      .extract_elements(query, idxs_qry[ord], name, drop = TRUE)
  } else stop("'", name, "' not available")
}

.matchedData <- function(query, target, matches, columns, ...) {
  cnms <- .colnames(query, target, matches)
  if (any(!columns %in% cnms))
    stop("column(s) ", paste0(columns[!columns %in% cnms],
                              collapse = ", "), " not available")
  not_mtchd <- setdiff(seq_len(.nelements(query)), matches$query_idx)
  idxs_qry <- c(matches$query_idx, not_mtchd)
  ord <- order(idxs_qry)
  from_target <- columns %in% .cnt(target)
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

#' @importMethodsFrom SummarizedExperiment rowData
#' 
#' @noRd
.objectToMatch <- function(x, assayname = character()) {
    #what <- as.character(sys.call()[-1])[1]
    if (is(x, "QFeatures")) {
        if (!(len_ass <- length(assayname)))
            stop ("assay name has to be provided when x is `QFeatures`.")
        if (len_ass != 1)
            stop ("assay name must be `character(1)`")
        if (!assayname %in% names(x))
            stop ("No assay in x with name \"", assayname, "\"")
        else x <- x[[assayname]] 
    }
    if(is(x, "SummarizedExperiment"))
        x <- rowData(x)
    x
}

#' @rdname Matched
#'
#' @importFrom methods validObject
#'
#' @export
pruneTarget <- function(object) {
    keep <- whichTarget(object)
    object@target <- .subset_qt(object@target, object@targetAssay, keep)
    object@matches$target_idx <- match(object@matches$target_idx, keep)
    validObject(object)
    object
}

.findMatchesIdxs <- function(query, target, matches, queryValue = integer(),
                         targetValue = integer(),
                         queryColname = character(),
                         targetColname = character()) {
    if (length(queryValue) != length(targetValue))
        stop("'queryValue' and 'targetValue' must have the same length")
    if (length(dim(query)) == 2) {
        if (length(queryColname) == 0)
            stop("\"", queryColname,"\" must be set when 'query' is 2-dimensional")
        if (!queryColname %in% colnames(query))
            stop("\"", queryColname,"\" is not a column of 'query'")
    }
    targetColname <- sub("target_", "", targetColname)
    if (length(dim(target)) == 2) {
        if(length(targetColname) == 0)
            stop("\"", targetColname,"\" must be set when 'target' is 2-dimensional")
        if (!targetColname %in% colnames(target))
            stop("\"", targetColname,"\" is not a column of 'target'")
    }
    mq <- .extract_elements(query, matches$query_idx, queryColname,
                            drop = TRUE)
    mt <- .extract_elements(target, matches$target_idx, targetColname,
                                               drop = TRUE)
    unlist(sapply(seq_along(queryValue), function(i)
        which(mq == queryValue[i] & mt == targetValue[i])))
}

#' @rdname Matched
#'
#' @importFrom methods validObject
#'
#' @export
setMethod("filterMatches", "Matched", function (object, queryValue = integer(),
                                                targetValue = integer(),
                                                queryColname = character(),
                                                targetColname = character(),
                                                index = integer(),
                                                keep = TRUE, ...) {
    if (length(index) && any(!index %in% seq_len(nrow(object@matches))))
        stop("some indexes in \"index\" are out-of-bounds")
    if (!length(index) && length(queryValue))
        index  <- .findMatchesIdxs(.objectToMatch(object@query, object@queryAssay),
                                   .objectToMatch(object@target, object@targetAssay),
                                   object@matches, queryValue,
                                   targetValue, queryColname,
                                   targetColname)
    if (keep) to_keep <- seq_len(nrow(object@matches)) %in% index
    else to_keep <- !seq_len(nrow(object@matches)) %in% index
    object@matches <- object@matches[to_keep, , drop = FALSE]
    validObject(object)
    object
})

#' @importFrom MsCoreUtils rbindFill
.addMatches <- function(query, target, matches, queryValue = integer(),
                        targetValue = integer(), queryColname = character(),
                        targetColname = character(),
                        score = rep(NA_real_, length(queryValue)),
                        isIndex = FALSE) {
    if (!is.data.frame(score))
        score <- data.frame(score = score)
    if (!is.data.frame(score))
        stop("'score' needs to be either a 'data.frame' or a numeric")
    if (length(queryValue) != length(targetValue) ||
        length(queryValue) != nrow(score))
        stop("'queryValue', 'targetValue' and 'score' must have the same length")
    if (isIndex) {
        if (!is.integer(queryValue) || !is.integer(targetValue))
            stop(paste0("'queryValue' and 'targetValue' must be integer ",
                        "vectors when 'isIndex = TRUE'"))
        if (any(!queryValue %in% seq_len(.nelements(query))))
            stop("Provided indices in 'queryValue' are out-of-bounds.")
        if (any(!targetValue %in% seq_len(.nelements(target))))
            stop("Provided indices in 'queryValue' are out-of-bounds.")
        query_idx <- queryValue
        target_idx <- targetValue
    } else {
        if (length(dim(query)) == 2)
            if (!queryColname %in% colnames(query))
                stop("\"", queryColname,"\" is not a column of 'query'")
        if (length(dim(target)) == 2)
            if (!targetColname %in% colnames(target))
                stop("\"", targetColname,"\" is not a column of 'target'")
        mq <- match(queryValue, .extract_elements(query, j = queryColname,
                                                  drop = TRUE))
        mt <- match(targetValue, .extract_elements(target, j = targetColname,
                                                   drop = TRUE))
        to_keep <- !is.na(mq) & !is.na(mt)
        query_idx <- mq[to_keep]
        target_idx <- mt[to_keep]
        score <- score[to_keep, ]
    }
    new_matches <- cbind(
        data.frame(query_idx = query_idx, target_idx = target_idx),
        score)
    new_matches <- rbindFill(matches, new_matches)
    ## remove possible matches that were already in matches
    new_matches[!duplicated(new_matches[, c("query_idx", "target_idx")]), ]
}

#' @rdname Matched
#'
#' @importFrom methods validObject
#'
#' @export
setMethod("addMatches", "Matched",
          function(object, queryValue = integer(), targetValue = integer(),
                   queryColname = character(), targetColname = character(),
                   score = rep(NA_real_, length(queryValue)), isIndex = FALSE) {
              object@matches <- .addMatches(.objectToMatch(object@query, object@queryAssay),
                                            .objectToMatch(object@target, object@targetAssay),
                                            object@matches, queryValue,
                                            targetValue, queryColname,
                                            targetColname, score, isIndex)
              validObject(object)
              object
          })
