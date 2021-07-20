#' @title Representation of SummarizedExperiment matches
#'
#' @name MatchedSummarizedExperiment
#'
#' @aliases MatchedSummarizedExperiment MatchedSummarizedExperiment-class
#'
#' @description
#'
#' Matches between a query `SummarizedExperiment` and target `data.frame` can be 
#' represented by the`MatchedSummarizedExperiment` object. Functions like the 
#' [matchMz()] function can return this type of object. 
#' By default, all data accessors work as *left joins* between the *query* 
#' `SummarizedExperiment` and the *target* `data.frame`, i.e. values are 
#' returned for each row in `rowData` of *query*  with eventual duplicated 
#' entries (values) if such row matches more than one target row.
#'
#' @section Creation and subsetting:
#'
#' `MatchedSummarizedExperiment` objects can be created with the 
#' `MatchedSummarizedExperiment` function providing the `query` 
#' `SummarizedExperiment` and `target` `data.frame` as well as a `data.frame`
#'  with two columns of integer indices defining which elements from `rowData` 
#'  of *query* match which element from *target*.
#'
#' - `[` subsets the rows to keep (specified throught parameter `i`) in 
#'   `rowData` of `query` . The `target` `data.frame` will by default be 
#'   returned as-is.
#'
#' - `pruneTarget` *cleans* the `MatchedSummarizedExperiment` object by removing 
#'   non-matched rows in target `data.frame`.
#'
#' @section Extracting data:
#'
#' - `$` extracts a single variable from the `MatchedSummarizedExperiment` `x`. 
#'   The variables that can be extracted can be listed using `colnames(x)`. 
#'   These variables can belong to `rowData` of *query*, *target* or be related 
#'   to the matches (e.g. the score of each match available as *variable* 
#'   `"score"`). Prefix` "target_"` is used for columns in the *target* object.
#'   Similar to a left join between the query `rowData` rows and target rows,
#'   this function returns a value for each query row, with eventual
#'   duplicated values for query rows matching more than one target element.
#'   If variables from the target `data.frame` are extracted, an `NA` is
#'   reported for the entries corresponding to *query* `rowData` rows that don't 
#'   match any target row See examples below for more details.
#'
#' - `length` returns the number of rows in **query** `rowData`.
#'
#' - `matchedData` allows to extract multiple variables contained in the
#'   `MatchedSummarizedExperiment` object as a `DataFrame`. Parameter `columns` 
#'   allows to define which columns (or variables) should be returned (defaults 
#'   to `columns = colnames(object)`). Each single column in the returned
#'   `DataFrame` is constructed in the same way as in `$`. That is, like `$`,
#'   this function performs a *left join* of variables from the *query* and
#'   *target* objects returning all values for rows in query `rowData` 
#'   (eventually returning duplicated elements for query elements matching 
#'   multiple targetelements) and the values for the target elements matched to 
#'   the respective query elements (or `NA` if the target element is not matched 
#'   to any query element).
#'
#' - `target` returns the *target* `data.frame`.
#'
#' - `query` returns the *query* `SummarizedExperiment`.
#'
#' - `whichTarget` returns an `integer` with the indices of the rows in
#'   *target* that match at least one row in `rowData` of *query*.
#'
#' - `whichQuery` returns an `integer` with the indices of the rows in `rowData` 
#'   of *query* that match at least one row in *target*.
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
#'   mapping between the rows of `rowData` of query` and the `target` rows.
#'
#' @param name for `$`: the name of the column (or variable) to extract.
#'
#' @param object a `MatchedSummarizedExperiment` object.
#'
#' @param target `data.frame` with the elements against which `query` has been
#'   matched.
#'
#' @param query `SummarizedExperiment` object having as `rowData` a `DFrame` 
#'   with rows to be matched to those in `target`
#'
#' @param x `MatchedSummarizedExperiment` object.
#'
#' @param ... additional parameters.
#'
#' @return See individual method description above for details.
#' 
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
#' @exportClass MatchedSummarizedExperiment
#'
#' @author Andrea Vicini, Johannes Rainer
#'
#' @rdname MatchedSummarizedExperiment
#'
#' @examples
#'
#' ## Creating a MatchedSummarizedExperiment object.
#' library(SummarizedExperiment)
#' q1 <- SummarizedExperiment(
#'   assays = data.frame(matrix(NA, 5, 2)), 
#'   rowData = data.frame(col1 = 1:5, col2 = 6:10),
#'   colData = data.frame(cD1 = c(NA, NA), cD2 = c(NA, NA)))
#' t1 <- data.frame(col1 = 11:16, col2 = 17:22)
#' ## Define matches between row 1 in rowData(q1) with target row 2 and, 
#' ## rowData(q1) row 2 with target rows 2,3,4 and rowData(q1) row 5 with target 
#' ## row 5.
#' mo <- MatchedSummarizedExperiment(
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
#' res <- matchedData(mo, columns = c("col1", "col2", "target_col1", "target_col2"))
#' res
#' res$col1
#' res$target_col1
#' 
#' ## The example MatchedSummarizedExperiment object contains all rows in the 
#' ## rowData of the SummarizedExperiment and all target rows. Below we subset the 
#' ## object keeping only rows that are matched to at least one target row.
#' mo_sub <- mo[whichQuery(mo)]
#' 
#' ## mo_sub contains now a SummarizedExperiment with only 3 rows:
#' nrow(query(mo_sub))
#' 
#' ## while the original object contains a SummarizedExperiment with all 5 rows:
#' nrow(query(mo))
#' 
#' ## Both objects contain however still the full target object:
#' nrow(target(mo))
#' nrow(target(mo_sub))
#' 
#' ## With the `pruneTarget` we can however reduce also the target rows to
#' ## only those that match at least one in the rowData of query
#' mo_sub <- pruneTarget(mo_sub)
#' nrow(target(mo_sub))
NULL

setClass(
  "MatchedSummarizedExperiment",
  contains = "Matched"
)

#' @export
#'
#' @rdname MatchedSummarizedExperiment
MatchedSummarizedExperiment <- function(query = SummarizedExperiment(), 
                                        target = data.frame(),
                                        matches = data.frame(
                                          query_idx = integer(),
                                          target_idx = integer(),
                                          score = numeric())) {
  new("MatchedSummarizedExperiment", query = query, target = target,
      matches = matches)
}

setValidity("MatchedSummarizedExperiment", function(object) {
  msg <- NULL
  if(!is(object@query, "SummarizedExperiment"))
    msg <- c(msg, "query must be a SummarizedExperiment")
  if (length(msg)) return(msg)
  TRUE
})

#' @importMethodsFrom SummarizedExperiment rowData
#' 
#' @importFrom BiocGenerics colnames
#'
#' @exportMethod colnames
#'
#' @rdname MatchedSummarizedExperiment

setMethod("colnames", "MatchedSummarizedExperiment", function(x) {
  .colnames(rowData(x@query), x@target, x@matches)
})

#' @importMethodsFrom SummarizedExperiment rowData
#' 
#' @importMethodsFrom S4Vectors $
#'
#' @rdname MatchedSummarizedExperiment
#'
#' @export

setMethod("$", "MatchedSummarizedExperiment", function(x, name) {
  .dollar2(rowData(x@query), x@target, x@matches, name)
})

#' @importMethodsFrom S4Vectors cbind
#' 
#' @importMethodsFrom SummarizedExperiment rowData
#'
#' @importFrom S4Vectors DataFrame
#'
#' @rdname MatchedSummarizedExperiment
#'
#' @export

setMethod("matchedData", "MatchedSummarizedExperiment",
          function(object, columns = colnames(object), ...) {
  .matchedData2(rowData(object@query), object@target, object@matches, columns, 
               ...)
})
