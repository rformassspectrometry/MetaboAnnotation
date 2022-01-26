setClass(
  "MatchedSummarizedExperiment",
  contains = "Matched",
  slots = c(query = "SummarizedExperiment")
)

#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
#' @exportClass MatchedSummarizedExperiment
#'
#' @export
#'
#' @rdname Matched
MatchedSummarizedExperiment <- function(query = SummarizedExperiment(),
                                        target = data.frame(),
                                        matches = data.frame(
                                          query_idx = integer(),
                                          target_idx = integer(),
                                          score = numeric()),
                                        metadata = list()) {
  new("MatchedSummarizedExperiment", query = query, target = target,
      matches = matches, metadata = metadata)
}

setValidity("MatchedSummarizedExperiment", function(object) {
  msg <- NULL
  if (length(msg)) return(msg)
  TRUE
})

#' @importMethodsFrom SummarizedExperiment rowData
#'
#' @importFrom BiocGenerics colnames
#'
#' @exportMethod colnames
#'
#' @rdname hidden_aliases
setMethod("colnames", "MatchedSummarizedExperiment", function(x) {
  .colnames(rowData(x@query), x@target, x@matches)
})

#' @importMethodsFrom SummarizedExperiment rowData
#'
#' @importMethodsFrom S4Vectors $
#'
#' @rdname hidden_aliases
#'
#' @export
setMethod("$", "MatchedSummarizedExperiment", function(x, name) {
  .dollar(rowData(x@query), x@target, x@matches, name)
})

#' @importMethodsFrom S4Vectors cbind
#'
#' @importMethodsFrom SummarizedExperiment rowData
#'
#' @importFrom S4Vectors DataFrame
#'
#' @rdname hidden_aliases
#'
#' @export
setMethod("matchedData", "MatchedSummarizedExperiment",
          function(object, columns = colnames(object), ...) {
  .matchedData(rowData(object@query), object@target, object@matches, columns,
               ...)
})

#' @importFrom methods validObject
#'
#' @rdname hidden_aliases
#'
#' @export
setMethod("filterMatches", "MatchedSummarizedExperiment",
          function(object, queryValue = integer(), targetValue = integer(),
                   queryColname = character(), targetColname = character(),
                   index = integer(), keep = TRUE, ...) {
            if(length(index) && any(!index%in%seq_len(nrow(object@matches))))
              stop("some indexes in \"index\" are out of bounds")
            if(!length(index) && length(queryValue))
              index  <- .findMatchesIdxs(rowData(object@query), object@target,
                                        object@matches, queryValue,
                                        targetValue, queryColname,
                                        targetColname)
            if(keep) to_keep <- seq_len(nrow(object@matches)) %in% index
            else to_keep <- !seq_len(nrow(object@matches)) %in% index
            object@matches <- object@matches[to_keep, , drop = FALSE]
            validObject(object)
            object
          })

#' @importFrom methods validObject
#'
#' @rdname hidden_aliases
#'
#' @export
setMethod("addMatches", "MatchedSummarizedExperiment",
          function(object, queryValue = integer(), targetValue = integer(),
                   queryColname = character(), targetColname = character(),
                   score = data.frame(), isIndex = FALSE) {
            object@matches <- .addMatches(rowData(object@query), object@target,
                                          object@matches, queryValue,
                                          targetValue, queryColname,
                                          targetColname, score, isIndex)
            validObject(object)
            object
          })
