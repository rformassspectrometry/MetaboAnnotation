setClass(
  "MatchedSummarizedExperiment",
  contains = "Matched"
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
