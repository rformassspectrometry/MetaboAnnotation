#' @title formula matching
#'
#' @name matchFormula
#'
#' @description
#'
#' The `matchFormula` method matches chemical formula from different inputs
#' (parameter `query` and `target`). Before comparison all formulas are
#' normalized using [MetaboCoreUtils::standardizeFormula()]. Inputs can be either
#' a `character` or `data.frame` containing a column with formulas. In case of
#' a `data.frame` as input the parameter `formulaColName` needs to specified.
#'
#' @param query chemical formulas to search
#'
#' @param target compound table to compare against
#'
#' @param formulaColname name of the column contain chemical formula. Only used
#'     if a `data.frame` is supplied.
#'
#' @param ... currently ignored
#'
#' @return [Matched] object representing the result
#'
#' @author Michael Witting
#'
#' @examples
#'
#' ## input formula
#' query <- c("H12C6O6", "C11H12O2", "HN3")
#' target <- c("HCl", "C2H4O", "C6H12O6")
#'
#' query_df <- data.frame(
#'     formula = c("H12C6O6", "C11H12O2", "HN3"),
#'     name = c("A", "B", "C")
#' )
#' target_df <- data.frame(
#'     formula = c("HCl", "C2H4O", "C6H12O6"),
#'     name = c("D", "E", "F")
#' )
#'
#' ## character vs character
#' matches <- matchFormula(query, target)
#' matchedData(matches)
#'
#' ## data.frame vs data.frame
#' matches <- matchFormula(query_df, target_df)
#' matchedData(matches)
#' ## data.frame vs character
#' matches <- matchFormula(query_df, target)
#' matchedData(matches)
#' ## character vs data.frame
#' matches <- matchFormula(query, target_df)
#' matchedData(matches)
NULL

#' @rdname matchFormula
#'
#' @export
setGeneric("matchFormula", function(query, target, ...)
  standardGeneric("matchFormula"))

#' @rdname matchFormula
#'
#' @importFrom MetaboCoreUtils standardizeFormula
setMethod("matchFormula",
          signature = c(query = "character",
                        target = "character"),
          function(query, target, BPPARAM = SerialParam()) {
            query <- unname(standardizeFormula(query))
            target <- unname(standardizeFormula(target))
            target_formula <- data.frame(index = rep(seq_along(target)),
                                         formula = target)
            matches <- do.call(
              rbind, bpmapply(seq_along(query), query, FUN = .getFormulaMatches,
                              MoreArgs = list(target = target_formula),
                              BPPARAM = BPPARAM, SIMPLIFY = FALSE))
            Matched(query = query, target = target, matches = matches)

          })
#' @rdname matchFormula
setMethod("matchFormula",
          signature = c(query = "data.frameOrSimilar",
                        target = "data.frameOrSimilar"),
          function(query, target, formulaColname = c("formula", "formula"),
                   BPPARAM = SerialParam()) {
            if(length(formulaColname) == 1)
              formulaColname <- rep(formulaColname, 2)
            if(!formulaColname[1] %in% colnames(query))
              stop("Missing column \"", formulaColname[1], "\" in query")
            if(!formulaColname[2] %in% colnames(target))
              stop("Missing column \"", formulaColname[2], "\" in target")
            res <- matchFormula(query[, formulaColname[1]],
                                target[, formulaColname[2]])
            res@query <- query
            res@target <- target
            res
          })

#' @rdname matchFormula
setMethod("matchFormula",
          signature = c(query = "character",
                        target = "data.frameOrSimilar"),
          function(query, target, formulaColname = "formula",
                   BBPARAM = SerialParam()) {
            if(!formulaColname %in% colnames(target))
              stop("Missing column \"", formulaColname, "\" in target")
            res <- matchFormula(query, target[, formulaColname])
            res@query <- query
            res@target <- target
            res
          })

#' @rdname matchFormula
setMethod("matchFormula",
          signature = c(query = "data.frameOrSimilar",
                        target = "character"),
          function(query, target, formulaColname = "formula",
                   BBPARAM = SerialParam()) {
            if(!formulaColname %in% colnames(query))
              stop("Missing column \"", formulaColname, "\" in query")
            res <- matchFormula(query[, formulaColname], target)
            res@query <- query
            res@target <- target
            res
          })


#' @author Michael Witting
#'
#' @param queryIndex `integer(1)` with the index of the query
#'
#' @param queryFormula `character(1)` with the formula of the query
#'
#' @param target `data.frame` with columns `"index"` and  `"formula"`
#'
#' @noRd
.getFormulaMatches <- function(queryIndex, queryFormula, target) {

  cls <- which(queryFormula == target$formula)

  if(length(cls)) {
    data.frame(query_idx = queryIndex,
               target_idx = cls,
               score = 1)
  } else {
    data.frame(query_idx = integer(),
              target_idx = integer(),
              score = numeric())
  }
}

