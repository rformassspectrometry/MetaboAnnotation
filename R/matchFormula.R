#' @title Chemical Formula Matching
#'
#' @name matchFormula
#'
#' @description
#'
#' The `matchFormula` method matches chemical formulas from different inputs
#' (parameter `query` and `target`). Before comparison all formulas are
#' normalized using [MetaboCoreUtils::standardizeFormula()]. Inputs can be
#' either a `character` or `data.frame` containing a column with formulas.
#' In case of `data.frame`s parameter `formulaColname` needs to be used to
#' specify the name of the column containing the chemical formulas.
#'
#' @param BPPARAM parallel processing setup. See `BiocParallel::bpparam()` for
#'     details.
#'
#' @param formulaColname `character` with the name of the column containing
#'     chemical formulas. Can be of length 1 if both `query` and `target` are
#'     `data.frame`s and the name of the column with chemical formulas is the
#'     same for both. If different columns are used, `formulaColname[1]` can be
#'     used to define the column name in `query` and `formulaColname[2]` the
#'     one of `target`.
#'
#' @param target `character` or `data.frame` with chemical formulas to compare
#'     against.
#'
#' @param query `character` or `data.frame` with chemical formulas to search.
#'
#' @param ... currently ignored
#'
#' @return [Matched] object representing the result.
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
#' 
#' @importFrom BiocParallel bpmapply SerialParam
setMethod(
    "matchFormula",
    signature = c(query = "character",
                  target = "character"),
    function(query, target, BPPARAM = SerialParam()) {
        matches <- do.call(
            rbind, bpmapply(seq_along(query), standardizeFormula(query),
                            FUN = .getFormulaMatches,
                            MoreArgs =list(target = standardizeFormula(target)),
                            BPPARAM = BPPARAM, SIMPLIFY = FALSE))
        Matched(query = query, target = target, matches = matches)
    })

#' @rdname matchFormula
setMethod("matchFormula",
          signature = c(query = "data.frameOrSimilar",
                        target = "data.frameOrSimilar"),
          function(query, target, formulaColname = c("formula", "formula"),
                   BPPARAM = SerialParam()) {
            if(length(formulaColname) == 1L)
              formulaColname <- rep(formulaColname, 2L)
            if(!formulaColname[1] %in% colnames(query))
              stop("Missing column \"", formulaColname[1L], "\" in query")
            if(!formulaColname[2L] %in% colnames(target))
              stop("Missing column \"", formulaColname[2L], "\" in target")
            res <- matchFormula(query[, formulaColname[1L]],
                                target[, formulaColname[2L]],
                                BPPARAM = BPPARAM)
            res@query <- query
            res@target <- target
            res
          })

#' @rdname matchFormula
setMethod("matchFormula",
          signature = c(query = "character",
                        target = "data.frameOrSimilar"),
          function(query, target, formulaColname = "formula",
                   BPPARAM = SerialParam()) {
            if(!formulaColname %in% colnames(target))
              stop("Missing column \"", formulaColname, "\" in target")
            res <- matchFormula(query, target[, formulaColname],
                                BPPARAM = BPPARAM)
            res@query <- query
            res@target <- target
            res
          })

#' @rdname matchFormula
setMethod("matchFormula",
          signature = c(query = "data.frameOrSimilar",
                        target = "character"),
          function(query, target, formulaColname = "formula",
                   BPPARAM = SerialParam()) {
            if(!formulaColname %in% colnames(query))
              stop("Missing column \"", formulaColname, "\" in query")
            res <- matchFormula(query[, formulaColname], target,
                                BPPARAM = BPPARAM)
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
#' @param target `character` with target formulas
#'
#' @noRd
.getFormulaMatches <- function(queryIndex, queryFormula, target) {
    cls <- which(target == queryFormula)
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
