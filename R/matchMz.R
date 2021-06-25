#' @noRd
setClass("MzParam",
         slots = c(
           tolerance = "numeric",
           ppm = "numeric"),
         contains = "Param",
         prototype = prototype(
           tolerance = 0,
           ppm = 5),
         validity = function(object) {
           msg <- NULL
           if (length(object@tolerance) != 1 || object@tolerance < 0)
             msg <- c("'tolerance' has to be a positive number of length 1")
           if (length(object@ppm) != 1 || object@ppm < 0)
             msg <- c("'ppm' has to be a positive number of length 1")
           msg
         })

#' @rdname matchMz
#'
#' @importFrom methods new
#'
#' @export
MzParam <- function(tolerance = 0, ppm = 5) {
  new("MzParam", tolerance = tolerance, ppm = ppm)
}

#' @importFrom MetaboCoreUtils adductNames
#'
#' @noRd
setClass("Mass2MzParam",
         slots = c(
           adducts = "character"),
         contains = "MzParam",
         prototype = prototype(
           adducts = c("[M+H]+")),
         validity = function(object) {
           msg <- NULL
           if (!all(object@adducts %in% c(adductNames("positive"),
                                          adductNames("negative"))))
             msg <- paste0("Unknown adducts, please check MetaboCoreUtils",
                           " for valid adducts")
           msg
         })

#' @rdname matchMz
#'
#' @importFrom methods new
#'
#' @export
Mass2MzParam <- function(adducts = c("[M+H]+"), tolerance = 0, ppm = 5) {
  new("Mass2MzParam", adducts = adducts, tolerance = tolerance,
      ppm = ppm)
}

#' @importFrom MetaboCoreUtils adductNames
#'
#' @noRd
setClass("Mass2MzRtParam",
         slots = c(
           toleranceRt = "numeric"),
         contains = "Mass2MzParam",
         prototype = prototype(
           toleranceRt = 0),
         validity = function(object) {
           msg <- NULL
           if (length(object@toleranceRt) != 1 || object@toleranceRt < 0)
             msg <- c("'toleranceRt' has to be a positive number of length 1")
           msg
         })

#' @rdname matchMz
#'
#' @importFrom methods new
#'
#' @export
Mass2MzRtParam <- function(adducts = c("[M+H]+"), tolerance = 0, ppm = 5,
                                 toleranceRt = 0) {
  new("Mass2MzRtParam", adducts = adducts, tolerance = tolerance,
      ppm = ppm, toleranceRt = toleranceRt)
}

#' @noRd
setClass("MzRtParam",
         slots = c(
           toleranceRt = "numeric"),
         contains = "MzParam",
         prototype = prototype(
           toleranceRt = 0),
         validity = function(object) {
           msg <- NULL
           if (length(object@toleranceRt) != 1 || object@toleranceRt < 0)
             msg <- c("'toleranceRt' has to be a positive number of length 1")
           msg
         })

#' @rdname matchMz
#'
#' @importFrom methods new
#'
#' @export
MzRtParam <- function(tolerance = 0, ppm = 0, toleranceRt = 0) {
  new("MzRtParam", tolerance = tolerance, ppm = ppm, toleranceRt = toleranceRt)
}

#' @title m/z matching
#'
#' @name matchMz
#'
#' @description
#'
#' The `matchMz` method matches (compares) m/z values from a MS1 data table
#' (parameter `query`) with theoretical m/z values for (reference) compounds
#' (parameter `target`) considering or not also retention times. The approach
#' which is used to perform the comparison and additional settings for this
#' matching can be defined with `param`.
#'
#' Available matching approaches and respective `param` objects are:
#'
#' - `Mass2MzParam`: match m/z values against reference compounds for
#'   which the (exact) mass is known. Before matching, m/z values are calculated
#'   from the mass of the compounds in the *target* table for the specified
#'   adducts (parameter `adduct` in the parameter object).
#'   `query` must be a `data.frame` with a column containing m/z values (column
#'   name can be specified with parameter `mzColumn` which defaults to `"mz"`)
#'   or a `numeric` with the m/z values of the features.
#'   If `target` is a `data.frame` it must
#'   contain a column with the (monoisotopic) mass of the compounds (the column
#'   name can be specified with parameter `massColumn` which defaults to
#'   `"exactmass"`) `Mass2MzParam`'s parameter `adducts` allows to define
#'   the expected adducts (defaults to `adducts = "[M+H]+"` but any adducts
#'   available in [MetaboCoreUtils::adducts()] are supported). Parameter
#'   `tolerance` and `ppm` allow to define the maximal acceptable (constant or
#'   m/z relative) difference between query and target m/z values.
#'
#' - `Mass2MzRtParam`: match m/z and retention time values against
#'   reference compounds for which the (exact) mass **and** retention time are
#'   known. Before matching, m/z values are calculated from the mass of the
#'   compounds in the *target* table for the specified adducts (parameter
#'   `adduct` in the parameter object). The retention time is considered equal
#'   for different adducts of the same compound. `query` must be a `data.frame`
#'   with m/z and retention time values of the features. The names of the
#'   columns containing these information can be specified with parameters
#'   `mzColumn` and `rtColumn` which default to `"mz"` and `"rt"`,
#'   respectively). `target` must be a `data.frame` with the (monoisotopic)
#'   mass and retention time for each compound. The names of the columns
#'   containing these information can be specified with parameters `massColumn`
#'   and `rtColumn` which default to `"exactmass"` and `"rt"`, respectively.
#'   `Mass2MzRtParam`'s parameter `adducts` allows to define the expected
#'   adducts (defaults to `adducts = "[M+H]+"` but any adducts available in
#'   [MetaboCoreUtils::adducts()] are supported). Parameter `tolerance` and
#'   `ppm` allow to define the maximal acceptable (constant or m/z relative)
#'   difference between query and target m/z values; parameter `toleranceRt`
#'   allows to specify the maximal acceptable difference between query and
#'   target retention time values.
#'
#' - `MzParam`: match m/z values against reference compounds for which the m/z
#'   is known. `query` must be either a `data.frame` with a column containing
#'   m/z values (which name can be specified with parameter `mzColumn`,
#'   default being `mzColumn = "mz"`) or a `numeric` with the m/z values of
#'   the features. The same holds for `target`. `MzParam` parameters
#'   `tolerance` and `ppm` allow to define the maximal acceptable (constant or
#'   m/z relative) difference between query and target m/z values.
#'
#' - `MzRtParam`: match m/z and retention time values against reference
#'   compounds for which m/z and retention time are known. `query` must be a
#'   `data.frame` or a `SummarizedExperiment`. The `data.frame` in one case or 
#'   the `SummarizedExperiment` `rowData` in the other must have columns 
#'   containing the m/z and retention times of the
#'   features. The names of the respective columns can be specified with
#'   parameters `mzColumn` and `rtColumn` which default to `"mz"` and `"rt"`,
#'   respectively.`target` must be a `data.frame` again with the information on 
#'   m/z and retention time.`MzRtParam` parameters `tolerance` and `ppm` allow 
#'   to define the maximal acceptable (constant or m/z relative) difference 
#'   between query and target m/z values; `MzRtParam` parameter `toleranceRt` 
#'   allows to specify the maximal acceptable difference between query and 
#'   target retention time values.
#'
#' @param adducts for `Mass2MzParam` or `Mass2MzRtParam`:
#'     `character` with the names of adducts to calculate m/z from target
#'     compounds' masses. Use `MetaboCoreUtils::adductNames("positive")` and
#'     `MetaboCoreUtils::adductNames("negative")` for valid names.
#'
#' @param BPPARAM parallel processing setup. See `BiocParallel::bpparam()` for
#'     details.
#'
#' @param massColumn `character(1)` with the name of the column containing the
#'     mass of compounds. Defaults to `massColumn = "exactmass"`.
#'
#' @param mzColumn `character(1)` with the name of the column containing the
#'     m/z values. Defaults to `mzColumn = "mz"`.
#'
#' @param param parameter object defining the matching approach and containing
#'     the settings for that approach. See description above for details.
#'
#' @param ppm for any `param` object: `numeric(1)` defining the maximal
#'     acceptable m/z-dependent difference (in parts-per-million) in m/z values
#'     to consider them to be *matching*.
#'
#' @param rtColumn `character(1)` with the name of the column containing the
#'     retention times of the compounds. Defaults to `rtColumn = "rt"`.
#'
#' @param target compound table with metabolites to compare against.
#'
#' @param tolerance for any `param` object: `numeric(1)` defining the maximal
#'     acceptable absolute difference in m/z values to consider them *matching*.
#'
#' @param toleranceRt for `Mass2MzRtParam` or `MzRtParam`: `numeric(1)`
#'     defining the maximal acceptable absolute difference in retention time
#'     values to consider them them *matching*.
#'
#' @param query feature table containing information on MS1 features. Can be
#'     a `data.frame` (with mandatory column names `"mz"`) or a `numeric` with
#'     the m/z values. A matching based on both m/z and retention time can be
#'     performed when a column `"rt"` is present in both `query` and `target`.
#'
#' @param ... currently ignored.
#'
#' @return [Matched] object representing the result
#'
#' @author Andrea Vicini, Michael Witting
#'
#' @seealso [matchSpectra] or [CompareSpectraParam()] for spectra data matching
#'
#' @examples
#'
#' library(MetaboCoreUtils)
#' ## Create a simple "target/reference" compound table
#' target_df <- data.frame(
#'    name = c("Tryptophan", "Leucine", "Isoleucine"),
#'    formula = c("C11H12N2O2", "C6H13NO2", "C6H13NO2"),
#'    exactmass = c(204.089878, 131.094629, 131.094629)
#' )
#'
#' ## Create a "feature" table with m/z of features. We calculate m/z for
#' ## certain adducts of some of the compounds in the reference table.
#' fts <- data.frame(
#'     feature_id = c("FT001", "FT002", "FT003"),
#'     mz = c(mass2mz(204.089878, "[M+H]+"),
#'            mass2mz(131.094629, "[M+H]+"),
#'            mass2mz(204.089878, "[M+Na]+") + 1e-6))
#'
#' ## Define the parameters for the matching
#' parm <- Mass2MzParam(
#'     adducts = c("[M+H]+", "[M+Na]+"),
#'     tolerance = 0,
#'     ppm = 20)
#' res <- matchMz(fts, target_df, parm)
#' res
#'
#' ## Get the full matching result:
#' matchedData(res)
#'
#' ## We have thus matches of FT002 to two different compounds (but with the
#' ## same mass).
#'
#' ## We repeat the matching requiring an exact match
#' parm <- Mass2MzParam(
#'     adducts = c("[M+H]+", "[M+Na]+"),
#'     tolerance = 0,
#'     ppm = 0)
#' res <- matchMz(fts, target_df, parm)
#' res
#'
#' matchedData(res)
#'
#' ## The last feature could thus not be matched to any compound.
#'
#' ## At last we use also different adduct definitions.
#' parm <- Mass2MzParam(
#'     adducts = c("[M+K]+", "[M+Li]+"),
#'     tolerance = 0,
#'     ppm = 20)
#' res <- matchMz(fts, target_df, parm)
#' res
#'
#' matchedData(res)
#'
#' ## No matches were found.
#'
#' ## We can also match a "feature" table with a target data.frame taking into
#' ## account both m/z and retention time values.
#' target_df <- data.frame(
#'   name = c("Tryptophan", "Leucine", "Isoleucine"),
#'   formula = c("C11H12N2O2", "C6H13NO2", "C6H13NO2"),
#'   exactmass = c(204.089878, 131.094629, 131.094629),
#'   rt = c(150, 140, 140)
#' )
#'
#' fts <- data.frame(
#'   feature_id = c("FT001", "FT002", "FT003"),
#'   mz = c(mass2mz(204.089878, "[M+H]+"),
#'          mass2mz(131.094629, "[M+H]+"),
#'          mass2mz(204.089878, "[M+Na]+") + 1e-6),
#'   rt = c(150, 140, 150.1)
#' )
#'
#' ## Define the parameters for the matching
#' parm <- Mass2MzRtParam(
#'   adducts = c("[M+H]+", "[M+Na]+"),
#'   tolerance = 0,
#'   ppm = 20,
#'   toleranceRt = 0)
#'
#' res <- matchMz(fts, target_df, parm)
#' res
#'
#' ## Get the full matching result:
#' matchedData(res)
#'
#' ## FT003 could not be matched to any compound, FT002 was matched to two
#' ## different compounds (but with the same mass).
#'
#' ## We repeat the matching allowing a positive tolerance for the matches
#' ## between rt values
#'
#' ## Define the parameters for the matching
#' parm <- Mass2MzRtParam(
#'   adducts = c("[M+H]+", "[M+Na]+"),
#'   tolerance = 0,
#'   ppm = 20,
#'   toleranceRt = 0.1)
#'
#' res <- matchMz(fts, target_df, parm)
#' res
#'
#' ## Get the full matching result:
#' matchedData(res)
#'
#' ## Also FT003 was matched in this case
#'
#' ## It is also possible to match directly m/z values
#' mz1 <- c(12, 343, 23, 231)
#' mz2 <- mz1 + rnorm(4, sd = 0.001)
#'
#' res <- matchMz(mz1, mz2, MzParam(tolerance = 0.001))
#'
#' matchedData(res)
#' 
#' ## Should I add here an example where `query` is a `SummarizedExperiment`? 
NULL

#' @rdname matchMz
#'
#' @export
setGeneric("matchMz", function(query, target, param, ...)
  standardGeneric("matchMz"))
#' @rdname matchMz
#'
#' @importFrom BiocParallel bpmapply SerialParam
setMethod("matchMz",
          signature = c(query = "numeric",
                        target = "numeric",
                        param = "Mass2MzParam"),
          function(query, target, param, BPPARAM = SerialParam()) {
            target_mz <- .mass_to_mz_df(target, param@adducts)
            matches <- do.call(
              rbind, bpmapply(seq_along(query), query, FUN = .getMatches,
                              MoreArgs = list(target = target_mz,
                                              tolerance = param@tolerance,
                                              ppm = param@ppm),
                              BPPARAM = BPPARAM, SIMPLIFY = FALSE))
            Matched(query = query, target = target, matches = matches)
          })
#' @rdname matchMz
setMethod("matchMz",
          signature = c(query = "numeric",
                        target = "data.frame",
                        param = "Mass2MzParam"),
          function(query, target, param, massColumn = "exactmass",
                   BPPARAM = SerialParam()) {
            if (!massColumn %in% colnames(target))
              stop("Missing column \"", massColumn, "\" in target")
            res <- matchMz(query, target[, massColumn], param)
            res@target <- target
            res
          })
#' @rdname matchMz
setMethod("matchMz",
          signature = c(query = "data.frame",
                        target = "numeric",
                        param = "Mass2MzParam"),
          function(query, target, param, BPPARAM = SerialParam()) {
            if (!"mz" %in% colnames(query))
              stop("Missing column \"mz\" in query")
            res <- matchMz(query$mz, target, param)
            res@query <- query
            res
          })
#' @rdname matchMz
setMethod("matchMz",
          signature = c(query = "data.frame",
                        target = "data.frame",
                        param = "Mass2MzParam"),
          function(query, target, param, mzColumn = "mz",
                   massColumn = "exactmass", BPPARAM = SerialParam()) {
            if (!mzColumn %in% colnames(query))
              stop("Missing column \"", mzColumn, "\" in query")
            if (!massColumn %in% colnames(target))
              stop("Missing column \"", massColumn, "\" in target")
            res <- matchMz(query[, mzColumn], target[, massColumn], param)
            res@query <- query
            res@target <- target
            res
          })
#' @rdname matchMz
#'
#' @importFrom BiocParallel bpmapply SerialParam
setMethod("matchMz",
          signature = c(query = "numeric",
                        target = "numeric",
                        param = "MzParam"),
          function(query, target, param, BPPARAM = SerialParam()) {
            target_mz <- data.frame(index = seq_along(target),
                                    mz = target)
            matches <- do.call(
              rbind, bpmapply(seq_along(query), query, FUN = .getMatches,
                              MoreArgs = list(target = target_mz,
                                              tolerance = param@tolerance,
                                              ppm = param@ppm),
                              BPPARAM = BPPARAM, SIMPLIFY = FALSE))
            Matched(query = query, target = target, matches = matches)
          })
#' @rdname matchMz
setMethod("matchMz",
          signature = c(query = "numeric",
                        target = "data.frame",
                        param = "MzParam"),
          function(query, target, param, mzColumn = "mz",
                   BPPARAM = SerialParam()) {
            if (!mzColumn %in% colnames(target))
              stop("Missing column \"", mzColumn, "\" in target")
            res <- matchMz(query, target[, mzColumn], param)
            res@target <- target
            res
          })
#' @rdname matchMz
setMethod("matchMz",
          signature = c(query = "data.frame",
                        target = "numeric",
                        param = "MzParam"),
          function(query, target, param, mzColumn = "mz",
                   BPPARAM = SerialParam()) {
            if (!mzColumn %in% colnames(query))
              stop("Missing column \"", mzColumn, "\" in query")
            res <- matchMz(query[, mzColumn], target, param)
            res@query <- query
            res
          })
#' @rdname matchMz
setMethod("matchMz",
          signature = c(query = "data.frame",
                        target = "data.frame",
                        param = "MzParam"),
          function(query, target, param, mzColumn = "mz",
                   BPPARAM = SerialParam()) {
            if (!mzColumn %in% colnames(query))
              stop("Missing column \"", mzColumn, "\" in query")
            if (!mzColumn %in% colnames(target))
              stop("Missing column \"", mzColumn, "\" in target")
            res <- matchMz(query[, mzColumn], target[, mzColumn], param)
            res@query <- query
            res@target <- target
            res
          })
#' @rdname matchMz
#'
#' @importFrom BiocParallel bpmapply SerialParam
setMethod("matchMz",
          signature = c(query = "data.frame",
                        target = "data.frame",
                        param = "Mass2MzRtParam"),
          function(query, target, param, massColumn = "exactmass",
                   mzColumn = "mz", rtColumn = "rt", BPPARAM = SerialParam()) {
            if (!mzColumn %in% colnames(query))
              stop("Missing column \"", mzColumn, "\" in query")
            if (!rtColumn %in% colnames(query))
              stop("Missing column \"", rtColumn, "\" in query")
            if (!massColumn %in% colnames(target))
              stop("Missing column \"", massColumn, "\" in target")
            if (!rtColumn %in% colnames(target))
              stop("Missing column \"", rtColumn, "\" in target")
            target_mz <- .mass_to_mz_df(target[, massColumn], param@adducts)
            target_mz$rt <- rep(target[, rtColumn], length(param@adducts))
            matches <- do.call(
                rbind, bpmapply(seq_len(nrow(query)), query[, mzColumn],
                                query[, rtColumn],
                                FUN = .getMatchesMzRt,
                                MoreArgs = list(target = target_mz,
                                                tolerance = param@tolerance,
                                                ppm = param@ppm,
                                                toleranceRt = param@toleranceRt),
                                BPPARAM = BPPARAM, SIMPLIFY = FALSE))
            Matched(query = query, target = target, matches = matches)
          })
#' @rdname matchMz
#'
#' @importFrom BiocParallel bpmapply SerialParam
setMethod("matchMz",
          signature = c(query = "data.frame",
                        target = "data.frame",
                        param = "MzRtParam"),
          function(query, target, param, mzColumn = "mz", rtColumn = "rt",
                   BPPARAM = SerialParam()) {
            if (!mzColumn %in% colnames(query))
              stop("Missing column \"", mzColumn, "\" in query")
            if (!rtColumn %in% colnames(query))
              stop("Missing column \"", rtColumn, "\" in query")
            target_mz <- data.frame(index = seq_len(nrow(target)),
                                    mz = target[, mzColumn],
                                    rt = target[, rtColumn])
            matches <- do.call(
                rbind, bpmapply(seq_len(nrow(query)), query[, mzColumn],
                                query[, rtColumn],
                                FUN = .getMatchesMzRt,
                                MoreArgs = list(target = target_mz,
                                                tolerance = param@tolerance,
                                                ppm = param@ppm,
                                                toleranceRt = param@toleranceRt),
                                BPPARAM = BPPARAM, SIMPLIFY = FALSE))
            Matched(query = query, target = target, matches = matches)
          })
#' @rdname matchMz
#'
#' @importFrom BiocParallel bpmapply SerialParam
setMethod("matchMz",
          signature = c(query = "SummarizedExperiment",
                        target = "data.frame",
                        param = "MzRtParam"),
          function(query, target, param, mzColumn = "mz", rtColumn = "rt",
                   BPPARAM = SerialParam()) {
            rD_query <- rowData(query)
            if (!mzColumn %in% colnames(rD_query))
              stop("Missing column \"", mzColumn, "\" in rowData(query)")
            if (!rtColumn %in% colnames(rD_query))
              stop("Missing column \"", rtColumn, "\" in rowData(query)")
            target_mz <- data.frame(index = seq_len(nrow(target)),
                                    mz = target[, mzColumn],
                                    rt = target[, rtColumn])
            matches <- do.call(
              rbind, bpmapply(seq_len(nrow(rD_query)), rD_query[, mzColumn],
                              rD_query[, rtColumn],
                              FUN = .getMatchesMzRt,
                              MoreArgs = list(target = target_mz,
                                              tolerance = param@tolerance,
                                              ppm = param@ppm,
                                              toleranceRt = param@toleranceRt),
                              BPPARAM = BPPARAM, SIMPLIFY = FALSE))
            MatchedSummarizedExperiment(query = query, target = target, 
                                        matches = matches)
          })

#' @author Andrea Vicini
#'
#' @noRd
.getMatches <- function(queryIndex, queryMz, target, tolerance, ppm){
  diffs <- abs(queryMz - target$mz)
  cls <- which(diffs <= (tolerance + ppm(queryMz, ppm)))
  if ("adduct" %in% colnames(target)){
    if (length(cls))
      data.frame(query_idx = queryIndex,
                 target_idx = target$index[cls],
                 adduct = target$adduct[cls],
                 score = diffs[cls])
    else data.frame(query_idx = integer(),
                    target_idx = integer(),
                    adduct = character(),
                    score = numeric())
  } else {
    if (length(cls))
      data.frame(query_idx = queryIndex,
                 target_idx = target$index[cls],
                 score = diffs[cls])
    else data.frame(query_idx = integer(),
                    target_idx = integer(),
                    score = numeric())
  }
}

#' @noRd
.getMatchesMzRt <- function(queryIndex, queryMz, queryRt, target, tolerance,
                            ppm, toleranceRt){
  diffs_rt <- abs(queryRt - target$rt)
  cls_rt <- which(abs(queryRt - target$rt) <= toleranceRt)
  diffs <- abs(queryMz - target$mz[cls_rt])
  cls <- which(abs(queryMz - target$mz[cls_rt]) <= (tolerance + ppm(queryMz, ppm)))
  if ("adduct" %in% colnames(target)){
    if (length(cls))
      data.frame(query_idx = queryIndex,
                 target_idx = target$index[cls_rt[cls]],
                 adduct = target$adduct[cls_rt[cls]],
                 score = diffs[cls],
                 score_rt = diffs_rt[cls_rt[cls]])
    else data.frame(query_idx = integer(),
                    target_idx = integer(),
                    adduct = character(),
                    score = numeric(),
                    score_rt = numeric())
  } else {
    if (length(cls))
      data.frame(query_idx = queryIndex,
                 target_idx = target$index[cls_rt[cls]],
                 score = diffs[cls],
                 score_rt = diffs_rt[cls_rt[cls]])
    else data.frame(query_idx = integer(),
                    target_idx = integer(),
                    score = numeric(),
                    score_rt = numeric())
  }
}

#' creates a `data.frame` with columns `"index"`, `"adduct"` and `"mz"` from
#' provided mass values `x` and adduct definition `adducts`. `"index"` is the
#' index of the mass in the input vector `x`.
#'
#' @importFrom MetaboCoreUtils mass2mz
#'
#' @noRd
.mass_to_mz_df <- function(x, adducts) {
  mz <- mass2mz(x, adducts)
  data.frame(index = rep(seq_along(x), length(adducts)),
             adduct = rep(colnames(mz), each = length(x)),
             mz = as.numeric(mz))
}
