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
           adducts = "adductClass"),
         contains = "MzParam",
         prototype = prototype(
           adducts = c("[M+H]+")),
         validity = function(object) {
           msg <- NULL
           if (is(object@adducts, "data.frame")) {
             if(any(!c("mass_add", "mass_multi") %in% colnames(object@adducts)))
               msg <- paste0("Columns \"mass_add\" and \"mass_multi\" must be ",
                             "present when adducts is a data.frame")
           } else {
             if (!all(object@adducts %in% c(adductNames("positive"),
                                            adductNames("negative"))))
               msg <- paste0("Unknown adducts, please check MetaboCoreUtils",
                             " for valid adducts")
           }
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
#'   name can be specified with parameter `mzColname` which defaults to `"mz"`)
#'   or a `numeric` with the m/z values of the features.
#'   If `target` is a `data.frame` it must
#'   contain a column with the (monoisotopic) mass of the compounds (the column
#'   name can be specified with parameter `massColname` which defaults to
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
#'   `mzColname` and `rtColname` which default to `"mz"` and `"rt"`,
#'   respectively). `target` must be a `data.frame` with the (monoisotopic)
#'   mass and retention time for each compound. The names of the columns
#'   containing these information can be specified with parameters `massColname`
#'   and `rtColname` which default to `"exactmass"` and `"rt"`, respectively.
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
#'   m/z values (which name can be specified with parameter `mzColname`,
#'   default being `mzColname = "mz"`) or a `numeric` with the m/z values of
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
#'   parameters `mzColname` and `rtColname` which default to `"mz"` and `"rt"`,
#'   respectively.`target` must be a `data.frame` again with the information on
#'   m/z and retention time.`MzRtParam` parameters `tolerance` and `ppm` allow
#'   to define the maximal acceptable (constant or m/z relative) difference
#'   between query and target m/z values; `MzRtParam` parameter `toleranceRt`
#'   allows to specify the maximal acceptable difference between query and
#'   target retention time values.
#'
#' @param adducts for `Mass2MzParam` or `Mass2MzRtParam`:
#'     either `character` with the names of adducts or `data.frame` with the
#'     adduct definition. This parameter is used to calculate m/z from target
#'     compounds' masses. Custom adduct definitions can be passed to the adduct
#'     parameter in form of a `data.frame`. This `data.frame` is expected to
#'     have columns `"mass_add"` and `"mass_multi"` defining the *additive* and
#'     *multiplicative* part of the calculation. See
#'     [MetaboCoreUtils::adducts()] for the expected format or use
#'     `MetaboCoreUtils::adductNames("positive")` and
#'     `MetaboCoreUtils::adductNames("negative")` for valid adduct names.
#'
#' @param BPPARAM parallel processing setup. See `BiocParallel::bpparam()` for
#'     details.
#'
#' @param massColname `character(1)` with the name of the column containing the
#'     mass of compounds. Defaults to `massColname = "exactmass"`.
#'
#' @param mzColname `character(1)` with the name of the column containing the
#'     m/z values. Defaults to `mzColname = "mz"`.
#'
#' @param param parameter object defining the matching approach and containing
#'     the settings for that approach. See description above for details.
#'
#' @param ppm for any `param` object: `numeric(1)` defining the maximal
#'     acceptable m/z-dependent difference (in parts-per-million) in m/z values
#'     to consider them to be *matching*.
#'
#' @param rtColname `character(1)` with the name of the column containing the
#'     retention times of the compounds. Defaults to `rtColname = "rt"`.
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
#' @return [Matched] object representing the result. To evaluate each match the
#' object contains the m/z error in ppm (variable `"ppm_error"`) as well as the
#' difference between the target and query m/z (variable `"score"`). The
#' difference between the target and query retention time (variable `"score_rt"`
#' is also present if retention time is considered for the match.
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
            Matched(query = query, target = target, matches = matches,
                    metadata = list(param = param))
          })
#' @rdname matchMz
setMethod("matchMz",
          signature = c(query = "numeric",
                        target = "data.frameOrSimilar",
                        param = "Mass2MzParam"),
          function(query, target, param, massColname = "exactmass",
                   BPPARAM = SerialParam()) {
            if (!massColname %in% colnames(target))
              stop("Missing column \"", massColname, "\" in target")
            res <- matchMz(query, target[, massColname], param)
            res@target <- target
            res
          })
#' @rdname matchMz
setMethod("matchMz",
          signature = c(query = "data.frameOrSimilar",
                        target = "numeric",
                        param = "Mass2MzParam"),
          function(query, target, param, BPPARAM = SerialParam(),
                   mzColname = "mz") {
            if (!mzColname %in% colnames(query))
              stop("Missing column \"", mzColname, "\" in query")
            res <- matchMz(query$mz, target, param)
            res@query <- query
            res
          })
#' @rdname matchMz
setMethod("matchMz",
          signature = c(query = "data.frameOrSimilar",
                        target = "data.frameOrSimilar",
                        param = "Mass2MzParam"),
          function(query, target, param, mzColname = "mz",
                   massColname = "exactmass", BPPARAM = SerialParam()) {
            if (!mzColname %in% colnames(query))
              stop("Missing column \"", mzColname, "\" in query")
            if (!massColname %in% colnames(target))
              stop("Missing column \"", massColname, "\" in target")
            res <- matchMz(query[, mzColname], target[, massColname], param)
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
            queryl <- length(query)
            res <- vector("list", queryl)
            for (i in seq_len(queryl)) {
                res[[i]] <- .getMatches(i, query[i], target = target_mz,
                                        tolerance = param@tolerance,
                                        ppm = param@ppm)
            }
            Matched(query = query, target = target,
                    matches = do.call(rbind, res),
                    metadata = list(param = param))
          })
#' @rdname matchMz
setMethod("matchMz",
          signature = c(query = "numeric",
                        target = "data.frameOrSimilar",
                        param = "MzParam"),
          function(query, target, param, mzColname = "mz",
                   BPPARAM = SerialParam()) {
            if (!mzColname %in% colnames(target))
              stop("Missing column \"", mzColname, "\" in target")
            res <- matchMz(query, target[, mzColname], param)
            res@target <- target
            res
          })
#' @rdname matchMz
setMethod("matchMz",
          signature = c(query = "data.frameOrSimilar",
                        target = "numeric",
                        param = "MzParam"),
          function(query, target, param, mzColname = "mz",
                   BPPARAM = SerialParam()) {
            if (!mzColname %in% colnames(query))
              stop("Missing column \"", mzColname, "\" in query")
            res <- matchMz(query[, mzColname], target, param)
            res@query <- query
            res
          })
#' @rdname matchMz
setMethod("matchMz",
          signature = c(query = "data.frameOrSimilar",
                        target = "data.frameOrSimilar",
                        param = "MzParam"),
          function(query, target, param, mzColname = c("mz", "mz"),
                   BPPARAM = SerialParam()) {
            if(length(mzColname) == 1)
              mzColname <- rep(mzColname, 2)
            if (!mzColname[1] %in% colnames(query))
              stop("Missing column \"", mzColname[1], "\" in query")
            if (!mzColname[2] %in% colnames(target))
              stop("Missing column \"", mzColname[2], "\" in target")
            res <- matchMz(query[, mzColname[1]], target[, mzColname[2]], param)
            res@query <- query
            res@target <- target
            res
          })
#' @rdname matchMz
#'
#' @importFrom BiocParallel bpmapply SerialParam
setMethod("matchMz",
          signature = c(query = "data.frameOrSimilar",
                        target = "data.frameOrSimilar",
                        param = "Mass2MzRtParam"),
          function(query, target, param, massColname = "exactmass",
                   mzColname = "mz", rtColname = c("rt", "rt"),
                   BPPARAM = SerialParam()) {
            if(length(rtColname) == 1)
              rtColname <- rep(rtColname, 2)
            if (!mzColname %in% colnames(query))
              stop("Missing column \"", mzColname, "\" in query")
            if (!rtColname[1] %in% colnames(query))
              stop("Missing column \"", rtColname[1], "\" in query")
            if (!massColname %in% colnames(target))
              stop("Missing column \"", massColname, "\" in target")
            if (!rtColname[2] %in% colnames(target))
              stop("Missing column \"", rtColname[2], "\" in target")
            target_mz <- .mass_to_mz_df(target[, massColname], param@adducts)
            target_mz$rt <- rep(target[, rtColname[2]], .nelements(param@adducts))
            queryl <- nrow(query)
            matches <- vector("list", queryl)
            query_mz <- query[, mzColname]
            query_rt <- query[, rtColname[1L]]
            for (i in seq_len(queryl)) {
                matches[[i]] <- .getMatchesMzRt(i, query_mz[i], query_rt[i],
                                                target = target_mz,
                                                tolerance = param@tolerance,
                                                ppm = param@ppm,
                                                toleranceRt = param@toleranceRt)
            }
            Matched(query = query, target = target,
                    matches = do.call(rbind, matches),
                    metadata = list(param = param))
          })
#' @rdname matchMz
#'
#' @importFrom BiocParallel bpmapply SerialParam
setMethod("matchMz",
          signature = c(query = "data.frameOrSimilar",
                        target = "data.frameOrSimilar",
                        param = "MzRtParam"),
          function(query, target, param, mzColname = c("mz", "mz"),
                   rtColname = c("rt", "rt"), BPPARAM = SerialParam()) {
            if(length(mzColname) == 1)
              mzColname <- rep(mzColname, 2)
            if(length(rtColname) == 1)
              rtColname <- rep(rtColname, 2)
            if (!mzColname[1] %in% colnames(query))
              stop("Missing column \"", mzColname[1], "\" in query")
            if (!mzColname[2] %in% colnames(target))
              stop("Missing column \"", mzColname[2], "\" in target")
            if (!rtColname[1] %in% colnames(query))
              stop("Missing column \"", rtColname[1], "\" in query")
            if (!rtColname[2] %in% colnames(target))
              stop("Missing column \"", rtColname[2], "\" in target")
            target_mz <- data.frame(index = seq_len(nrow(target)),
                                    mz = target[, mzColname[2]],
                                    rt = target[, rtColname[2]])
            queryl <- nrow(query)
            matches <- vector("list", queryl)
            query_mz <- query[, mzColname[1L]]
            query_rt <- query[, rtColname[1L]]
            for (i in seq_len(queryl)) {
                matches[[i]] <- .getMatchesMzRt(i, query_mz[i],
                                                query_rt[i],
                                                target = target_mz,
                                                tolerance = param@tolerance,
                                                ppm = param@ppm,
                                                toleranceRt = param@toleranceRt)
            }
            Matched(query = query, target = target,
                    matches = do.call(rbind, matches),
                    metadata = list(param = param))
          })

#' @rdname matchMz
setMethod("matchMz",
          signature = c(query = "SummarizedExperiment",
                        target = "ANY",
                        param = "Param"),
          function(query, target, param, mzColname = "mz", rtColname = "rt",
                   BPPARAM = SerialParam()) {
            matches <- matchMz(data.frame(rowData(query)), target, param,
                               mzColname , rtColname, BPPARAM)@matches
            MatchedSummarizedExperiment(query, target, matches,
                                        metadata = list(param = param))
          })

#' @author Andrea Vicini
#'
#' @param queryIndex `integer(1)` with the index of the query.
#'
#' @param queryMz `numeric(1)` with the m/z of the query.
#'
#' @param target `data.frame` with columns `"index"`, `"mz"` and optionally
#'     `"adduct"`.
#'
#' @noRd
.getMatches <- function(queryIndex, queryMz, target, tolerance, ppm){
  diffs <- queryMz - target$mz
  absdiffs <- abs(diffs)
  cls <- which(absdiffs <= (tolerance + ppm(queryMz, ppm)))
  if ("adduct" %in% colnames(target)){
    if (length(cls))
      data.frame(query_idx = queryIndex,
                 target_idx = target$index[cls],
                 adduct = target$adduct[cls],
                 score = diffs[cls],
                 ppm_error = absdiffs[cls] / target[cls, "mz"] * 10^6)
    else data.frame(query_idx = integer(),
                    target_idx = integer(),
                    adduct = character(),
                    score = numeric(),
                    ppm = numeric())
  } else {
    if (length(cls))
      data.frame(query_idx = queryIndex,
                 target_idx = target$index[cls],
                 score = diffs[cls],
                 ppm_error = absdiffs[cls] / target[cls, "mz"] * 10^6)
    else data.frame(query_idx = integer(),
                    target_idx = integer(),
                    score = numeric(),
                    ppm = numeric())
  }
}

#' @noRd
.getMatchesMzRt <- function(queryIndex, queryMz, queryRt, target, tolerance,
                            ppm, toleranceRt){
  diffs_rt <- queryRt - target$rt
  cls_rt <- which(abs(diffs_rt) <= toleranceRt)
  diffs <- queryMz - target$mz[cls_rt]
  absdiffs <- abs(diffs)
  cls <- which(absdiffs <= (tolerance + ppm(queryMz, ppm)))
  if ("adduct" %in% colnames(target)){
    if (length(cls))
      data.frame(query_idx = queryIndex,
                 target_idx = target$index[cls_rt[cls]],
                 adduct = target$adduct[cls_rt[cls]],
                 score = diffs[cls],
                 ppm_error = absdiffs[cls] / target[cls_rt[cls], "mz"] * 10^6,
                 score_rt = diffs_rt[cls_rt[cls]])
    else data.frame(query_idx = integer(),
                    target_idx = integer(),
                    adduct = character(),
                    score = numeric(),
                    ppm = numeric(),
                    score_rt = numeric())
  } else {
    if (length(cls))
      data.frame(query_idx = queryIndex,
                 target_idx = target$index[cls_rt[cls]],
                 score = diffs[cls],
                 ppm_error = absdiffs[cls] / target[cls_rt[cls], "mz"] * 10^6,
                 score_rt = diffs_rt[cls_rt[cls]])
    else data.frame(query_idx = integer(),
                    target_idx = integer(),
                    score = numeric(),
                    ppm = numeric(),
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
  data.frame(index = rep(seq_along(x), .nelements(adducts)),
             adduct = rep(colnames(mz), each = length(x)),
             mz = as.numeric(mz))
}
