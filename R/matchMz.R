#' @importFrom MetaboCoreUtils adductNames
#'
#' @noRd
setClass("TargetMass2MzParam",
         slots = c(
             adducts = "character",
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
             if(!all(object@adducts %in% c(adductNames("positive"),
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
TargetMass2MzParam <- function(adducts = c("[M+H]+"), tolerance = 0, ppm = 5) {
  new("TargetMass2MzParam", adducts = adducts, tolerance = tolerance,
      ppm = ppm)
}


#' @importFrom MetaboCoreUtils adductNames
#'
#' @noRd
setClass("TargetMass2MzRtParam",
         slots = c(
           adducts = "character",
           tolerance = "numeric",
           ppm = "numeric",
           toleranceRt = "numeric"),
         contains = "Param",
         prototype = prototype(
           tolerance = 0,
           ppm = 5,
           toleranceRt = 0),
         validity = function(object) {
           msg <- NULL
           if (length(object@tolerance) != 1 || object@tolerance < 0)
             msg <- c("'tolerance' has to be a positive number of length 1")
           if (length(object@ppm) != 1 || object@ppm < 0)
             msg <- c("'ppm' has to be a positive number of length 1")
           if (length(object@toleranceRt) != 1 || object@toleranceRt < 0)
             msg <- c("'toleranceRt' has to be a positive number of length 1")
           if(!all(object@adducts %in% c(adductNames("positive"),
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
TargetMass2MzRtParam <- function(adducts = c("[M+H]+"), tolerance = 0, ppm = 5,
                                 toleranceRt = 0) {
  new("TargetMass2MzRtParam", adducts = adducts, tolerance = tolerance,
      ppm = ppm, toleranceRt = toleranceRt)
}


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

#' @noRd
setClass("MzRtParam",
         slots = c(
           tolerance = "numeric",
           ppm = "numeric",
           toleranceRt = "numeric"),
         contains = "Param",
         prototype = prototype(
           tolerance = 0,
           ppm = 5,
           toleranceRt = 0),
         validity = function(object) {
           msg <- NULL
           if (length(object@tolerance) != 1 || object@tolerance < 0)
             msg <- c("'tolerance' has to be a positive number of length 1")
           if (length(object@ppm) != 1 || object@ppm < 0)
             msg <- c("'ppm' has to be a positive number of length 1")
           if (length(object@toleranceRt) != 1 || object@toleranceRt < 0)
             msg <- c("'toleranceRt' has to be a positive number of length 1")
           msg
         })

#' @rdname matchMz
#'
#' @importFrom methods new
#'
#' @export
MzRtParam <- function(tolerance = 0, ppm = 5, toleranceRt = 0) {
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
#' (parameter `target`). The approach which is used to perform the comparison
#' and additional settings for this matching can be defined with `param`.
#'
#' Available matching approaches and respective `param` objects are:
#'
#' - `TargetMass2MzParam`: match m/z values against reference compounds for
#'   which the (exact) mass is known. Before matching, m/z values are calculated
#'   from the mass of the compounds in the *target* table for the specified
#'   adducts (parameter `adduct` in the parameter object).
#'   `query` must be a `data.frame` with a column named `"mz"` or a `numeric`
#'   with the m/z values of the features. If `target` is a `data.frame` it must
#'   contain a column `"exactmass"` with the exact monoisotopic mass for each
#'   compound. `TargetMass2MzParam`'s parameter `adducts` allows to define the
#'   expected adducts (defaults to `adducts = "[M+H]+" but any adducts
#'   available in [MetaboCoreUtils::adducts()] are supported). Parameter
#'   `tolerance` and `ppm` allow to define the maximal acceptable (constant or
#'   m/z relative) difference between query and target m/z values.
#' - `TargetMass2MzRtParam`: match m/z and retention time values against 
#'   reference compounds for which the (exact) mass and retention time are 
#'   known. Before matching, m/z values are calculated from the mass of the 
#'   compounds in the *target* table for the specified adducts (parameter 
#'   `adduct` in the parameter object). The retention time is considered equal 
#'   for different adducts of the same compound. `query` must be a `data.frame` 
#'   with columns `"mz"`with the m/z values of the features and `"rt"` with 
#'   their retention time. `target` must be a `data.frame` with columns 
#'   `"exactmass"` with the exact monoisotopic mass for each compound and `"rt"` 
#'   with the associated retention time. `TargetMass2MzRtParam`'s parameter 
#'   `adducts` allows to define the expected adducts (defaults to 
#'   `adducts = "[M+H]+" but any adducts available in 
#'   [MetaboCoreUtils::adducts()] are supported). Parameter `tolerance` and 
#'   `ppm` allow to define the maximal acceptable (constant or m/z relative) 
#'   difference between query and target m/z values; parameter `toleranceRt` 
#'   allows to specify the maximal acceptable difference between query and 
#'   target retention time values.
#' - `MzParam`: match m/z values against reference compounds for which the m/z 
#'   is known. `query` must be either a `data.frame` with a column named `"mz"` 
#'   or a `numeric` with the m/z values of the features. The same holds for 
#'   `target`. `MzParam` parameters `tolerance` and `ppm` allow to define the 
#'   maximal acceptable (constant or m/z relative) difference between query and 
#'   target m/z values.
#' - `MzRtParam`: match m/z and retention time values against reference 
#'   compounds for which m/z and retention time are known.`query` must be a 
#'   `data.frame` with columns `"mz"` with the m/z values of the features and 
#'   `"rt"` with their retention time. The same holds for `target`.`MzRtParam` 
#'   parameters `tolerance` and `ppm` allow to define the maximal acceptable 
#'   (constant or m/z relative) difference between query and target m/z values; 
#'   `MzRtParam` parameter `toleranceRt` allows to specify the maximal 
#'   acceptable difference between query and target retention time values.  
#'
#' @param adducts for `TargetMass2MzParam` or `TargetMass2MzRtParam`: `character` 
#'     with the names of adducts to calculate m/z from target compounds' masses. 
#'     Use `MetaboCoreUtils::adductNames("positive")` and 
#'     `MetaboCoreUtils::adductNames("negative")` for valid names.
#'
#' @param BPPARAM parallel processing setup. See `BiocParallel::bpparam()` for
#'     details.
#'
#' @param query feature table containing information on MS1 features. Can be
#'     a `data.frame` (with mandatory column names `"mz"`) or a `numeric` with
#'     the m/z values. A matching based on both m/z and retention time can be 
#'     performed when a column `"rt"` is present in both `query` and `target`.
#'
#' @param target compound table with metabolites to compare against.
#'
#' @param param parameter object defining the matching approach and containing
#'     the settings for that approach. See description above for details.
#'
#' @param ppm for any `param` object: `numeric(1)` defining the maximal
#'     acceptable m/z-dependent difference (in parts-per-million) in m/z values
#'     to consider them to be *matching*. See [MsCoreUtils::closest()] for
#'     more information.
#'
#' @param tolerance for any `param` object: `numeric(1)` defining the maximal
#'     acceptable absolute difference in m/z values to consider them *matching*.
#'     See [MsCoreUtils::closest()] for more information.
#'     
#' @param toleranceRt for `TargetMass2MzRtParam` or `MzRtParam`: `numeric(1)` 
#' defining the maximal acceptable absolute difference in retention time values 
#' to consider them them *matching*.
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
#' parm <- TargetMass2MzParam(
#'     adducts = c("[M+H]+", "[M+Na]+"),
#'     tolerance = 0,
#'     ppm = 20)
#' res <- matchMz(fts, target_df, parm)
#' res
#'
#' ## Get the full matching result:
#' data(res)
#'
#' ## We have thus matches of FT002 to two different compounds (but with the
#' ## same mass).
#'
#' ## We repeat the matching requiring an exact match
#' parm <- TargetMass2MzParam(
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
#' parm <- TargetMass2MzParam(
#'     adducts = c("[M+K]+", "[M+Li]+"),
#'     tolerance = 0,
#'     ppm = 20)
#' res <- matchMz(fts, target_df, parm)
#' res
#'
#' matchedData(res)
#'
#' ## No matches were found.
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
#' parm <- TargetMass2MzRtParam(
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
#' parm <- TargetMass2MzRtParam(
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
                        param = "TargetMass2MzParam"),
          function(query, target, param, BPPARAM = SerialParam()) {
              target_mz <- .mass_to_mz_df(target, param@adducts)
              target_mz <- target_mz[order(target_mz$mz), ]
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
                        param = "TargetMass2MzParam"),
          function(query, target, param, BPPARAM = SerialParam()) {
            if (!"exactmass" %in% colnames(target))
              stop("Missing column \"exactmass\" in target")
            res <- matchMz(query, target$exactmass, param)
            res@target <- target
          })
#' @rdname matchMz
setMethod("matchMz",
          signature = c(query = "data.frame",
                        target = "numeric",
                        param = "TargetMass2MzParam"),
          function(query, target, param, BPPARAM = SerialParam()) {
            if (!"mz" %in% colnames(query))
              stop("Missing column \"mz\" in query")
            res <- matchMz(query$mz, target, param)
            res@query <- query
          })
#' @rdname matchMz
setMethod("matchMz",
          signature = c(query = "data.frame",
                        target = "data.frame",
                        param = "TargetMass2MzParam"),
          function(query, target, param, BPPARAM = SerialParam()) {
              if (!"mz" %in% colnames(query))
                  stop("Missing column \"mz\" in query")
            if (!"exactmass" %in% colnames(target))
              stop("Missing column \"exactmass\" in target")
              res <- matchMz(query$mz, target$exactmass, param)
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
            target_mz <- target_mz[order(target_mz$mz), ]
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
          function(query, target, param, BPPARAM = SerialParam()) {
            if (!"mz" %in% colnames(target))
              stop("Missing column \"mz\" in target")
            res <- matchMz(query, target$mz, param)
            res@target <- target
            res
          })
#' @rdname matchMz
setMethod("matchMz",
          signature = c(query = "data.frame",
                        target = "numeric",
                        param = "MzParam"),
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
                        param = "MzParam"),
          function(query, target, param, BPPARAM = SerialParam()) {
            if (!"mz" %in% colnames(query))
              stop("Missing column \"mz\" in query")
            if (!"mz" %in% colnames(target))
              stop("Missing column \"mz\" in target")
            res <- matchMz(query$mz, target$mz, param)
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
                        param = "TargetMass2MzRtParam"),
          function(query, target, param, BPPARAM = SerialParam()) {
            if (!"mz" %in% colnames(query))
              stop("Missing column \"mz\" in query")
            if (!"rt" %in% colnames(query))
              stop("Missing column \"rt\" in query")
            if (!"exactmass" %in% colnames(target))
              stop("Missing column \"exactmass\" in target")
            if (!"rt" %in% colnames(target))
              stop("Missing column \"rt\" in target")
            target_mz <- .mass_to_mz_df(target$exactmass, param@adducts)
            target_mz$rt <- rep(target$rt, length(param@adducts))
            target_mz <- target_mz[order(target_mz$mz), ]
            matches <- do.call(
              rbind, bpmapply(seq_len(nrow(query)), query$mz, query$rt, 
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
          function(query, target, param, BPPARAM = SerialParam()) {
            if (!"mz" %in% colnames(query))
              stop("Missing column \"mz\" in query")
            if (!"rt" %in% colnames(query))
              stop("Missing column \"rt\" in query")
            target_mz <- data.frame(index = seq_len(nrow(target)), 
                                    mz = target$mz, rt = target$rt)
            target_mz <- target_mz[order(target_mz$mz), ]
            matches <- do.call(
              rbind, bpmapply(seq_len(nrow(query)), query$mz, query$rt, 
                              FUN = .getMatchesMzRt,
                              MoreArgs = list(target = target_mz,
                                              tolerance = param@tolerance,
                                              ppm = param@ppm,
                                              toleranceRt = param@toleranceRt),
                              BPPARAM = BPPARAM, SIMPLIFY = FALSE))
            Matched(query = query, target = target, matches = matches)
          })
#' @author Andrea Vicini
#'
#' @importFrom MsCoreUtils closest
#'
#' @noRd
.getMatches <- function(queryIndex, queryMz, target, tolerance, ppm){
  cls <- closest(target$mz, queryMz, tolerance = tolerance, ppm = ppm,
                 .check = FALSE)
  cls <- which(!is.na(cls))
  if ("adduct" %in% colnames(target)){
    if (length(cls))
      data.frame(query_idx = queryIndex,
                 target_idx = target$index[cls],
                 adduct = target$adduct[cls],
                 score = abs(queryMz - target$mz[cls]))
    else data.frame(query_idx = integer(),
                    target_idx = integer(),
                    adduct = character(),
                    score = numeric())
  } else {
    if (length(cls))
      data.frame(query_idx = queryIndex,
                 target_idx = target$index[cls],
                 score = abs(queryMz - target$mz[cls]))
    else data.frame(query_idx = integer(),
                    target_idx = integer(),
                    score = numeric())
  }
}

#' @noRd
.getMatchesMzRt <- function(queryIndex, queryMz, queryRt, target, tolerance, 
                            ppm, toleranceRt){
  cls_rt <- which(abs(queryRt - target$rt) <= toleranceRt)
  cls <- closest(target$mz[cls_rt], queryMz, tolerance = tolerance, ppm = ppm,
                 .check = FALSE)
  cls <- which(!is.na(cls))
  if ("adduct" %in% colnames(target)){
    if (length(cls))
      data.frame(query_idx = queryIndex,
                 target_idx = target$index[cls_rt[cls]],
                 adduct = target$adduct[cls_rt[cls]],
                 score = abs(queryMz - target$mz[cls_rt[cls]]),
                 score_rt = abs(queryRt - target$rt[cls_rt[cls]]))
    else data.frame(query_idx = integer(),
                    target_idx = integer(),
                    adduct = character(),
                    score = numeric(),
                    score_rt = numeric())
  } else {
    if (length(cls))
      data.frame(query_idx = queryIndex,
                 target_idx = target$index[cls_rt[cls]],
                 score = abs(queryMz - target$mz[cls_rt[cls]]),
                 score_rt = abs(queryRt - target$rt[cls_rt[cls]]))
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
