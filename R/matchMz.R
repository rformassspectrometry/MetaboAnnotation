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
#'
#' @param adducts for `TargetMass2MzParam`: `character` with the names of
#'     adducts to calculate m/z from target compounds' masses. Use
#'     `MetaboCoreUtils::adductNames("positive")` and
#'     `MetaboCoreUtils::adductNames("negative")` for valid names.
#'
#' @param BPPARAM parallel processing setup. See `BiocParallel::bpparam()` for
#'     details.
#'
#' @param query feature table containing information on MS1 features. Can be
#'     a `data.frame` (with mandatory column names `"mz"`) or a `numeric` with
#'     the m/z values.
#'
#' @param target compound table with metabolites to compare against.
#'
#' @param param parameter object defining the matching approach and containing
#'     the settings for that approach. See description above for details.
#'
#' @param ppm for `TargetMass2MzParam`: `numeric(1)` defining the maximal
#'     acceptable m/z-dependent difference (in parts-per-million) in m/z values
#'     to consider them to be *matching*. See [MsCoreUtils::closest()] for
#'     more information.
#'
#' @param tolerance for `TargetMass2MzParam`: `numeric(1)` defining the maximal
#'     acceptable absolute difference in m/z values to consider them *matching*.
#'     See [MsCoreUtils::closest()] for more information.
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
#' ## same mass.
#'
#' ## We repeat the matching requiring an exact match
#' parm <- TargetMass2MzParam(
#'     adducts = c("[M+H]+", "[M+Na]+"),
#'     tolerance = 0,
#'     ppm = 0)
#' res <- matchMz(fts, target_df, parm)
#' res
#'
#' data(res)
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
#' data(res)
#'
#' ## No matches were found.
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
                        target = "data.frame",
                        param = "TargetMass2MzParam"),
          function(query, target, param, BPPARAM = SerialParam()) {
              if (!"exactmass" %in% colnames(target))
                  stop("Missing column \"exactmass\" in target")
              target_mz <- .mass_to_mz_df(target$exactmass, param@adducts)
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
          signature = c(query = "data.frame",
                        target = "data.frame",
                        param = "TargetMass2MzParam"),
          function(query, target, param, BPPARAM = SerialParam()) {
              if (!"mz" %in% colnames(query))
                  stop("Missing column \"mz\" in query")
              res <- matchMz(query$mz, target, param)
              res@query <- query
              res
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
    if (length(cls))
        data.frame(query_idx = queryIndex,
                   target_idx = target$index[cls],
                   adduct = target$adduct[cls],
                   score = abs(queryMz - target$mz[cls]))
    else data.frame(query_idx = integer(),
                    target_idx = integer(),
                    adduct = character(),
                    score = numeric())
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
