#' @noRd
setClass("ValueParam",
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
ValueParam <- function(tolerance = 0, ppm = 5) {
    new("ValueParam", tolerance = tolerance, ppm = ppm)
}

#' @noRd
setClass("MzParam", contains = "ValueParam")

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
             targetAdducts = "adductClass"),
         contains = "ValueParam",
         prototype = prototype(
             targetAdducts = c("[M+H]+")),
         validity = function(object) {
             .valid_adduct(object@targetAdducts, "`targetAdducts`")
         })

#' @rdname matchMz
#'
#' @importFrom methods new
#'
#' @export
Mass2MzParam <- function(adducts = c("[M+H]+"), tolerance = 0, ppm = 5) {
  new("Mass2MzParam", targetAdducts = adducts, tolerance = tolerance,
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
  new("Mass2MzRtParam", targetAdducts = adducts, tolerance = tolerance,
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

#' @importFrom MetaboCoreUtils adductNames
#'
#' @noRd
setClass("Mz2MassParam",
         slots = c(
           queryAdducts = "adductClass",
           targetAdducts = "adductClass"),
         contains = "ValueParam",
         prototype = prototype(
           queryAdducts = c("[M+H]+"),
           targetAdducts = c("[M-H]-")),
         validity = function(object) {
           c(.valid_adduct(object@queryAdducts, "`queryAdducts`"),
             .valid_adduct(object@targetAdducts, "`targetAdducts`"))
         })

#' @rdname matchMz
#'
#' @importFrom methods new
#'
#' @export
Mz2MassParam <- function(queryAdducts = c("[M+H]+"),
                             targetAdducts = c("[M-H]-"),
                             tolerance = 0, ppm = 5) {
  new("Mz2MassParam", queryAdducts = queryAdducts,
      targetAdducts = targetAdducts, tolerance = tolerance, ppm = ppm)
}


#' @noRd
setClass("Mz2MassRtParam",
         slots = c(
             toleranceRt = "numeric"),
         contains = "Mz2MassParam",
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
Mz2MassRtParam <- function(queryAdducts = c("[M+H]+"),
                           targetAdducts = c("[M+H]+"),
                           tolerance = 0, ppm = 5, toleranceRt = 0) {
    new("Mz2MassRtParam", queryAdducts = queryAdducts,
        targetAdducts = targetAdducts, tolerance = tolerance, ppm = ppm,
        toleranceRt = toleranceRt)
}

#' @title m/z matching
#'
#' @name matchMz
#'
#' @description
#' 
#' The `matchMz` method matches elements from `query` with those in `target`
#' using different matching approaches depending on parameter `param`.
#' Generally, `query` is expected to contain MS experimental values
#' (m/z and possibly retention time) while `target` theoretical values from
#' reference compounds. `target` can be `numeric` or `data.frameOrSimilar`
#' (i.e. `data.frame`, `matrix` or `DataFrame`). The same is true for `query`
#' which in addition can also be `SummarizedExperiment`.
#'
#' Available `param` objects and corresponding matching approaches are:
#'
#' - `Mass2MzParam`: match m/z values against reference compounds for
#'   which the (exact) mass is known. Before matching, m/z values are calculated
#'   from the compounds masses in the *target* table using the adducts
#'   specified via `Mass2MzParam` `adducts` parameter (defaults to
#'   `adducts = "[M+H]+"`). `query` must be a `data.frameOrSimilar` with a
#'   column containing m/z values (column name can be specified with parameter
#'   `mzColname` which defaults to `"mz"`) or a `numeric` with the m/z values of
#'   the features. `target` must be a `data.frameOrSimilar` with a column
#'   containing (monoisotopic) masses of the reference compounds (the column
#'   name can be specified with parameter `massColname` which defaults
#'   to `"exactmass"`) or a `numeric` with the mass values.
#'   Parameters `tolerance` and `ppm` allow to define the maximal acceptable
#'   (constant or m/z relative) difference between m/zs in `query` and those
#'   obtained from masses in `target`.
#'
#' - `Mass2MzRtParam`: match m/z and retention time values against
#'   reference compounds for which the (exact) mass **and** retention time are
#'   known. Before matching, m/z values are calculated from the compounds masses
#'   in the *target* table using the adducts specified via
#'   `Mass2MzParam` `adducts` parameter (defaults to `adducts = "[M+H]+"`).
#'   The retention time is considered equal for different adducts of the 
#'   same compound. `query` must be a `data.frameOrSimilar` with m/z and
#'   retention time values of the features. The names of the columns containing
#'   these information can be specified with parameters `mzColname` and
#'   `rtColname` which default to `"mz"` and `"rt"`, respectively.
#'   `target` must be a `data.frameOrSimilar` with the (monoisotopic)
#'   masses and retention times for each reference compound. The names of the
#'   columns containing this information can be specified with parameters
#'   `massColname` and `rtColname` which default to `"exactmass"` and `"rt"`,
#'   respectively. `Mass2MzRtParam` parameters `tolerance` and `ppm` have
#'   the same meaning as in `Mass2MzParam`; parameter `toleranceRt` allows to
#'   specify the maximal acceptable difference between query and target
#'   retention time values.
#'
#' - `MzParam`: match m/z values against reference compounds for which the m/z
#'   is known. `query` must be either a `data.frameOrSimilar` with a column
#'   containing m/z values (whose name can be specified with parameter
#'   `mzColname`, default being `mzColname = "mz"`) or a `numeric` with the m/z
#'   values of the features. The same holds for `target`. `MzParam` parameters
#'   `tolerance` and `ppm` allow to define the maximal acceptable (constant or
#'   m/z relative) difference between query and target m/z values.
#'
#' - `MzRtParam`: match m/z and retention time values against reference
#'   compounds for which m/z and retention time are known. `query` must be a
#'   `data.frameOrSimilar` or a `SummarizedExperiment` `query` in the first
#'   case or its `rowData` in the second must have columns containing the m/zs
#'   and retention times of the features. `target` must be a
#'   `data.frameOrSimilar` again with columns containing m/zs and retention
#'   times. The names of such columns can be specified with parameters
#'   `mzColname` and `rtColname` (defaulting to `c("mz", "mz")` and
#'   `c("rt", "rt")`, respectively). `MzRtParam` parameters `tolerance` and
#'   `ppm` have the same meaning as in `MzParam`; `MzRtParam` parameter
#'   `toleranceRt` allows to specify the maximal acceptable difference between
#'   query and target retention time values.
#'
#' - `Mz2MassParam`: first convert to masses the m/z values of `query` and
#'   `target` based respectively on `Mz2MassParam` parameters `queryAdducts`
#'   (defaults to `queryAdducts = "[M+H]+"`) and `targetAdducts` (defaults to
#'   `targetAdducts = "[M-H]-"`) and then match the obtained masses. 
#'   `query` must be either a `data.frameOrSimilar` with a column containing
#'   m/z values (whose name can be specified with parameter `mzColname`,
#'   default being `mzColname = "mz"`) or a `numeric` with the m/z values.
#'   The same holds for `target`. Parameters `tolerance` and `ppm` allow to
#'   define the maximal acceptable (constant or mass relative) difference
#'   between masses derived from `query` and `target`.
#'
#' - `Mz2MassRtParam`: first convert to masses the m/z values of `query` and
#'   `target` based respectively on `Mz2MassRtParam` parameters `queryAdducts`
#'   (defaults to `queryAdducts = "[M+H]+"`) and `targetAdducts` (defaults to
#'   `targetAdducts = "[M-H]-"`) and then match the obtained masses and
#'   corresponding retention times (to each mass value corresponds the same
#'   retention time associated to the m/z from which the mass value was
#'   obtained). `query` must be a `data.frameOrSimilar` with a column containing
#'   m/z values (whose name can be specified with parameter `mzColname`,
#'   default being `mzColname = "mz"`). The same holds for `target`.
#'   Parameters `tolerance` and `ppm` have the same meaning as in
#'   `Mz2MassParam`. Parameter `toleranceRt` allows to specify the maximal
#'   acceptable difference between retention time values.
#'
#' - `ValueParam`: matches elements from from `query` (if `numeric`) or else
#'   from one of its columns with those from `target` (if `numeric`) or else
#'   from one of its columns. The name(s) of the colum(s) used for the matching
#'   (if `query` and `target` are not bot numeric has/have to be specified via
#'   parameter `valueColname`.
#'
#' @param adducts for `Mass2MzParam` or `Mass2MzRtParam`: either `character`
#'     with adduct names from [MetaboCoreUtils::adducts()] or `data.frame` with
#'     a custom adduct definition. This parameter is used to calculate m/z from
#'     target compounds' masses. Custom adduct definitions can be passed to the
#'     adduct parameter in form of a `data.frame`. This `data.frame` is expected
#'     to have columns `"mass_add"` and `"mass_multi"` defining the *additive*
#'     and *multiplicative* part of the calculation. See
#'     [MetaboCoreUtils::adducts()] for the expected format or use
#'     `MetaboCoreUtils::adductNames("positive")` and
#'     `MetaboCoreUtils::adductNames("negative")` for valid adduct names.
#'
#' @param massColname `character(1)` with the name of the column in `target`
#'     containing the mass of compounds. To be used when `param` is
#'     `Mass2MzParam` or `Mass2MzRtParam` (and target is not already `numeric`
#'     with the masses). Defaults to `massColname = "exactmass"`.
#'
#' @param mzColname `character` specifying the name(s) of the column(s) in
#'     `query` or/and `target`with the m/z values. If one among `query` and
#'     `target` is `numeric` (and therefore there is no need to specify the
#'     column name) or `query` is not `numeric` and `param` is `Mass2MzParam`
#'     or `Mass2MzRtParam` (and therefore the name of the column with m/z needs
#'     only to be specified for `query`) then `mzColname` is expected to be
#'     `character(1)`. If both `query` and `target` are not numeric `mzColname`
#'     is expected to be `character(2)` (or `character(1)` and in this last case
#'     the two column names are assumed to be the same). If not specified the
#'     assumed default name for coulmns with m/z values is `"mz"`. 
#'
#' @param param parameter object defining the matching approach and containing
#'     the settings for that approach. See description above for details.
#'
#' @param ppm for any `param` object: `numeric(1)` defining the maximal
#'     acceptable m/z-dependent (or mass-dependent for `Mz2MassParam`)
#'     difference (in parts-per-million) in m/z values to consider them to
#'     be *matching*.
#'     
#' @param query feature table containing information on MS1 features. Can be
#'     a `data.frame`, `DataFrame`, `matrix`, `SummarizedExperiment` or
#'     `numeric`. It is expected to contain m/z values and can contain also
#'     other variables. Matchings based on both m/z and retention time can be
#'     performed when a column with retention times is present in both `query`
#'     and `target`.
#'     
#' @param queryAdducts for `Mz2MassParam`. Adducts used to derive mass
#'     values from query m/z values. The expected format is the same as that
#'     for parameter `adducts`.
#'
#' @param rtColname `character(2)` with the name of the column containing
#'     the compounds retention times in `query` and the name for the one in
#'     `target`. It can also be `character(1)` if the two names are the same.
#'     To be used when `param` is `MzRtParam` or `Mass2MzRtParam`.  
#'     Defaults to `rtColname = c("rt", "rt")`.
#'
#' @param target compound table with metabolites to compare against.
#' 
#' @param targetAdducts for `Mz2MassParam`. Adducts used to derive mass
#'     values from target m/z values. The expected format is the same as that
#'     for parameter `adducts`.
#'
#' @param tolerance for any `param` object: `numeric(1)` defining the maximal
#'     acceptable absolute difference in m/z (or in mass for `Mz2MassParam`) 
#'     to consider them *matching*.
#'
#' @param toleranceRt for `Mass2MzRtParam` or `MzRtParam`: `numeric(1)`
#'     defining the maximal acceptable absolute difference in retention time
#'     values to consider them them *matching*.
#'
#' @param valueColname `character` specifying the name of the column in
#'     `query` or/and the one in `target`with the desired values for the
#'     matching. This parameter should only be used when `param` is
#'     `valueParam` and in this case it must be provided (unless both `query`
#'     and `target` are `numeric`). It can be `character(1)` or `character(2)`
#'     in a similar way to `mzColname`.
#'
#' @param ... currently ignored.
#'
#' @return [Matched] object representing the result (or more precisely
#' [MatchedSummarizedExperiment] if `query` is `SummarizedExperiment`).
#' To evaluate each match the object contains, depending on which `param` is
#' used, the m/z or mass error in ppm (variable `"ppm_error"`) as well as
#' the m/z or mass difference (variable `"score"`). When `param` is
#' `Mz2MassParam` such difference and ppm error are obtained using mass values
#' internally computed from m/z values in `query` and `target` (see above
#' for more detail on `Mz2MassParam`). When `parm` is either `MzParam`,
#' `MzRtParam`, `Mass2MzParam` or `Mass2MzRtParam` `"ppm_error"` and `"score"`
#' are computed from m/z values. In the first two cases both the m/z values are
#' contained in `query` and `target` parameters. In the last two m/z values are
#' internally computed from reference mass values in `target`. See more details
#' for each individual parameter above. Additionally, if `param` is either
#' `MzRtParam` or `Mass2MzRtParam` also retention time is used for the matching
#' and the difference in retention time for each of the matched elements is also
#' added to the returned `Matched` object. Note that, for a match, a negative
#' value of `"score"` (or `"score_rt"`) indicates that the m/z or mass (or
#' retention time) associated to the query element is smaller than that
#' associated to the target element.
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
setMethod("matchMz",
          signature = c(query = "numeric",
                        target = "numeric",
                        param = "ValueParam"),
          function(query, target, param) {
              target_ <- data.frame(index = seq_along(target), value = target)
              queryl <- length(query)
              res <- vector("list", queryl)
              for (i in seq_len(queryl)) {
                  res[[i]] <- .getMatches(i, query[i], target = target_,
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
                        param = "ValueParam"),
          function(query, target, param, valueColname = character()) {
              if(!length(valueColname))
                  stop("`valueColname` has to be provided.")
              if (!valueColname %in% colnames(target))
                  stop("Missing column \"", valueColname, "\" in target")
              res <- matchMz(query, target[, valueColname], param)
              res@target <- target
              res
          })

#' @rdname matchMz
setMethod("matchMz",
          signature = c(query = "data.frameOrSimilar",
                        target = "numeric",
                        param = "ValueParam"),
          function(query, target, param, valueColname = character()) {
              if(!length(valueColname))
                  stop("`valueColname` has to be provided.")
              if (!valueColname %in% colnames(query))
                  stop("Missing column \"", valueColname, "\" in query")
              res <- matchMz(query[, valueColname], target, param)
              res@query <- query
              res
          })

#' @rdname matchMz
setMethod("matchMz",
          signature = c(query = "data.frameOrSimilar",
                        target = "data.frameOrSimilar",
                        param = "ValueParam"),
          function(query, target, param, valueColname = character()) {
              if(!length(valueColname))
                  stop("`valueColname` has to be provided.")
              if(length(valueColname) == 1)
                  valueColname <- rep(valueColname, 2)
              if (!valueColname[1] %in% colnames(query))
                  stop("Missing column \"", valueColname[1], "\" in query")
              if (!valueColname[2] %in% colnames(target))
                  stop("Missing column \"", valueColname[2], "\" in target")
              res <- matchMz(query[, valueColname[1]],
                             target[, valueColname[2]], param)
              res@query <- query
              res@target <- target
              res
          })

#' @rdname matchMz
setMethod("matchMz",
          signature = c(query = "numeric",
                        target = "numeric",
                        param = "Mass2MzParam"),
          function(query, target, param) {
              target_mz <- .mass_to_mz_df(target, param@targetAdducts)
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
                        param = "Mass2MzParam"),
          function(query, target, param, massColname = "exactmass") {
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
          function(query, target, param, mzColname = "mz") {
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
                   massColname = "exactmass") {
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
setMethod("matchMz",
          signature = c(query = "numeric",
                        target = "data.frameOrSimilar",
                        param = "MzParam"),
          function(query, target, param, mzColname = "mz") {
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
          function(query, target, param, mzColname = "mz") {
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
          function(query, target, param, mzColname = c("mz", "mz")) {
              if(length(mzColname) == 1)
                  mzColname <- rep(mzColname, 2)
              if (!mzColname[1] %in% colnames(query))
                  stop("Missing column \"", mzColname[1], "\" in query")
              if (!mzColname[2] %in% colnames(target))
                  stop("Missing column \"", mzColname[2], "\" in target")
              res <- matchMz(query[, mzColname[1]],
                             target[, mzColname[2]], param)
              res@query <- query
              res@target <- target
              res
          })

#' @rdname matchMz
setMethod("matchMz",
          signature = c(query = "data.frameOrSimilar",
                        target = "data.frameOrSimilar",
                        param = "Mass2MzRtParam"),
          function(query, target, param, massColname = "exactmass",
                   mzColname = "mz", rtColname = c("rt", "rt")) {
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
              target_mz <- .mass_to_mz_df(target[, massColname],
                                          param@targetAdducts)
              target_mz$rt <- rep(target[, rtColname[2]],
                                  .nelements(param@targetAdducts))
              queryl <- nrow(query)
              matches <- vector("list", queryl)
              query_mz <- query[, mzColname]
              query_rt <- query[, rtColname[1L]]
              for (i in seq_len(queryl)) {
                  matches[[i]] <-
                      .getMatchesMzRt(i, query_mz[i], query_rt[i],
                                      target = target_mz,
                                      tolerance = param@tolerance,
                                      ppm = param@ppm,
                                      toleranceRt =param@toleranceRt)
              }
              Matched(query = query, target = target,
                      matches = do.call(rbind, matches),
                      metadata = list(param = param))
          })

#' @rdname matchMz
setMethod("matchMz",
          signature = c(query = "data.frameOrSimilar",
                        target = "data.frameOrSimilar",
                        param = "MzRtParam"),
          function(query, target, param, mzColname = c("mz", "mz"),
                   rtColname = c("rt", "rt")) {
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
                  matches[[i]] <-
                      .getMatchesMzRt(i, query_mz[i],
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
          signature = c(query = "numeric",
                        target = "numeric",
                        param = "Mz2MassParam"),
          function(query, target, param) {
              query_mass <- .mz_to_mass_df(query, param@queryAdducts)
              target_mass<- .mz_to_mass_df(target, param@targetAdducts)
              queryl <- nrow(query_mass)
              matches <- vector("list", queryl)
              for(i in seq_len(queryl)) {
                  matches[[i]] <- .getMatches(query_mass$index[i],
                                              query_mass$mass[i],
                                              target = target_mass,
                                              tolerance = param@tolerance,
                                              ppm = param@ppm)
                  if (nrow(matches[[i]]))
                      matches[[i]]$query_adduct <- query_mass$adduct[i]
                  else matches[[i]]$query_adduct <- character()
              }
              matches <- do.call(rbind, matches)
              colnames(matches)[3] <- "target_adduct"
              Matched(query = query, target = target,
                      matches = matches[, c(1, 2, 6, 3, 4, 5)],
                      metadata = list(param = param))
              
          })

#' @rdname matchMz
setMethod("matchMz",
          signature = c(query = "numeric",
                        target = "data.frameOrSimilar",
                        param = "Mz2MassParam"),
          function(query, target, param, mzColname = "mz") {
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
                        param = "Mz2MassParam"),
          function(query, target, param, mzColname = "mz") {
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
                        param = "Mz2MassParam"),
          function(query, target, param, mzColname = c("mz", "mz")) {
              if(length(mzColname) == 1)
                  mzColname <- rep(mzColname, 2)
              if (!mzColname[1] %in% colnames(query))
                  stop("Missing column \"", mzColname[1], "\" in query")
              if (!mzColname[2] %in% colnames(target))
                  stop("Missing column \"", mzColname[2], "\" in target")
              res <- matchMz(query[, mzColname[1]], target[, mzColname[2]],
                             param)
              res@query <- query
              res@target <- target
              res
          })

#' @rdname matchMz
setMethod("matchMz",
          signature = c(query = "data.frameOrSimilar",
                        target = "data.frameOrSimilar",
                        param = "Mz2MassRtParam"),
          function(query, target, param, mzColname = c("mz", "mz"),
                   rtColname = c("rt", "rt")) {
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
              query_mass <- .mz_to_mass_df(query[, mzColname[1]],
                                           param@queryAdducts)
              query_mass$rt <- rep(query[, rtColname[1]],
                                   .nelements(param@queryAdducts))
              target_mass<- .mz_to_mass_df(target[, mzColname[2]],
                                           param@targetAdducts)
              target_mass$rt <- rep(target[, rtColname[2]],
                                    .nelements(param@targetAdducts))
              queryl <- nrow(query_mass)
              matches <- vector("list", queryl)
              for (i in seq_len(queryl)) {
                  matches[[i]] <-
                      .getMatchesMzRt(query_mass$index[i],
                                      query_mass$mass[i],
                                      query_mass$rt[i],
                                      target = target_mass,
                                      tolerance = param@tolerance,
                                      ppm = param@ppm,
                                      toleranceRt = param@toleranceRt)
                  if (nrow(matches[[i]]))
                      matches[[i]]$query_adduct <- query_mass$adduct[i]
                  else matches[[i]]$query_adduct <- character()
              }
              matches <- do.call(rbind, matches)
              colnames(matches)[3] <- "target_adduct"
              Matched(query = query, target = target,
                      matches = matches[, c(1, 2, 7, 3, 4, 5, 6)],
                      metadata = list(param = param))
          })

#' @rdname matchMz
setMethod("matchMz",
          signature = c(query = "SummarizedExperiment",
                        target = "ANY",
                        param = "Param"),
          function(query, target, param, mzColname = "mz",
                   rtColname = c("rt", "rt")) {
              matches <- matchMz(data.frame(rowData(query)), target, param,
                                 mzColname , rtColname)@matches
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
    diffs <- queryMz - target[, 2]
    absdiffs <- abs(diffs)
    cls <- which(absdiffs <= (tolerance + ppm(queryMz, ppm)))
    if ("adduct" %in% colnames(target)){
        if (length(cls))
            data.frame(query_idx = queryIndex,
                       target_idx = target$index[cls],
                       adduct = target$adduct[cls],
                       score = diffs[cls],
                       ppm_error = absdiffs[cls] / target[cls, 2] * 10^6)
        else data.frame(query_idx = integer(),
                        target_idx = integer(),
                        adduct = character(),
                        score = numeric(),
                        ppm_error = numeric())
    } else {
        if (length(cls))
            data.frame(query_idx = queryIndex,
                       target_idx = target$index[cls],
                       score = diffs[cls],
                       ppm_error = absdiffs[cls] / target[cls, 2] * 10^6)
        else data.frame(query_idx = integer(),
                        target_idx = integer(),
                        score = numeric(),
                        ppm_error = numeric())
    }
}

#' @noRd
.getMatchesMzRt <- function(queryIndex, queryMz, queryRt, target, tolerance,
                            ppm, toleranceRt){
    diffs_rt <- queryRt - target$rt
    cls_rt <- which(abs(diffs_rt) <= toleranceRt)
    diffs <- queryMz - target[cls_rt, 2]
    absdiffs <- abs(diffs)
    cls <- which(absdiffs <= (tolerance + ppm(queryMz, ppm)))
    if ("adduct" %in% colnames(target)){
        if (length(cls))
            data.frame(query_idx = queryIndex,
                       target_idx = target$index[cls_rt[cls]],
                       adduct = target$adduct[cls_rt[cls]],
                       score = diffs[cls],
                       ppm_error = absdiffs[cls] /
                           target[cls_rt[cls], 2] * 10^6,
                       score_rt = diffs_rt[cls_rt[cls]])
        else data.frame(query_idx = integer(),
                        target_idx = integer(),
                        adduct = character(),
                        score = numeric(),
                        ppm_error = numeric(),
                        score_rt = numeric())
    } else {
        if (length(cls))
            data.frame(query_idx = queryIndex,
                       target_idx = target$index[cls_rt[cls]],
                       score = diffs[cls],
                       ppm_error = absdiffs[cls] /
                           target[cls_rt[cls], 2] * 10^6,
                       score_rt = diffs_rt[cls_rt[cls]])
        else data.frame(query_idx = integer(),
                        target_idx = integer(),
                        score = numeric(),
                        ppm_error = numeric(),
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
               mz = as.numeric(mz),
               adduct = rep(colnames(mz), each = length(x)))
}

#' creates a `data.frame` with columns `"index"`, `"adduct"` and `"mass"` from
#' provided m/z values `x` and adduct definition `adducts`. `"index"` is the
#' index of the m/z in the input vector `x`.
#'
#' @importFrom MetaboCoreUtils mz2mass
#'
#' @noRd
.mz_to_mass_df <- function(x, adducts) {
    mass <- mz2mass(x, adducts)
    data.frame(index = rep(seq_along(x), .nelements(adducts)),
               mass = as.numeric(mass),
               adduct = rep(colnames(mass), each = length(x)))
}

.valid_adduct <- function(adducts, name = "`adducts`") {
    msg <- NULL
    if (is(adducts, "data.frame")) {
        if(any(!c("mass_add", "mass_multi") %in% colnames(adducts)))
            msg <- paste0("Columns \"mass_add\" and \"mass_multi\" must be ",
                          "present when ", name, " is a data.frame")
    } else {
        if (!all(adducts %in% c(adductNames("positive"),
                                adductNames("negative"))))
            msg <- paste0("Unknown adducts in ", name, " please check",
                          " MetaboCoreUtils for valid adducts")
    }
    msg
}
