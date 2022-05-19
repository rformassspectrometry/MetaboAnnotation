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

#' @rdname matchValues
#'
#' @importFrom methods new
#'
#' @export
ValueParam <- function(tolerance = 0, ppm = 5) {
    new("ValueParam", tolerance = tolerance, ppm = ppm)
}

#' @noRd
setClass("MzParam", contains = "ValueParam")

#' @rdname matchValues
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

#' @rdname matchValues
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
                 msg <- "'toleranceRt' has to be a positive number of length 1"
             msg
         })

#' @rdname matchValues
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
                 msg <- "'toleranceRt' has to be a positive number of length 1"
             msg
         })

#' @rdname matchValues
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

#' @rdname matchValues
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
                 msg <- "'toleranceRt' has to be a positive number of length 1"
             msg
         })

#' @rdname matchValues
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

#' @title Matching of numeric values
#'
#' @aliases matchMz
#'
#' @name matchValues
#'
#' @export matchMz
#'
#' @description
#'
#' The `matchValues` method matches elements from `query` with those in `target`
#' using different matching approaches depending on parameter `param`.
#' Generally, `query` is expected to contain MS experimental values
#' (m/z and possibly retention time) while `target` reference values. `query`
#' and `target` can be `numeric`, a two dimensional array (such as a
#' `data.frame`, `matrix` or `DataFrame`), a `SummarizedExperiment`
#' or a `QFeatures`. For `SummarizedExperiment`, the information for
#' the matching is expected to be in the object's `rowData`. For `QFeatures`
#' matching is performed for values present in the `rowData` of one of the
#' object's assays (which needs to be specified with the `assayQuery`
#' parameter - if a `QFeatures` is used as `target` the name of the assay needs
#' to be specified with parameter `assayTarget`). `matchMz` is an alias for
#' `matchValues` to allow backward compatibility.
#'
#' Available `param` objects and corresponding matching approaches are:
#'
#' - `ValueParam`: generic matching between values in `query` and `target` given
#'   acceptable differences expressed in `ppm` and `tolerance`. If `query` or
#'   `target` are not numeric, parameter `valueColname` has to be used to
#'   specify the name of the column that contains the values to be matched.
#'   The function returns a [Matched()] object.
#'
#' - `MzParam`: match query m/z values against reference compounds for which
#'   also m/z are known. Matching is performed similarly to the `ValueParam`
#'   above. If `query` or `target` are not numeric, the column name containing
#'   the values to be compared must be defined with `matchValues`' parameter
#'   `mzColname`, which defaults to `"mz"`. `MzParam` parameters `tolerance`
#'   and `ppm` allow to define the maximal acceptable (constant or m/z relative)
#'   difference between query and target m/z values.
#'
#' - `MzRtParam`: match m/z **and** retention time values between `query` and
#'   `target`. Parameters `mzColname` and `rtColname` of the `matchValues`
#'   function allow to define the columns in `query` and `target` containing
#'   these values (defaulting to `c("mz", "mz")` and `c("rt", "rt")`,
#'   respectively). `MzRtParam` parameters `tolerance` and
#'   `ppm` have the same meaning as in `MzParam`; `MzRtParam` parameter
#'   `toleranceRt` allows to specify the maximal acceptable difference between
#'   query and target retention time values.
#'
#' - `Mass2MzParam`: match m/z values against reference compounds for
#'   which only the (exact) mass is known. Before matching, m/z values are
#'   calculated from the compounds masses in the *target* table using the
#'   adducts specified via `Mass2MzParam` `adducts` parameter (defaults to
#'   `adducts = "[M+H]+"`). After conversion of adduct masses to m/z values,
#'   matching is performed similarly to `MzParam` (i.e. the same parameters
#'   `ppm` and `tolerance` can be used). If `query` is not `numeric`,
#'   parameter `mzColname` of `matchValues` can be used to specify the column
#'   containing the query's m/z values (defaults to `"mz"`). If `target` is a
#'   is not `numeric`, parameter `massColname` can be used to define the
#'   column containing the reference compound's masses (defaults to
#'   `"exactmass"`).
#'
#' - `Mass2MzRtParam`: match m/z **and** retention time values against
#'   reference compounds for which the (exact) mass **and** retention time are
#'   known. Before matching, exact masses in `target` are converted to m/z
#'   values as for `Mass2MzParam`. Matching is then performed similarly to
#'   `MzRtParam`, i.e. m/z and retention times of entities are compared. With
#'   `matchValues`' parameters `mzColname`, `rtColname` and `massColname` the
#'   columns containing m/z values (in `query`), retention time values (in
#'   `query` and `target`) and exact masses (in `target`) can be specified.
#'
#' - `Mz2MassParam`: input values for `query` and `target` are expected to be
#'   m/z values but matching is performed on exact masses calculated from these
#'   (based on the provided adduct definitions). In detail, m/z values in
#'   `query` are first converted to masses with the [mz2mass()] function based
#'   on the adducts defined with `queryAdducts` (defaults to `"[M+H]+"`). The
#'   same is done for m/z values in `target` (adducts can be defined with
#'   `targetAdducts` which defaults to `"[M-H-]"). Matching is then performed
#'   on these converted values similarly to `ValueParam`. If `query` or `target`
#'   are not numeric, the column containing the m/z values can be
#'   specified with `matchValues`' parameter `mzColname` (defaults to `"mz"`).
#'
#' - `Mz2MassRtParam`: same as `Mz2MassParam` but with additional comparison of
#'   retention times between `query` and `target`. Parameters `rtColname` and
#'   `mzColname` of `matchValues` allow to specify which columns contain the
#'   retention times and m/z values, respectively.
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
#'     assumed default name for columns with m/z values is `"mz"`.
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
#'     a `numeric`, `data.frame`, `DataFrame`, `matrix`, `SummarizedExperiment`
#'     or `QFeatures`. It is expected to contain m/z values and can contain also
#'     other variables. Matchings based on both m/z and retention time can be
#'     performed when a column with retention times is present in both `query`
#'     and `target`.
#'
#' @param queryAdducts for `Mz2MassParam`. Adducts used to derive mass
#'     values from query m/z values. The expected format is the same as that
#'     for parameter `adducts`.
#'
#' @param queryAssay `character(1)` specifying the name of the assay of the
#'     provided `QFeatures` that should be used for the matching (values from
#'     this assay's `rowData` will be used for matching). Only used if `query`
#'     is an instance of a `QFeatures` object.
#'
#' @param rtColname `character(2)` with the name of the column containing
#'     the compounds retention times in `query` and the name for the one in
#'     `target`. It can also be `character(1)` if the two names are the same.
#'     To be used when `param` is `MzRtParam` or `Mass2MzRtParam`.
#'     Defaults to `rtColname = c("rt", "rt")`.
#'
#' @param target compound table with metabolites to compare against. The
#'     expected types are the same as those for `query`.
#'
#' @param targetAdducts for `Mz2MassParam`. Adducts used to derive mass
#'     values from target m/z values. The expected format is the same as that
#'     for parameter `adducts`.
#'
#' @param targetAssay `character(1)` specifying the name of the assay of the
#'     provided `QFeatures` that should be used for the matching (values from
#'     this assay's `rowData` will be used for matching). Only used if `target`
#'     is an instance of a `QFeatures` object.
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
#' @return [Matched] object representing the result.
#'
#' Depending on the `param` object different *scores* representing the quality
#' of the match are provided. This comprises absolute as well as relative
#' differences (column/variables `"score"` and `"ppm_error"` respectively).
#' If `param` is a `Mz2MassParam`, `"score"` and `"ppm_error"` represent
#' differences of the compared masses (calculated from the provided m/z values).
#' If `param` an `MzParam`, `MzRtParam`, `Mass2MzParam` or `Mass2MzRtParam`,
#' `"score"` and `"ppm_error"` represent absolute and relative differences of
#' m/z values.
#' Additionally, if `param` is either an `MzRtParam` or `Mass2MzRtParam`
#' differences between query and target retention times for each matched
#' element is available in the column/variable `"score_rt"` in the returned
#' `Matched` object.
#' Negative values of `"score"` (or `"score_rt"`) indicate that the m/z or mass
#' (or retention time) of the query element is smaller than that of the target
#' element.
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
#' res <- matchValues(fts, target_df, parm)
#' res
#'
#' ## List the available variables/columns
#' colnames(res)
#'
#' ## feature_id and mz are from the query data frame, while target_name,
#' ## target_formula and target_exactmass are from the query object (columns
#' ## from the target object have a prefix *target_* added to the original
#' ## column name. Columns adduct, score and ppm_error represent the results
#' ## of the matching: adduct the adduct/ion of the original compound for which
#' ## the m/z matches, score the absolute difference of the query and target
#' ## m/z and ppm_error the relative difference in m/z values.
#'
#' ## Get the full matching result:
#' matchedData(res)
#'
#' ## We have thus matches of FT002 to two different compounds (but with the
#' ## same mass).
#'
#' ## Individual columns can also be accessed with the $ operator:
#' res$feature_id
#' res$target_name
#' res$ppm_error
#'
#'
#' ## We repeat the matching requiring an exact match
#' parm <- Mass2MzParam(
#'     adducts = c("[M+H]+", "[M+Na]+"),
#'     tolerance = 0,
#'     ppm = 0)
#' res <- matchValues(fts, target_df, parm)
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
#' res <- matchValues(fts, target_df, parm)
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
#' res <- matchValues(fts, target_df, parm)
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
#' res <- matchValues(fts, target_df, parm)
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
#' res <- matchValues(mz1, mz2, MzParam(tolerance = 0.001))
#'
#' matchedData(res)
#'
#' ## Matching with a SummarizedExperiment or a QFeatures work analogously,
#' ## only that the matching is performed on the object's `rowData`.
#'
#' ## Below we create a simple SummarizedExperiment with some random assay data.
#' ## Note that results from a data preprocessing with the `xcms` package could
#' ## be extracted as a `SummarizedExperiment` with the `quantify` method from
#' ## the `xcms` package.
#' library(SummarizedExperiment)
#' se <- SummarizedExperiment(
#'     assays = matrix(rnorm(12), nrow = 3, ncol = 4),
#'     rowData = fts)
#'
#' ## We can now perform the matching of this SummarizedExperiment against the
#' ## target_df as before.
#' res <- matchValues(se, target_df,
#'     param = Mass2MzParam(adducts = c("[M+H]+", "[M+Na]+"),
#'         tolerance = 0, ppm = 20))
#' res
#'
#' ## Getting the available columns
#' colnames(res)
#'
#' ## The query columns represent the columns of the object's `rowData`
#' rowData(se)
#'
#' ## matchedData also returns the query object's rowData along with the
#' ## matching entries in the target object.
#' matchedData(res)
#'
#' ## While `query` will return the full SummarizedExperiment.
#' query(res)
#'
#' ## To illustrate use with a QFeatures object we first create a simple
#' ## QFeatures object with two assays, `"ions"` representing the full feature
#' ## data.frame and `"compounds"` a subset of it.
#' library(QFeatures)
#' qf <- QFeatures(list(ions = se, compounds = se[2,]))
#'
#' ## We can perform the same matching as before, but need to specify which of
#' ## the assays in the QFeatures should be used for the matching. Below we
#' ## perform the matching using the "ions" assay.
#' res <- matchValues(qf, target_df, queryAssay = "ions",
#'     param = Mass2MzParam(adducts = c("[M+H]+", "[M+Na]+"),
#'         tolerance = 0, ppm = 20))
#' res
#'
#' ## colnames returns now the colnames of the `rowData` of the `"ions"` assay.
#' colnames(res)
#'
#' matchedData(res)
NULL

#' @rdname matchValues
#'
#' @export
setGeneric("matchValues", function(query, target, param, ...)
    standardGeneric("matchValues"))

matchMz <- matchValues

#' @rdname matchValues
setMethod("matchValues",
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

#' @rdname matchValues
setMethod("matchValues",
          signature = c(query = "numeric",
                        target = "data.frameOrSimilar",
                        param = "ValueParam"),
          function(query, target, param, valueColname = character(),
                   targetAssay = character()) {
              if(!length(valueColname))
                  stop("`valueColname` has to be provided.")
              target_ <- .objectToMatch(target, targetAssay, valueColname)
              res <- matchValues(query, target_, param)
              res@target <- target
              res@targetAssay <- targetAssay
              res
          })

#' @rdname matchValues
setMethod("matchValues",
          signature = c(query = "data.frameOrSimilar",
                        target = "numeric",
                        param = "ValueParam"),
          function(query, target, param, valueColname = character(),
                   queryAssay = character()) {
              if(!length(valueColname))
                  stop("`valueColname` has to be provided.")
              query_ <- .objectToMatch(query, queryAssay, valueColname)
              res <- matchValues(query_, target, param)
              res@query <- query
              res@queryAssay <- queryAssay
              res
          })

#' @rdname matchValues
setMethod("matchValues",
          signature = c(query = "data.frameOrSimilar",
                        target = "data.frameOrSimilar",
                        param = "ValueParam"),
          function(query, target, param, valueColname = character(),
                   queryAssay = character(), targetAssay = character()) {
              if(!length(valueColname))
                  stop("`valueColname` has to be provided.")
              if(length(valueColname) == 1)
                  valueColname <- rep(valueColname, 2)
              query_ <- .objectToMatch(query, queryAssay, valueColname[1])
              target_ <- .objectToMatch(target, targetAssay, valueColname[2])
              res <- matchValues(query_, target_, param)
              res@query <- query
              res@queryAssay <- queryAssay
              res@target <- target
              res@targetAssay <- targetAssay
              res
          })

#' @rdname matchValues
setMethod("matchValues",
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

#' @rdname matchValues
setMethod("matchValues",
          signature = c(query = "numeric",
                        target = "data.frameOrSimilar",
                        param = "Mass2MzParam"),
          function(query, target, param, massColname = "exactmass",
                   targetAssay = character()) {
              target_ <- .objectToMatch(target, targetAssay, massColname)
              res <- matchValues(query, target_, param)
              res@target <- target
              res@targetAssay <- targetAssay
              res
          })

#' @rdname matchValues
setMethod("matchValues",
          signature = c(query = "data.frameOrSimilar",
                        target = "numeric",
                        param = "Mass2MzParam"),
          function(query, target, param, mzColname = "mz",
                   queryAssay = character()) {
              query_ <- .objectToMatch(query, queryAssay, mzColname)
              res <- matchValues(query_, target, param)
              res@query <- query
              res@queryAssay <- queryAssay
              res
          })

#' @rdname matchValues
setMethod("matchValues",
          signature = c(query = "data.frameOrSimilar",
                        target = "data.frameOrSimilar",
                        param = "Mass2MzParam"),
          function(query, target, param, mzColname = "mz",
                   massColname = "exactmass", queryAssay = character(0),
                   targetAssay = character(0)) {
              query_ <- .objectToMatch(query, queryAssay, mzColname)
              target_ <- .objectToMatch(target, targetAssay, massColname)
              res <- matchValues(query_, target_, param)
              res@query <- query
              res@queryAssay <- queryAssay
              res@target <- target
              res@targetAssay <- targetAssay
              res
          })

#' @rdname matchValues
setMethod("matchValues",
          signature = c(query = "numeric",
                        target = "data.frameOrSimilar",
                        param = "MzParam"),
          function(query, target, param, mzColname = "mz",
                   targetAssay = character()) {
              target_ <- .objectToMatch(target, targetAssay, mzColname)
              res <- matchValues(query, target_, param)
              res@target <- target
              res@targetAssay <- targetAssay
              res
          })

#' @rdname matchValues
setMethod("matchValues",
          signature = c(query = "data.frameOrSimilar",
                        target = "numeric",
                        param = "MzParam"),
          function(query, target, param, mzColname = "mz",
                   queryAssay = character()) {
              query_ <- .objectToMatch(query, queryAssay, mzColname)
              res <- matchValues(query_, target, param)
              res@query <- query
              res@queryAssay <- queryAssay
              res
          })

#' @rdname matchValues
setMethod("matchValues",
          signature = c(query = "data.frameOrSimilar",
                        target = "data.frameOrSimilar",
                        param = "MzParam"),
          function(query, target, param, mzColname = c("mz", "mz"),
                   queryAssay = character(), targetAssay = character()) {
              if(length(mzColname) == 1)
                  mzColname <- rep(mzColname, 2)
              query_ <- .objectToMatch(query, queryAssay, mzColname[1])
              target_ <- .objectToMatch(target, targetAssay, mzColname[2])
              res <- matchValues(query_, target_, param)
              res@query <- query
              res@queryAssay <- queryAssay
              res@target <- target
              res@targetAssay <- targetAssay
              res
          })

#' @rdname matchValues
setMethod("matchValues",
          signature = c(query = "data.frameOrSimilar",
                        target = "data.frameOrSimilar",
                        param = "Mass2MzRtParam"),
          function(query, target, param, massColname = "exactmass",
                   mzColname = "mz", rtColname = c("rt", "rt"),
                   queryAssay = character(), targetAssay = character()) {
              if(length(rtColname) == 1)
                  rtColname <- rep(rtColname, 2)
              query_ <- .objectToMatch(query, queryAssay,
                                       c(mzColname, rtColname[1]))
              target_ <- .objectToMatch(target, targetAssay,
                                        c(massColname, rtColname[2]))
              target_mz <- .mass_to_mz_df(target_[, massColname],
                                          param@targetAdducts)
              target_mz$rt <- rep(target_[, rtColname[2]],
                                  .nelements(param@targetAdducts))
              queryl <- nrow(query_)
              matches <- vector("list", queryl)
              query_mz <- query_[, mzColname]
              query_rt <- query_[, rtColname[1L]]
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
                      queryAssay = queryAssay, targetAssay = targetAssay,
                      metadata = list(param = param))
          })

#' @rdname matchValues
setMethod("matchValues",
          signature = c(query = "data.frameOrSimilar",
                        target = "data.frameOrSimilar",
                        param = "MzRtParam"),
          function(query, target, param, mzColname = c("mz", "mz"),
                   rtColname = c("rt", "rt"), queryAssay = character(),
                   targetAssay = character()) {
              if(length(mzColname) == 1)
                  mzColname <- rep(mzColname, 2)
              if(length(rtColname) == 1)
                  rtColname <- rep(rtColname, 2)
              query_ <- .objectToMatch(query, queryAssay,
                                       c(mzColname[1], rtColname[1]))
              target_ <- .objectToMatch(target, targetAssay,
                                        c(mzColname[2], rtColname[2]))
              target_mz <- data.frame(index = seq_len(nrow(target_)),
                                      mz = target_[, 1], rt = target_[, 2])
              queryl <- nrow(query_)
              matches <- vector("list", queryl)
              query_mz <- query_[, 1]
              query_rt <- query_[, 2]
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
                      queryAssay = queryAssay, targetAssay = targetAssay,
                      metadata = list(param = param))
          })

#' @rdname matchValues
setMethod("matchValues",
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

#' @rdname matchValues
setMethod("matchValues",
          signature = c(query = "numeric",
                        target = "data.frameOrSimilar",
                        param = "Mz2MassParam"),
          function(query, target, param, mzColname = "mz",
                   targetAssay = character()) {
              target_ <- .objectToMatch(target, targetAssay, mzColname)
              res <- matchValues(query, target_, param)
              res@target <- target
              res@targetAssay <- targetAssay
              res
          })

#' @rdname matchValues
setMethod("matchValues",
          signature = c(query = "data.frameOrSimilar",
                        target = "numeric",
                        param = "Mz2MassParam"),
          function(query, target, param, mzColname = "mz",
                   queryAssay = character()) {
              query_ <- .objectToMatch(query, queryAssay, mzColname)
              res <- matchValues(query_, target, param)
              res@query <- query
              res@queryAssay <- queryAssay
              res
          })
#' @rdname matchValues
setMethod("matchValues",
          signature = c(query = "data.frameOrSimilar",
                        target = "data.frameOrSimilar",
                        param = "Mz2MassParam"),
          function(query, target, param, mzColname = c("mz", "mz"),
                   queryAssay = character(), targetAssay = character()) {
              if(length(mzColname) == 1)
                  mzColname <- rep(mzColname, 2)
              query_ <- .objectToMatch(query, queryAssay, mzColname[1])
              target_ <- .objectToMatch(target, targetAssay, mzColname[2])
              res <- matchValues(query_, target_, param)
              res@query <- query
              res@queryAssay <- queryAssay
              res@target <- target
              res@targetAssay <- targetAssay
              res
          })

#' @rdname matchValues
setMethod("matchValues",
          signature = c(query = "data.frameOrSimilar",
                        target = "data.frameOrSimilar",
                        param = "Mz2MassRtParam"),
          function(query, target, param, mzColname = c("mz", "mz"),
                   rtColname = c("rt", "rt"), queryAssay = character(),
                   targetAssay = character()) {
              if(length(mzColname) == 1)
                  mzColname <- rep(mzColname, 2)
              if(length(rtColname) == 1)
                  rtColname <- rep(rtColname, 2)
              query_ <- .objectToMatch(query, queryAssay,
                                       c(mzColname[1], rtColname[1]))
              target_ <- .objectToMatch(target, targetAssay,
                                        c(mzColname[2], rtColname[2 ]))
              query_mass <- .mz_to_mass_df(query_[, mzColname[1]],
                                           param@queryAdducts)
              query_mass$rt <- rep(query_[, rtColname[1]],
                                   .nelements(param@queryAdducts))
              target_mass<- .mz_to_mass_df(target_[, mzColname[2]],
                                           param@targetAdducts)
              target_mass$rt <- rep(target_[, rtColname[2]],
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
                      queryAssay = queryAssay, targetAssay = targetAssay,
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
