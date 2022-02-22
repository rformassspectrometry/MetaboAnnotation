#' @noRd
setClass("IonModeParam",
         slots = c(
           tolerance = "numeric",
           ppm = "numeric",
           queryAdducts = "adductClass",
           targetAdducts = "adductClass"),
         contains = "Param",
         prototype = prototype(
           tolerance = 0,
           ppm = 5,
           queryAdducts = "[M+H]+",
           targetAdducts = "[M-H]-"),
         validity = function(object) {
           msg <- NULL
           if (length(object@tolerance) != 1 || object@tolerance < 0) {
             msg <- c("'tolerance' has to be a positive number of length 1")
           }

           if (length(object@ppm) != 1 || object@ppm < 0) {
             msg <- c("'ppm' has to be a positive number of length 1")
           }

           if (is(object@queryAdducts, "data.frame")) {
             if(any(!c("mass_add", "mass_multi") %in% colnames(object@adducts)))
               msg <- paste0("Columns \"mass_add\" and \"mass_multi\" must be ",
                             "present when adducts is a data.frame")
           } else {
             if (!all(object@queryAdducts %in% c(adductNames("positive"),
                                            adductNames("negative"))))
               msg <- paste0("Unknown adducts, please check MetaboCoreUtils",
                             " for valid adducts")
           }

           if (is(object@targetAdducts, "data.frame")) {
             if(any(!c("mass_add", "mass_multi") %in% colnames(object@adducts)))
               msg <- paste0("Columns \"mass_add\" and \"mass_multi\" must be ",
                             "present when adducts is a data.frame")
           } else {
             if (!all(object@targetAdducts %in% c(adductNames("positive"),
                                                 adductNames("negative"))))
               msg <- paste0("Unknown adducts, please check MetaboCoreUtils",
                             " for valid adducts")
           }

           msg

         })

#' @title m/z matching
#'
#' @name matchIonMode
#'
#' @description
#'
#' The `matchIonMode` function is a prototype.
NULL

#' @rdname matchIonMode
#'
#' @export
setGeneric("matchIonMode", function(query, target, param, ...)
  standardGeneric("matchIonMode"))

#' @rdname matchIonMode
#'
#' @importFrom BiocParallel bpmapply SerialParam
setMethod("matchIonMode",
          signature = c(query = "numeric",
                        target = "numeric",
                        param = "IonModeParam"),
          function(query, target, param, BPPARAM = SerialParam()) {

            query_mass <- .mz_to_mass_df(query, param@queryAdducts)
            target_mass<- .mz_to_mass_df(target, param@targetAdducts)

            queryl <- nrow(query_mass)

            matches <- vector("list", queryl)

            for(i in seq_len(queryl)) {

              matches[[i]] <- .getMassMatches(query_mass$index[i],
                                              query_mass$mass[i],
                                              query_mass$adduct[i],
                                              target = target_mass,
                                              tolerance = param@tolerance,
                                              ppm = param@ppm)

            }

            Matched(query = query, target = target,
                    matches = do.call(rbind, matches),
                    metadata = list(param = param))

          })

#' @author Michael Witting
#'
#' @param queryIndex `integer(1)` with the index of the query.
#'
#' @param queryMass `numeric(1)` with the mass of the query.
#'
#' @param target `data.frame` with columns `"index"`, `"mass"` and  `"adduct"`.
#'
#' @noRd
.getMassMatches <- function(queryIndex, queryMass, queryAdduct, target, tolerance, ppm) {

  diffs <- abs(queryMass - target$mass)
  cls <- which(diffs <= (tolerance + ppm(queryMass, ppm)))

  if (length(cls)) {
    data.frame(query_idx = queryIndex,
               target_idx = target$index[cls],
               adduct = paste0(queryAdduct, " / ", target$adduct[cls]),
               score = diffs[cls])
  } else {
    data.frame(query_idx = integer(),
               target_idx = integer(),
               adduct = character(),
               score = numeric())
  }
}

#' creates a `data.frame` with columns `"index"`, `"adduct"` and `"mz"` from
#' provided mass values `x` and adduct definition `adducts`. `"index"` is the
#' index of the mass in the input vector `x`.
#'
#' @importFrom MetaboCoreUtils mz2mass
#'
#' @noRd
.mz_to_mass_df <- function(x, adducts) {
  mass <- mz2mass(x, adducts)
  data.frame(index = rep(seq_along(x), .nelements(adducts)),
             adduct = rep(colnames(mass), each = length(x)),
             mass = as.numeric(mass))
}
