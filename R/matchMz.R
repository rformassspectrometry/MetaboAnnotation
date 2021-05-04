#' @title m/z matching
#'
#' @name matchMz
#'
#' @description
#'
#' The `matchMz` method matches (compares) m/z values from a MS1 data table
#'    with theoretical m/z values from compounds. `adducts` specifies the allowed
#'    adducts based on adduct names from [MetaboCoreUtils::adductNames()]. The
#'    maximum allowed error is defined by `tolerance` or `ppm`
#'
#' @param x Feature table containing information on MS1 features. Must contain
#'    column called `mz`
#'
#' @param y Compound table with metabolites to compare against. Must contain
#'    columns called `exactmass` and `name`
#'
#' @param adducts Allowed adducts named accordingly to [MetaboCoreUtils::adductNames()]
#'
#' @param tolerance `numeric(1)` for an absolute maximal accepted difference
#'   between m/z values.
#'
#' @param ppm `numeric(1)` for a relative, m/z-dependent, maximal accepted
#'   difference between m/z values.
#'
#' @return `list` of the same length as x with data frames containing the
#'   annotations
#'
#' @author Michael Witting
#'
#' @export
setGeneric("matchMz", function(x, y, adducts, tolerance, ppm)
  standardGeneric("matchMz")
)


#' @rdname matchMz
#'
#' @importFrom MetaboCoreUtils adductNames
#'
#' @export
setMethod("matchMz",
          signature = c(x = "data.frame", y = "data.frame"),
          function(x, y, adducts = c("[M+H]+"), tolerance = 0, ppm = 0) {

            # some sanity checks
            if(!any(adducts %in% c(adductNames("positive"), adductNames("negative")))) {

              stop("Unknown adducts, please check MetaboCoreUtils for valid adducts")

            }

            if(!"mz" %in% colnames(x)) {

              stop("Missing mz column in x")

            }

            cmpds <- y

            if(!all(c("name", "exactmass") %in% colnames(cmpds))) {

              stop("Missing name and exactmass column in y")

            }

            mz <- x$mz

            .matchMz(mz, cmpds, adducts, tolerance = tolerance, ppm = ppm)

          })


#' @rdname matchMz
#'
#' @importFrom MetaboCoreUtils adductNames
#'
#' @importClassesFrom CompoundDb CompDb
#'
#' @importMethodsFrom CompoundDb compounds
#'
#' @export
setMethod("matchMz", signature = c(x = "data.frame", y = "CompDb"),
          function(x, y, adducts = c("[M+H]+"), tolerance = 0, ppm = 0) {

            # some sanity checks
            if(!any(adducts %in% c(adductNames("positive"), adductNames("negative")))) {

              stop("Unknown adducts, please check MetaboCoreUtils for valid adducts")

            }

            if(!"mz" %in% colnames(x)) {

              stop("Missing mz column in x")

            }

            cmpds <- compounds(y)

            cmpds <- y

            if(!all(c("name", "exactmass") %in% colnames(cmpds))) {

              stop("Missing name and exactmass column in y")

            }

            mz <- x$mz

            .matchMz(mz, cmpds, adducts, tolerance = tolerance, ppm = ppm)

          })


.matchMz <- function(mz,
                     cmpds,
                     adducts,
                     tolerance = 0,
                     ppm = 0) {
  
  # create data.frame will all ion m/z
  ionDf <- .createIonDf(cmpds, adducts)
  
  # perform annotation of m/z values
  lapply(mz, .closestFeature, ionDf = ionDf, tolerance = tolerance, ppm = ppm)
  
}


#' @importFrom MsCoreUtils closest
#'
#' @noRd
.closestFeature <- function(mz,
                            ionDf,
                            tolerance = 0,
                            ppm = 0) {
  
  # perform search based on mz
  matches <- closest(ionDf$mz,
                     mz,
                     tolerance = tolerance,
                     ppm = ppm) == 1
  
  matches[is.na(matches)] <- FALSE
  
  if(any(matches)) {
    
    return(ionDf[matches,])
    
  } else {
    
    return(NA)
    
  }
}

############################################################################
############################################################################
############################################################################
############################################################################


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
             msg <- "Unknown adducts, please check MetaboCoreUtils for valid adducts"
           msg
         })



#' @rdname TargetMass2MzParam
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
#'    with theoretical m/z values from compounds. `adducts` specifies the allowed
#'    adducts based on adduct names from [MetaboCoreUtils::adductNames()]. The
#'    maximum allowed error is defined by `tolerance` or `ppm` in `param`
#'
#' @param query Feature table containing information on MS1 features. Must contain
#'    column called `mz`
#'
#' @param target Compound table with metabolites to compare against. Must contain
#'    columns called `exactmass` and `name`
#'
#' @param param object containing information on adducts (named accordingly to 
#' [MetaboCoreUtils::adductNames()]), tolerance (`numeric(1)` for an absolute 
#' maximal accepted difference between m/z values), ppm (`numeric(1)` for a 
#' relative, m/z-dependent, maximal accepted difference between m/z values).
#'
#' @return `Matched` object representing the annotations
#'
#' @author Andrea Vicini, Michael Witting
#'
#' @export

setGeneric("matchMz", function(query, target, param, ...)
  standardGeneric("matchMz")
)


setMethod("matchMz", 
          signature = c(query = "data.frame", 
                        target = "data.frame", 
                        param = "TargetMass2MzParam"),
          function(query, target, param, BPPARAM = SerialParam()){
            
            if(!"mz" %in% colnames(query))  stop("Missing mz column in query")
            if(!all(c("name", "exactmass") %in% colnames(target)))
              stop("Missing name and exactmass column in target")
              
            mz <- query$mz
            ionDf <- .createIonDf(target, param@adducts)
            matches <- do.call(rbind, bplapply(seq_along(mz), .getMatches, mz, 
                                               ionDf, param@tolerance, param@ppm, 
                                               BPPARAM = BPPARAM))
            Matched(query = query, target = target, matches = matches)
          })


#' @author Andrea Vicini
#'
#' @noRd
.getMatches <- function(index, mzquery, ionDf, tolerance, ppm){
  mz <- mzquery[index]
  cls <- closest(ionDf$mz,
                 mz,
                 tolerance = tolerance,
                 ppm = ppm)
  
  mtchd <- which(!is.na(cls))
  matches <- data.frame(query_idx = rep(index, length(mtchd)),
                        target_idx = (as.integer(rownames(ionDf))
                                      [mtchd] - 1L) %% nrow(cmpds) + 1L,
                        adduct = ionDf$adduct[mtchd],
                        score = mz - ionDf$mz[mtchd])
  matches
}


#' @importFrom MetaboCoreUtils mass2mz
#'
#' @importFrom reshape2 melt
#'
#' @noRd
.createIonDf <- function(cmpds,
                         adducts) {

  adductMzs <- mass2mz(cmpds$exactmass, adducts)

  ionDf <- melt(cbind.data.frame(cmpds, adductMzs, stringsAsFactors = FALSE),
                id.vars = colnames(cmpds),
                value.name = "mz",
                variable.name = "adduct")

  ionDf$adduct <- as.character(ionDf$adduct)

  ionDf[order(ionDf$mz),]

}
