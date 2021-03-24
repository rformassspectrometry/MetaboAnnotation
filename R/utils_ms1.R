
setGeneric("annotateMz", function(x, y, adducts, tolerance, ppm)
  standardGeneric("annotateMz")
  )

#'
#'
#' @importFrom MetaboCoreUtils adductNames
#'
#' @exportMethod annotateMz
setMethod("annotateMz", signature = c(x = "data.frame", y = "data.frame"),
          function(x, y, adducts = c("[M+H]+"), tolerance = 0, ppm = 0) {

            # some sanity checks
            if(!any(adducts %in% c(adductNames("positive"), adductNames("negative")))) {

              stop("Unknown adducts, please check MetaboCoreUtils for valid adducts")

            }

            if(!"mz" %in% colnames(x)) {

              stop("Missing mz column in x")

            }

            cmpds <- y
            mz <- x$mz

            .annotateMz(mz, cmpds, adducts, tolerance = tolerance, ppm = ppm)

          })

#'
#'
#' @importFrom MetaboCoreUtils adductNames
#' @importClassesFrom CompoundDb CompDb
#'
#' @exportMethod annotateMz
setMethod("annotateMz", signature = c(x = "data.frame", y = "CompDb"),
          function(x, y, adducts = c("[M+H]+"), tolerance = 0, ppm = 0) {

            # some sanity checks
            if(!any(adducts %in% c(adductNames("positive"), adductNames("negative")))) {

              stop("Unknown adducts, please check MetaboCoreUtils for valid adducts")

            }

            if(!"mz" %in% colnames(x)) {

              stop("Missing mz column in x")

            }

            cmpds <- compounds(y)
            mz <- x$mz

            .annotateMz(mz, cmpds, adducts, tolerance = tolerance, ppm = ppm)

          })



#'
#'
.annotateMz <- function(mz,
                        cmpds,
                        adducts,
                        tolerance = 0,
                        ppm = 0) {

  # create data.frame will all ion m/z
  ionDf <- .createIonDf(cmpds, adducts)

  # perform annotation of m/z values
  lapply(mz, .closestFeature, ionDf = ionDf, tolerance = tolerance, ppm = ppm)

}


#'
#'
#' @importFrom MetaboCoreUtils mass2mz
#' @importFrom reshape2 melt
.createIonDf <- function(cmpds, adducts) {

  adductMzs <- mass2mz(cmpds$exactmass, adducts)

  ionDf <- melt(cbind.data.frame(cmpds, adductMzs, stringsAsFactors = FALSE),
                id.vars = colnames(cmpds),
                value.name = "mz",
                variable.name = "adduct")

  ionDf$adduct <- as.character(ionDf$adduct)

  ionDf[order(ionDf$mz),]

}

#'
#'
#' @importFrom MsCoreUtils closest
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

