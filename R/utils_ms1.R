#' @title Perform annotation on MS1 level
#'
#' @description
#'
#' `annotateMz` calculates the m/z value from a neutral mass and an adduct
#' definition.
#'
#' @param mz `numeric` neutral mass for which the adduct m/z shall be calculated.
#'
#' @param compdb `data.frame` data frame containing compounds used for annoation.
#'     Should minimally contain $name, $exactmass as input.
#'
#' @param adducts `character` specifying the name of the adduct;
#'     supported values are returned by [MetaboCoreUtils::adductNames()]
#'
#' @param tolerance `numeric` absolute annotation error
#'
#' @param ppm `numeric` ppm annotation error
#'
#' @return `list` List of data frames with the same length as mz.
#'
#' @author Michael Witting
#'
#' @export
annotateMz <- function(mz,
                       compdb,
                       adducts,
                       tolerance = 0,
                       ppm = 0) {

  # create data.frame will all ion m/z
  ionDf <- .createIonDf(compdb, adducts)

  # perform annotation of m/z values
  lapply(mz, .closestFeature, ionDf = ionDf, tolerance = tolerance, ppm = ppm)

}


#'
#'
#' @importFrom MetaboCoreUtils mass2mz
#' @importFrom reshape2 melt
.createIonDf <- function(compdb, adducts) {

  adductMzs <- mass2mz(compdb$exactmass, adducts)

  ionDf <- melt(cbind.data.frame(compdb, adductMzs),
                id.vars = colnames(compdb),
                value.name = "mz",
                variable.name = "adduct")

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

