forwardDp <- function(x, y, ...) {

  compareSpectra(x,
                 y,
                 FUN = .mydotproduct,
                 tolerance = tolerance,
                 ppm = ppm,
                 type = "outer")

}

backwardDp <- function(x, y, ...) {

  compareSpectra(x,
                 y,
                 FUN = .mydotproduct,
                 tolerance = tolerance,
                 ppm = ppm,
                 type = "right")

}

matchingPeaks <- function(x, y, ...) {

  compareSpectra(x,
                 y,
                 FUN = .specNrow,
                 tolerance = tolerance,
                 ppm = ppm,
                 type = "inner")

}


# costum matching function
.specNrow <-  function(x, y, ...) {
  nrow(x)
}

.mydotproduct <- function(x, y, type = NA, ...) {

  MsCoreUtils::ndotproduct(x, y, ...)

}
