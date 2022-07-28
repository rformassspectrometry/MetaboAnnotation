#' @rdname Matched
#'
#' @exportMethod addMatches
setGeneric("addMatches", function(object, ...)
  standardGeneric("addMatches"))

#' @rdname Matched
#'
#' @exportMethod filterMatches
setGeneric("filterMatches", function(object, param, ...)
  standardGeneric("filterMatches"))

#' @rdname Matched
#'
#' @exportMethod matchedData
setGeneric("matchedData", function(object, ...)
    standardGeneric("matchedData"))

#' @rdname Matched
#'
#' @exportMethod endoapply
setGeneric("endoapply", function(object, ...)
    standardGeneric("endoapply"))