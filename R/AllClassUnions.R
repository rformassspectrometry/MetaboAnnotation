#' @importClassesFrom S4Vectors DataFrame
#' @importClassesFrom QFeatures QFeatures
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setClassUnion("data.frameOrSimilar", c("data.frame", "matrix", "DataFrame",
                                       "SummarizedExperiment", "QFeatures"))
setClassUnion("adductClass", c("character", "data.frame"))
setClassUnion("characterOrNumeric", c("character", "numeric"))
