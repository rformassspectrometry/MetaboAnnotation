#' @importClassesFrom S4Vectors DataFrame
#' @importClassesFrom QFeatures QFeatures
setClassUnion("data.frameOrSimilar", c("data.frame", "matrix", "DataFrame",
                                       "SummarizedExperiment", "QFeatures"))
setClassUnion("adductClass", c("character", "data.frame"))
