#' @importClassesFrom S4Vectors DataFrame
setClassUnion("data.frameOrSimilar", c("data.frame", "matrix", "DataFrame", "SummarizedExperiment"))
setClassUnion("adductClass", c("character", "data.frame"))
