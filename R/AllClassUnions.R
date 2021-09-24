#' @importClassesFrom S4Vectors DataFrame
setClassUnion("data.frameOrSimilar", c("data.frame", "matrix", "DataFrame"))
setClassUnion("adductClass", c("character", "data.frame"))
