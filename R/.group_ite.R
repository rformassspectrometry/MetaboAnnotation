#' @title Iterate through a table of standard compounds to group them.
#' 
#' @description
#' The `.group_ite` groups rows of a data frame of standard compounds based on 
#' their similarity and a specified number of standards per group. 
#'
#' @param x `numeric` matrix with row names representing the compounds and 
#' columns representing different adducts.
#' 
#' @param nstd `numeric` number of maximum standards per group.
#' 
#' @param min_diff `numeric` Minimum difference for considering two values as 
#' distinct. Default is 2.
#'
#' @return `list` where each element is a vector of row names representing a 
#' group of standards.
#' 
#' @details
#' Users should be aware that because the function iterates through `x`, the 
#' compound at the bottom of the list is more complicated to group, and there's 
#' a possibility that some compounds will not be grouped with others. We advise 
#' increasing `nstd` and randomizing `x` to see if it can resolve this problem. 
#' 
#' @examples
#' x <- data.frame(
#'   row.names = c("Malic Acid", "Pyridoxic Acid", "Thiamine", "Uric acid",
#'                 "dUTP", "N-Formyl-L-methionine"),
#'   "adduct_1" = c(135.0288, 184.0604, 265.1118, 169.0356, 468.9809, 178.0532),
#'   "adduct_2" = c(157.0107, 206.0424, 287.0937, 191.0176, 490.9628, 200.0352)
#' )
#' 
#' x <- as.matrix(x)
#' ## Group standards with a maximum of 3 per group and a minimum difference
#' ## of 2.
#' result <- .group_ite(x, nstd = 3, min_diff = 2)
#' result
#' 
#' @noRd
#' 

.group_ite <- function(x, nstd, min_diff = 2) {   
  
  output <- vector("list")
  g <- 0
  
  while (nrow(x) > 1) {
    g <- g + 1
    i <- 1
    group <- row.names(x)[i]
    
    while(length(group) < nstd & i < nrow(x)) {
      i <- i + 1
      diff_table <- abs(outer(as.vector(x[group, ]), as.vector(x[i,]), "-"))
      
      if (all(diff_table > min_diff, na.rm= TRUE)) { 
        group <- c(group, row.names(x)[i])
      }
    } 
    x <- x[!(rownames(x) %in% group), , drop = FALSE] 
    output[[g]] <- group
  }
  if (nrow(x))  
    output[[g + 1]] <- row.names(x)
    
  output
} 
