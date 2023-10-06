#' @title Iterate through a table of standard compounds to group them.
#' 
#' @description
#' The `.group_standards_iteration` function groups rows of a matrix of standard 
#' compounds based on their similarity and a user-defined specified number of 
#' standards per group. The `.randomize_grouping` function utilize the 
#' `.group_standards_iteration` but tries for all combination of rows of `x`. 
#' This is great if the user is not satisfied by the results of the
#' `.group_standards_iteration` function.
#'
#' @param x `numeric` matrix with row names representing the compounds and 
#' columns representing different adducts.
#' 
#' @param max_nstd `numeric` number of maximum of standards per group.
#' 
#' @param min_nstd `numeric` number of minimum of standards per group.
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
#' to test the `.randomize_grouping`  function if that happens.  
#' 
#' @examples
#' 
#' ## `.group_standards_iteration` only
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
#' result <- .group_ite(x, max_nstd = 3, min_diff = 2)
#' result
#'
#'
#'
#' ## Comparing results with using `.randomize_grouping`.
#' set.seed(123)
#' ## Create a matrix with compound names and ion masses
#' x <- matrix(c(349.0544, 371.0363, 325.0431, 347.0251, 581.0416, 603.0235,
#'               167.0564, 189.0383, 150.0583, 172.0403, 171.0053, 192.9872,
#'               130.0863, 152.0682, 768.1225, 790.1044),
#'             ncol = 2, byrow = TRUE)
#'
#' rownames(x) <- c("IMP", "UMP", "UDP-glucuronate", "1-Methylxanthine", 
#'                  "Methionine", "Dihydroxyacetone phosphate", "Pipecolic acid", 
#'                  "CoA")
#' colnames(x) <- c("[M+H]+", "[M+Na]+")
#' 
#' ## run using `.group_standards_iteration` 
#' standard_groups <- .group_standards_iteration(x, max_nstd = 4, min_diff = 2)
#' standard_groups
#' 
#' ## get incomplete groups, rescue this using the `.randomize_grouping`. 
#' standard_groups_r <- .randomize_grouping(x, max_nstd = 4, 
#'                                          min_nstd = 3,
#'                                          min_diff = 2)
#' standard_groups_r 
#' 
#' 
#' @author Philippine Louail
#' 
#' @noRd
#'

.group_standards_iteration <- function(x, max_nstd, min_diff = 2) {   
    output <- vector("list")
    g <- 0
    
    while (nrow(x) > 1) {
        g <- g + 1
        i <- 1
        group <- row.names(x)[i]
        
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


.randomize_grouping <- function(x, 
                                max_nstd, 
                                min_nstd, 
                                min_diff = 2) { 
    n <- nrow(x)
    standard_groups <- vector("list")
    i <- 0
    while (length(standard_groups) == 0 || any(lengths(standard_groups) < 
                                              min_nstd)) {
        
        i <- i +1
        x <- x[sample(n), ]
        standard_groups <- .group_standards_iteration(x, 
                                                      max_nstd = max_nstd,
                                                      min_diff = min_diff)
        if (i > n*n)
            stop("all combination were tested, no possibility to fit your input requirement")
    }  
    standard_groups
}
