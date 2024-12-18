#' @title Create Standard Mixes from a Matrix of Standard Compounds
#'
#' @description
#'
#' The `createStandardMixes` function defines groups (mixes) of compounds
#' (standards) with dissimilar m/z values. The expected size of the groups can
#' be defined with parameters `max_nstd` and `min_nstd` and the minimum required
#' difference between m/z values within each group with parameter `min_diff`.
#' The group assignment will be reported in an additional column in the result
#' data frame.
#'
#' @param x `numeric` matrix with row names representing the compounds and
#' columns representing different adducts. Such a matrix with m/z values for
#' different adducts for compounds could e.g. be created with the
#' [MetaboCoreUtils::mass2mz()] function.
#'
#' @param max_nstd `numeric` number of maximum standards per group.
#'
#' @param min_nstd `numeric` number of minimum standards per group. Only
#' needed when using `iterativeRandomization = TRUE`.
#'
#' @param min_diff `numeric` Minimum difference for considering two values as
#' distinct.
#'
#' @param iterativeRandomization `logical` default `FALSE`. If set to `TRUE`,
#' `createStandardMixes` will randomly rearrange the rows of `x` until the user
#' inputs are satisfied.
#'
#' @return `data.frame` created by adding a column `group` to the input `x`
#' matrix, comprising the group number for each compound.
#'
#' @details
#' Users should be aware that because the function iterates through `x`, the
#' compounds at the bottom of the matrix are more complicated to group, and
#' there is a possibility that some compounds will not be grouped with others.
#' We advise specifyiong `iterativeRandomization = TRUE` even if it takes more
#' time.
#'
#' @examples
#'
#' ## Iterative grouping only
#' x <- matrix(c(135.0288, 157.0107, 184.0604, 206.0424, 265.1118, 287.0937,
#'               169.0356, 191.0176, 468.9809, 490.9628, 178.0532, 200.0352),
#'             ncol = 2, byrow = TRUE,
#'             dimnames = list(c("Malic Acid", "Pyridoxic Acid", "Thiamine",
#'                                 "Uric acid", "dUTP", "N-Formyl-L-methionine"),
#'                              c("adduct_1", "adduct_2")))
#' result <- createStandardMixes(x, max_nstd = 3, min_diff = 2)
#'
#' ## Randomize grouping
#' set.seed(123)
#' x <- matrix(c(349.0544, 371.0363, 325.0431, 347.0251, 581.0416, 603.0235,
#'               167.0564, 189.0383, 150.0583, 172.0403, 171.0053, 192.9872,
#'               130.0863, 152.0682, 768.1225, 790.1044),
#'             ncol = 2, byrow = TRUE,
#'             dimnames = list(c("IMP", "UMP", "UDP-glucuronate",
#'                                 "1-Methylxanthine", "Methionine",
#'                                 "Dihydroxyacetone phosphate",
#'                                 "Pipecolic acid", "CoA"),
#'                              c("[M+H]+", "[M+Na]+")))
#' result <- createStandardMixes(x, max_nstd = 4, min_nstd = 3, min_diff = 2,
#'                                iterativeRandomization = TRUE)
#'
#' @author Philippine Louail
#'
#' @export
#'
#' @rdname createStandardMixes
createStandardMixes <- function(x, max_nstd = 10,
                                min_nstd =5, min_diff = 2,
                                iterativeRandomization = FALSE) {
    if (!is.matrix(x))
        stop("The input `x` needs to be a matrix, try as.matrix(x)")
    if (!iterativeRandomization)
        output <- .iterative_grouping(x = x, max_nstd = max_nstd,
                                      min_diff = min_diff)
    else
        output <- .randomize_grouping(x = x, max_nstd = max_nstd,
                                      min_nstd = min_nstd, min_diff = min_diff)

    find_index <- function(name) {
        which(vapply(output, function(X) any(X %in% name), logical(1)))
        }
    x <- as.data.frame(x)
    x$group <- vapply(rownames(x), find_index, integer(1))

    x
}

#'
#' @noRd
#'

.iterative_grouping <- function(x, max_nstd, min_diff) {
    output <- vector("list")
    g <- 0

    while (nrow(x) > 1) {
        g <- g + 1
        i <- 1
        group <- row.names(x)[i]

        while (length(group) < max_nstd & i < nrow(x)) {
            i <- i + 1
            diff_table <- abs(outer(as.vector(x[group, ]), as.vector(x[i,]),
                                    "-"))

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

#'
#' @noRd
#'
.randomize_grouping <- function(x,
                                max_nstd,
                                min_nstd,
                                min_diff) {
    n <- nrow(x)
    output <- vector("list")
    i <- 0
    while (length(output) == 0 || any(lengths(output) < min_nstd)) {

        i <- i +1
        x <- x[sample(n), , drop = FALSE]
        output <- .iterative_grouping(x, max_nstd = max_nstd,
                                      min_diff = min_diff)
        if (i > n*n)
            stop("all combination were tested, no possibility to fit your",
                 " input requirements")
    }
    output
}
