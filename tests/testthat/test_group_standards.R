#' define some variable for test 
x <- data.frame(
  row.names = c("Malic Acid", "Pyridoxic Acid", "Thiamine", "Uric acid",
                "dUTP", "N-Formyl-L-methionine"),
  "adduct_1" = c(135.0288, 184.0604, 265.1118, 169.0356, 468.9809, 178.0532),
  "adduct_2" = c(157.0107, 206.0424, 287.0937, 191.0176, 490.9628, 200.0352)
)
x <- as.matrix(x)

#' Expected list to get for maximum of 3 standard per group 
results <- list(c("Malic Acid", "Pyridoxic Acid", "Thiamine"), 
                c("Uric acid", "dUTP", "N-Formyl-L-methionine"))

test_that("essentially that .group_standards_iteration works", {
    expect_equal(.group_standards_iteration(x, max_nstd = 3), results)
    expect_is(.group_standards_iteration(x, max_nstd = 3), "list")
})


#' test for randomization
set.seed(123) 
#' Create a matrix with compound names and ion masses
x <- matrix(c(349.0544, 371.0363, 325.0431, 347.0251, 581.0416, 603.0235,
              167.0564, 189.0383, 150.0583, 172.0403, 171.0053, 192.9872,
              130.0863, 152.0682, 768.1225, 790.1044),
            ncol = 2, byrow = TRUE)

rownames(x) <- c("IMP", "UMP", "UDP-glucuronate", "1-Methylxanthine",
                 "Methionine", "Dihydroxyacetone phosphate", "Pipecolic acid",
                 "CoA")
colnames(x) <- c("[M+H]+", "[M+Na]+")

#' run using `.group_standards_iteration`
standard_groups <- .group_standards_iteration(x, max_nstd = 4, min_diff = 2)

#' run using `.randomisation_group`.
standard_groups_r <- .randomize_grouping(x, max_nstd = 4,
                                         min_nstd = 3,
                                         min_diff = 2)
min_nstd <- 3
test_that("randomization improves results", {
  expect_true(length(standard_groups) > length(standard_groups_r))
  expect_true(any(lengths(standard_groups_r) < min_nstd))
  expect_false(length(standard_groups_r) == 0 )
}
  
)
