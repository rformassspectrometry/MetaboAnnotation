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
    expect_equal(.group_standards_iteration(x, nstd = 3), results)
    expect_is(.group_standards_iteration(x, nstd = 3), "list")
})
