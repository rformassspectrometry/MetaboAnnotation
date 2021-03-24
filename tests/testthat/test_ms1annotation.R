test_that("MS1 annotation works", {

  cmpds <- data.frame(
    name = c("Tryptophan", "Leucine", "Isoleucine"),
    formula = c("C11H12N2O2", "C6H13NO2", "C6H13NO2"),
    exactmass = c(204.089878, 131.094629, 131.094629)
  )

  adducts <- c("[M+H]+", "[M+Na]+")

  x <- data.frame(
    mz = sort(c(205.097154, 132.101905))
  )

  result <- annotateMz(x, cmpds, adducts, tolerance = 0.005, ppm = 0)

  expect_equal(length(result), 2L)
  expect_equal(nrow(result[[1]]), 2)
  expect_equal(ncol(result[[1]]), 5)
  expect_equal(nrow(result[[2]]), 1)
  expect_equal(ncol(result[[2]]), 5)

})
