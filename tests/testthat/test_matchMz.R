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

  result <- matchMz(x, cmpds, adducts, tolerance = 0.005, ppm = 0)

  expect_equal(length(result), 2L)
  expect_equal(nrow(result[[1]]), 2)
  expect_equal(ncol(result[[1]]), 5)
  expect_equal(nrow(result[[2]]), 1)
  expect_equal(ncol(result[[2]]), 5)

})

############################################################################
############################################################################
############################################################################
############################################################################

test_that("TargetMass2MzParam works", {
  res <- TargetMass2MzParam()
  expect_true(is(res, "TargetMass2MzParam"))
  
  expect_error(TargetMass2MzParam(tolerance = 1:3), "positive number")
  expect_error(TargetMass2MzParam(ppm = -4), "positive number")
})

test_that("matchMz,TargetMass2MzParam works", {
  
  cmpds <- data.frame(
    name = c("Tryptophan", "Leucine", "Isoleucine"),
    formula = c("C11H12N2O2", "C6H13NO2", "C6H13NO2"),
    exactmass = c(204.089878, 131.094629, 131.094629)
  )
  
  adducts <- c("[M+H]+", "[M+Na]+")
  
  x <- data.frame(
    mz = c(mass2mz(204.089878, "[M+H]+"), 
           mass2mz(131.094629, "[M+H]+"), 
           mass2mz(204.089878, "[M+Na]+")+1e-6)
  )
  library(MetaboCoreUtils)
  library(MsCoreUtils)
  library(reshape2)
  ionDf <- .createIonDf(cmpds, adducts)
  
  par <- TargetMass2MzParam(tolerance = 0, ppm = 20)
 
  res <- matchMz(x, cmpds, par, adducts)
  expect_equal(query(res), x) 
  expect_equal(target(res), cmpds)
  expect_equal(res@matches$query_idx, c(1, 2, 2, 3))
  expect_equal(res@matches$target_idx, c(1, 2, 3, 1))
  expect_equal(res@matches$score, c(0, 0, 0, 1e-6))
})
