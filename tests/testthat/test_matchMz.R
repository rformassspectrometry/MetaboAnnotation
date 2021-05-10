test_that("TargetMass2MzParam works", {
  res <- TargetMass2MzParam()
  expect_true(is(res, "TargetMass2MzParam"))

  expect_error(TargetMass2MzParam(tolerance = 1:3), "positive number")
  expect_error(TargetMass2MzParam(ppm = -4), "positive number")
  expect_error(TargetMass2MzParam(adducts = c("[M+H]+", "adduct2")), "Unknown")
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

  par <- TargetMass2MzParam(adducts = adducts, tolerance = 0, ppm = 20)
  res <- matchMz(x, cmpds, par)
  expect_equal(query(res), x)
  expect_equal(target(res), cmpds)
  expect_equal(res@matches$query_idx, c(1, 2, 2, 3))
  expect_equal(res@matches$target_idx, c(1, 2, 3, 1))
  expect_equal(res@matches$score, c(0, 0, 0, 1e-6))

  par <- TargetMass2MzParam(adducts = adducts, tolerance = 0, ppm = 0)
  res <- matchMz(x, cmpds, par)
  expect_equal(query(res), x)
  expect_equal(target(res), cmpds)
  expect_equal(res@matches$query_idx, c(1, 2, 2))
  expect_equal(res@matches$target_idx, c(1, 2, 3))
  expect_equal(res@matches$score, c(0, 0, 0))

  ## no matches
  adducts <- c("[M+Li]+", "[M+K]+")
  par <- TargetMass2MzParam(adducts = adducts, tolerance = 0, ppm = 20)
  res <- matchMz(x, cmpds, par)
  expect_true(is(res, "Matched"))
  expect_equal(query(res), x)
  expect_equal(target(res), cmpds)
  expect_true(nrow(res@matches) == 0)
})

test_that(".getMatches works", {
    trgt <- data.frame(index = c(1, 2, 3, 1, 2, 3),
                       adduct = c("A", "A", "A", "B", "B", "B"),
                       mz = c(1, 1.2, 1.21, 1.2, 2, 2.1))
    trgt <- trgt[order(trgt$mz), ]

    res <- .getMatches(3, 1.2, trgt, tolerance = 0, ppm = 0)

    expect_true(is.data.frame(res))
    expect_equal(res$query_idx, c(3, 3))
    expect_equal(res$target_idx, c(2, 1))

    res <- .getMatches(4, 4, trgt, tolerance = 0, ppm = 0)
    expect_true(is.data.frame(res))
    expect_true(nrow(res) == 0)

    res <- .getMatches(3, 1.2, trgt, tolerance = 0.1, ppm = 0)
    expect_equal(res$target_idx, c(2, 1, 3))
})

test_that(".mass_to_mz_df works", {
    mass <- c(4, 2, 5, 10)
    adds <- MetaboCoreUtils::adductNames()
    res <- .mass_to_mz_df(mass, adducts = adds)
    expect_true(is.data.frame(res))
    expect_equal(colnames(res), c("index", "adduct", "mz"))
    expect_equal(res$index, rep(1:4, length(adds)))
    expect_equal(res$adduct, rep(adds, each = 4))
})
