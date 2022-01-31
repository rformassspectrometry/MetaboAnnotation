test_that("Mass2MzParam works", {
  res <- Mass2MzParam()
  expect_true(is(res, "Mass2MzParam"))

  expect_error(Mass2MzParam(tolerance = 1:3), "positive number")
  expect_error(Mass2MzParam(ppm = -4), "positive number")
  expect_error(Mass2MzParam(adducts = c("[M+H]+", "adduct2")), "Unknown")

  adds <- data.frame(mass_add = c(1, 2, 3), mass_multi = c(1, 2, 0.5))
  rownames(adds) <- c("a", "b", "c")
  res <- Mass2MzParam(adducts = adds)
  expect_true(is(res, "Mass2MzParam"))
})

test_that("Mass2MzRtParam works", {
  res <- Mass2MzRtParam()
  expect_true(is(res, "Mass2MzRtParam"))

  expect_error(Mass2MzRtParam(tolerance = 1:3), "positive number")
  expect_error(Mass2MzRtParam(toleranceRt = -1), "positive number")
  expect_error(Mass2MzRtParam(ppm = -4), "positive number")
  expect_error(Mass2MzRtParam(adducts = c("[M+H]+", "adduct2")), "Unknown")
})

test_that("MzParam works", {
  res <- MzParam()
  expect_true(is(res, "MzParam"))

  expect_error(MzParam(tolerance = 1:3), "positive number")
  expect_error(MzParam(ppm = -4), "positive number")
})

test_that("MzRtParam works", {
  res <- MzRtParam()
  expect_true(is(res, "MzRtParam"))

  expect_error(MzRtParam(tolerance = 1:3), "positive number")
  expect_error(MzRtParam(toleranceRt = -1), "positive number")
  expect_error(MzRtParam(ppm = -4), "positive number")
})

test_that("matchMz,Mass2MzParam works", {

  cmpds <- data.frame(
    name = c("Tryptophan", "Leucine", "Isoleucine"),
    formula = c("C11H12N2O2", "C6H13NO2", "C6H13NO2"),
    exactmass = c(204.089878, 131.094629, 131.094629)
  )

  adducts <- c("[M+H]+", "[M+Na]+")

  x <- data.frame(
    mz = c(mass2mz(204.089878, "[M+H]+"),
           mass2mz(131.094629, "[M+H]+"),
           mass2mz(204.089878, "[M+Na]+") + 1e-6)
  )

  par <- Mass2MzParam(adducts = adducts, tolerance = 0, ppm = 20)
  res <- matchMz(x, cmpds, par)
  expect_equal(query(res), x)
  expect_equal(target(res), cmpds)
  expect_equal(res@matches$query_idx, c(1, 2, 2, 3))
  expect_equal(res@matches$target_idx, c(1, 2, 3, 1))
  expect_equal(res@matches$score, c(0, 0, 0, 1e-6))
  expect_equal(res@matches$ppm, c(0, 0, 0, 1 / mass2mz(204.089878, "[M+Na]+")))

  par <- Mass2MzParam(adducts = adducts, tolerance = 0, ppm = 0)
  res <- matchMz(x, cmpds, par)
  expect_equal(query(res), x)
  expect_equal(target(res), cmpds)
  expect_equal(res@matches$query_idx, c(1, 2, 2))
  expect_equal(res@matches$target_idx, c(1, 2, 3))
  expect_equal(res@matches$score, c(0, 0, 0))
  expect_equal(res@matches$ppm, c(0, 0, 0))

  ## no matches
  adducts <- c("[M+Li]+", "[M+K]+")
  par <- Mass2MzParam(adducts = adducts, tolerance = 0, ppm = 20)
  res <- matchMz(x, cmpds, par)
  expect_true(is(res, "Matched"))
  expect_equal(query(res), x)
  expect_equal(target(res), cmpds)
  expect_true(nrow(res@matches) == 0)

  ## with a custom adduct definition
  adducts <- data.frame(mass_add = c(1, 3), mass_multi = c(1, 0.5))
  rownames(adducts) <- c("a", "b")
  x <- data.frame(
    mz = c(mass2mz(204.089878, adducts["a", ]),
           mass2mz(131.094629, adducts["a", ]),
           mass2mz(204.089878, adducts["b", ]) - 2e-6)
  )

  par <- Mass2MzParam(adducts = adducts, tolerance = 0, ppm = 20)
  res <- matchMz(x, cmpds, par)
  expect_equal(query(res), x)
  expect_equal(target(res), cmpds)
  expect_equal(res@matches$query_idx, c(1, 2, 2, 3))
  expect_equal(res@matches$target_idx, c(1, 2, 3, 1))
  expect_equal(res@matches$score, c(0, 0, 0, 2e-6))
  expect_equal(res@matches$ppm,
               c(0, 0, 0, -2 / mass2mz(204.089878, adducts["b", ])))

  par <- Mass2MzParam(adducts = adducts, tolerance = 0, ppm = 0)
  res <- matchMz(x, cmpds, par)
  expect_equal(query(res), x)
  expect_equal(target(res), cmpds)
  expect_equal(res@matches$query_idx, c(1, 2, 2))
  expect_equal(res@matches$target_idx, c(1, 2, 3))
  expect_equal(res@matches$score, c(0, 0, 0))
  expect_equal(res@matches$ppm, c(0, 0, 0))

  expect_error(matchMz(x, cmpds, par, massColname = "other"), "other")

  ## numeric, data.frame
  res <- matchMz(x$mz, cmpds, par)
  expect_equal(query(res), x$mz)
  expect_equal(target(res), cmpds)
  expect_equal(res@matches$query_idx, c(1, 2, 2))
  expect_equal(res@matches$target_idx, c(1, 2, 3))
  expect_equal(res@matches$score, c(0, 0, 0))
  expect_equal(res@matches$ppm, c(0, 0, 0))
  

  expect_error(matchMz(x$mz, cmpds, par, massColname = "other"), "other")

  ## data.frame, numeric
  res <- matchMz(x, cmpds$exactmass, par)
  expect_equal(query(res), x)
  expect_equal(target(res), cmpds$exactmass)
  expect_equal(res@matches$query_idx, c(1, 2, 2))
  expect_equal(res@matches$target_idx, c(1, 2, 3))
  expect_equal(res@matches$score, c(0, 0, 0))
  expect_equal(res@matches$ppm, c(0, 0, 0))
  

  expect_error(matchMz(x, cmpds$exactmass, par, mzColname = "other"), "other")

  adducts <- data.frame(mass_add = c(2, 4), mass_multi = c(2, 1))
  par <- Mass2MzParam(adducts = adducts, tolerance = 0, ppm = 20)
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

test_that("matchMz,Mass2MzRtParam works", {

  cmpds <- data.frame(
    name = c("Tryptophan", "Leucine", "Isoleucine"),
    formula = c("C11H12N2O2", "C6H13NO2", "C6H13NO2"),
    exactmass = c(204.089878, 131.094629, 131.094629),
    rt = c(150, 140, 140)
  )

  adducts <- c("[M+H]+", "[M+Na]+")

  x <- data.frame(
    mz = c(mass2mz(204.089878, "[M+H]+"),
           mass2mz(131.094629, "[M+H]+"),
           mass2mz(204.089878, "[M+Na]+") + 1e-6),
    rt = c(150, 140, 150.1)
  )

  par <- Mass2MzRtParam(adducts = adducts, tolerance = 0, ppm = 20,
                              toleranceRt = 0)
  res <- matchMz(x, cmpds, par)
  expect_equal(query(res), x)
  expect_equal(target(res), cmpds)
  expect_equal(res@matches$query_idx, c(1, 2, 2))
  expect_equal(res@matches$target_idx, c(1, 2, 3))
  expect_equal(res@matches$score, c(0, 0, 0))
  expect_equal(res@matches$ppm, c(0, 0, 0))
  expect_equal(res@matches$score_rt, c(0, 0, 0))

  par <- Mass2MzRtParam(adducts = adducts, tolerance = 0, ppm = 20,
                            toleranceRt = 0.2)
  res <- matchMz(x, cmpds, par)
  expect_equal(query(res), x)
  expect_equal(target(res), cmpds)
  expect_equal(res@matches$query_idx, c(1, 2, 2, 3))
  expect_equal(res@matches$target_idx, c(1, 2, 3, 1))
  expect_equal(res@matches$score, c(0, 0, 0, 1e-6))
  expect_equal(res@matches$ppm, c(0, 0, 0, 1 / mass2mz(204.089878, "[M+Na]+")))
  expect_equal(res@matches$score_rt, c(0, 0, 0, 0.1))

  par <- Mass2MzRtParam(adducts = adducts, tolerance = 0, ppm = 0,
                              toleranceRt = 0.2)
  res <- matchMz(x, cmpds, par)
  expect_equal(query(res), x)
  expect_equal(target(res), cmpds)
  expect_equal(res@matches$query_idx, c(1, 2, 2))
  expect_equal(res@matches$target_idx, c(1, 2, 3))
  expect_equal(res@matches$score, c(0, 0, 0))
  expect_equal(res@matches$ppm, c(0, 0, 0))
  expect_equal(res@matches$score_rt, c(0, 0, 0))


  ## no matches
  adducts <- c("[M+Li]+", "[M+K]+")
  par <- Mass2MzRtParam(adducts = adducts, tolerance = 0, ppm = 20,
                            toleranceRt = 0.2)
  res <- matchMz(x, cmpds, par)
  expect_true(is(res, "Matched"))
  expect_equal(query(res), x)
  expect_equal(target(res), cmpds)
  expect_true(nrow(res@matches) == 0)
})

test_that(".getMatchesMzRt works", {
  trgt <- data.frame(index = rep(1:7, 2),
                     adduct = c(rep("A", 7), rep("B", 7)),
                     mz = c(11:17, 12:18),
                     rt = rep(21:27, 2))
  trgt <- trgt[order(trgt$mz), ]

  res <- .getMatchesMzRt(queryIndex = 3, queryMz = 13, queryRt = 23,
                         target = trgt, tolerance = 0, ppm = 0, toleranceRt = 0)

  expect_true(is.data.frame(res))
  expect_equal(res$query_idx, c(3))
  expect_equal(res$target_idx, c(3))

  res <- .getMatchesMzRt(queryIndex = 3, queryMz = 13, queryRt = 24,
                         target = trgt, tolerance = 0, ppm = 0, toleranceRt = 0)
  expect_true(is.data.frame(res))
  expect_true(nrow(res) == 0)

  res <- .getMatchesMzRt(queryIndex = 3, queryMz = 13, queryRt = 24,
                         target = trgt, tolerance = 0, ppm = 0, toleranceRt = 1)
  expect_true(is.data.frame(res))
  expect_equal(res$query_idx, c(3))
  expect_equal(res$target_idx, c(3))

  res <- .getMatchesMzRt(queryIndex = 3, queryMz = 13, queryRt = 23,
                         target = trgt, tolerance = 1, ppm = 0, toleranceRt = 1)
  expect_true(is.data.frame(res))
  expect_equal(res$query_idx, c(3, 3, 3, 3, 3))
  expect_equal(res$target_idx, c(2, 3, 2, 4, 3))
  expect_equal(res$adduct, c("A", "A", "B", "A", "B"))
  expect_equal(res$score, c(1, 0, 0, 1, 1))
  expect_equal(res$score_rt, c(1, 0, 1, 1, 0))
})


test_that(".mass_to_mz_df works", {
    mass <- c(4, 2, 5, 10)
    adds <- MetaboCoreUtils::adductNames()
    res <- .mass_to_mz_df(mass, adducts = adds)
    expect_true(is.data.frame(res))
    expect_equal(colnames(res), c("index", "adduct", "mz"))
    expect_equal(res$index, rep(1:4, length(adds)))
    expect_equal(res$adduct, rep(adds, each = 4))

    adds <- MetaboCoreUtils::adducts()
    res1 <- .mass_to_mz_df(mass, adducts = adds)
    expect_true(is.data.frame(res))
    expect_equal(colnames(res), c("index", "adduct", "mz"))
    expect_equal(res$index, rep(1:4, nrow(adds)))
    expect_equal(res$adduct, rep(adds$name, each = 4))
})

test_that("matchMz, MzParam works", {
  qry <- data.frame(mz = c(150, 170, 179))
  trgt <- data.frame(mz = seq(110, 200, 10))
  # qry <- data.frame(mz = c(150, 170, 179))
  # trgt <- seq(110, 200, 10)
  # qry <- c(150, 170, 179)
  # trgt <- data.frame(mz = seq(110, 200, 10))
  # qry <- c(150, 170, 179)
  # trgt <- seq(110, 200, 10)

  par <- MzParam(tolerance = 0)
  res <- matchMz(qry, trgt, par)
  expect_equal(query(res), qry)
  expect_equal(target(res), trgt)
  expect_equal(res@matches$query_idx, c(1, 2))
  expect_equal(res@matches$target_idx, c(5, 7))
  expect_equal(res@matches$score, c(0, 0))
  expect_equal(res@matches$ppm, c(0, 0))
  
  ## no matches
  res <- matchMz(qry + 0.1, trgt, par)
  expect_true(is(res, "Matched"))
  expect_equal(query(res), qry + 0.1)
  expect_equal(target(res), trgt)
  expect_true(nrow(res@matches) == 0)

  # positive tolerance
  par <- MzParam(tolerance = 10)
  res <- matchMz(qry, trgt, par)
  expect_equal(query(res), qry)
  expect_equal(target(res), trgt)
  expect_equal(res@matches$query_idx, c(1, 1, 1, 2, 2, 2, 3, 3))
  expect_equal(res@matches$target_idx, c(4, 5, 6, 6, 7, 8, 7, 8))
  expect_equal(res@matches$score, c(10, 0, 10, 10, 0, 10, 9, 1))
  expect_equal(res@matches$ppm, c(10, 0, -10, 10, 0, -10, 9, -1) /
                 c(140, 150, 160, 160, 170, 180, 170, 180) * 10^6)
})

test_that("matchMz, MzRtParam works", {

  qry <- data.frame(mz = c(13, 14.1, 17, 18), rt = c(23, 24, 26.8, 23))
  trgt <- data.frame(mz = 11:17, rt = 21:27)

  par <- MzRtParam(tolerance = 0, ppm = 0, toleranceRt = 0)
  res <- matchMz(qry, trgt, par)
  expect_equal(query(res), qry)
  expect_equal(target(res), trgt)
  expect_equal(res@matches$query_idx, c(1))
  expect_equal(res@matches$target_idx, c(3))
  expect_equal(res@matches$score, c(0))
  expect_equal(res@matches$ppm, c(0))
  expect_equal(res@matches$score_rt, c(0))

  par <- MzRtParam(tolerance = 0.1, ppm = 0, toleranceRt = 0)
  res <- matchMz(qry, trgt, par)
  expect_equal(res@matches$query_idx, c(1, 2))
  expect_equal(res@matches$target_idx, c(3, 4))
  expect_equal(res@matches$score, c(0, 0.1))
  expect_equal(res@matches$ppm, c(0, 0.1 / 14 * 10^6))
  expect_equal(res@matches$score_rt, c(0, 0))

  par <- MzRtParam(tolerance = 0.1, ppm = 0, toleranceRt = 0.2)
  res <- matchMz(qry, trgt, par)
  expect_equal(res@matches$query_idx, c(1, 2, 3))
  expect_equal(res@matches$target_idx, c(3, 4, 7))
  expect_equal(res@matches$score, c(0, 0.1, 0))
  expect_equal(res@matches$ppm, c(0, 0.1 / 14 * 10^6, 0))
  expect_equal(res@matches$score_rt, c(0, 0, 0.2))

  ## no matches
  res <- matchMz(qry + 0.5, trgt, par)
  expect_true(is(res, "Matched"))
  expect_equal(query(res), qry + 0.5)
  expect_equal(target(res), trgt)
  expect_true(nrow(res@matches) == 0)
})

test_that("matchMz, MzRtParam works", {
  library(SummarizedExperiment)
  qry <- SummarizedExperiment(
    assays = data.frame(matrix(NA, 4, 2)),
    colData = data.frame(cD1 = c(NA, NA), cD2 = c(NA, NA)),
    rowData = data.frame(mz = c(13, 14.1, 17, 18), rt = c(23, 24, 26.8, 23)))
  trgt <- data.frame(mz = 11:17, rt = 21:27)

  par <- MzRtParam(tolerance = 0, ppm = 0, toleranceRt = 0)
  res <- matchMz(qry, trgt, par)
  expect_equal(query(res), qry)
  expect_equal(target(res), trgt)
  expect_equal(res@matches$query_idx, c(1))
  expect_equal(res@matches$target_idx, c(3))
  expect_equal(res@matches$score, c(0))
  expect_equal(res@matches$ppm, c(0))
  expect_equal(res@matches$score_rt, c(0))

  par <- MzRtParam(tolerance = 0.1, ppm = 0, toleranceRt = 0)
  res <- matchMz(qry, trgt, par)
  expect_equal(res@matches$query_idx, c(1, 2))
  expect_equal(res@matches$target_idx, c(3, 4))
  expect_equal(res@matches$score, c(0, 0.1))
  expect_equal(res@matches$ppm, c(0, 0.1 / 14 * 10^6))
  expect_equal(res@matches$score_rt, c(0, 0))

  par <- MzRtParam(tolerance = 0.1, ppm = 0, toleranceRt = 0.2)
  res <- matchMz(qry, trgt, par)
  expect_equal(res@matches$query_idx, c(1, 2, 3))
  expect_equal(res@matches$target_idx, c(3, 4, 7))
  expect_equal(res@matches$score, c(0, 0.1, 0))
  expect_equal(res@matches$ppm, c(0, 0.1 / 14 * 10^6, 0))
  expect_equal(res@matches$score_rt, c(0, 0, 0.2))

  ## no matches
  res <- matchMz(qry , trgt + 0.5, par)
  expect_true(is(res, "Matched"))
  expect_equal(query(res), qry)
  expect_equal(target(res), trgt + 0.5)
  expect_true(nrow(res@matches) == 0)
})
