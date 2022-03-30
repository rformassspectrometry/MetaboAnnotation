test_that("ValueParam works", {
    res <- ValueParam()
    expect_true(is(res, "ValueParam"))

    expect_error(ValueParam(tolerance = 1:3), "positive number")
    expect_error(ValueParam(ppm = -4), "positive number")
})

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

test_that("Mz2MassParam works", {
    res <- Mz2MassParam()
    expect_true(is(res, "Mz2MassParam"))

    expect_error(Mz2MassParam(tolerance = 1:3), "positive number")
    expect_error(Mz2MassParam(ppm = -4), "positive number")

    adds <- data.frame(mass_add = c(1, 2, 3), mass_multi = c(1, 2, 0.5))
    rownames(adds) <- c("a", "b", "c")
    res <- Mz2MassParam(queryAdducts = adds)
    expect_true(is(res, "Mz2MassParam"))
})

test_that("Mz2MassRtParam works", {
    res <- Mz2MassRtParam()
    expect_true(is(res, "Mz2MassRtParam"))

    expect_error(Mz2MassRtParam(tolerance = 1:3), "positive number")
    expect_error(Mz2MassRtParam(ppm = -4), "positive number")
    expect_error(Mz2MassRtParam(toleranceRt = -1), "positive number")

    adds <- data.frame(mass_add = c(1, 2, 3), mass_multi = c(1, 2, 0.5))
    rownames(adds) <- c("a", "b", "c")
    res <- Mz2MassRtParam(queryAdducts = adds)
    expect_true(is(res, "Mz2MassRtParam"))
})

test_that("matchValues, ValueParam works", {
    qry <- data.frame(value = c(150, 170, 179))
    trgt <- data.frame(value = seq(110, 200, 10))
    par <- ValueParam(tolerance = 0)

    ## numeric, data.frame
    expect_error(matchValues(qry[, "value"], trgt, par), "provided")
    expect_error(matchValues(qry[, "value"], trgt, par, valueColname = "some"),
                 "Missing")
    res <- matchValues(qry[, "value"], trgt, par, valueColname = "value")
    expect_s4_class(res, "Matched")
    expect_equal(res@metadata$param, par)
    expect_equal(res@query, qry[, "value"])
    expect_equal(res@target, trgt)

    ## data.frame, numeric
    expect_error(matchValues(qry, trgt[, "value"], par), "provided")
    expect_error(matchValues(qry, trgt[, "value"], par, valueColname = "some"),
                 "Missing")
    res <- matchValues(qry, trgt[, "value"], par, valueColname = "value")
    expect_s4_class(res, "Matched")
    expect_equal(res@metadata$param, par)
    expect_equal(res@query, qry)
    expect_equal(res@target, trgt[, "value"])

    ## data.frame, data.frame
    expect_error(matchValues(qry, trgt, par), "`valueColname` has to be provided.")
    expect_error(matchValues(qry, trgt, par,
                         valueColname = c("notcolumn", "value")), "Missing")
    expect_error(matchValues(qry, trgt, par,
                         valueColname = c("value", "notcolumn")), "Missing")
    res <- matchValues(qry, trgt, par, valueColname = "value")
    expect_equal(query(res), qry)
    expect_equal(target(res), trgt)
    expect_equal(res@matches$query_idx, c(1, 2))
    expect_equal(res@matches$target_idx, c(5, 7))
    expect_equal(res@matches$score, c(0, 0))
    expect_equal(res@matches$ppm_error, c(0, 0))

    ## no matches
    res <- matchValues(qry + 0.1, trgt, par, valueColname = "value")
    expect_true(is(res, "Matched"))
    expect_equal(query(res), qry + 0.1)
    expect_equal(target(res), trgt)
    expect_true(nrow(res@matches) == 0)

    # positive tolerance
    par <- ValueParam(tolerance = 10)
    res <- matchValues(qry, trgt, par, valueColname = "value")
    expect_equal(query(res), qry)
    expect_equal(target(res), trgt)
    expect_equal(res@matches$query_idx, c(1, 1, 1, 2, 2, 2, 3, 3))
    expect_equal(res@matches$target_idx, c(4, 5, 6, 6, 7, 8, 7, 8))
    expect_equal(res@matches$score, c(10, 0, -10, 10, 0, -10, 9, -1))
    expect_equal(res@matches$ppm_error, c(10, 0, 10, 10, 0, 10, 9, 1) /
                     c(140, 150, 160, 160, 170, 180, 170, 180) * 10^6)
})

test_that("matchValues,Mass2MzParam works", {

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

    ## numeric, data.frame

    ## data.frame, data.frame
    expect_error(matchValues(x, cmpds, par, mzColname = "mzz"), "Missing")
    res <- matchValues(x, cmpds, par)
    expect_equal(query(res), x)
    expect_equal(target(res), cmpds)
    expect_equal(res@matches$query_idx, c(1, 2, 2, 3))
    expect_equal(res@matches$target_idx, c(1, 2, 3, 1))
    expect_equal(res@matches$score, c(0, 0, 0, 1e-6))
    expect_equal(res@matches$ppm_error,
                 c(0, 0, 0, 1 / mass2mz(204.089878, "[M+Na]+")))

    par <- Mass2MzParam(adducts = adducts, tolerance = 0, ppm = 0)
    res <- matchValues(x, cmpds, par)
    expect_equal(query(res), x)
    expect_equal(target(res), cmpds)
    expect_equal(res@matches$query_idx, c(1, 2, 2))
    expect_equal(res@matches$target_idx, c(1, 2, 3))
    expect_equal(res@matches$score, c(0, 0, 0))
    expect_equal(res@matches$ppm_error, c(0, 0, 0))

    ## no matches
    adducts <- c("[M+Li]+", "[M+K]+")
    par <- Mass2MzParam(adducts = adducts, tolerance = 0, ppm = 20)
    res <- matchValues(x, cmpds, par)
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
    res <- matchValues(x, cmpds, par)
    expect_equal(query(res), x)
    expect_equal(target(res), cmpds)
    expect_equal(res@matches$query_idx, c(1, 2, 2, 3))
    expect_equal(res@matches$target_idx, c(1, 2, 3, 1))
    expect_equal(res@matches$score, c(0, 0, 0, -2e-6))
    expect_equal(res@matches$ppm_error,
                 c(0, 0, 0, 2 / mass2mz(204.089878, adducts["b", ])))

    par <- Mass2MzParam(adducts = adducts, tolerance = 0, ppm = 0)
    res <- matchValues(x, cmpds, par)
    expect_equal(query(res), x)
    expect_equal(target(res), cmpds)
    expect_equal(res@matches$query_idx, c(1, 2, 2))
    expect_equal(res@matches$target_idx, c(1, 2, 3))
    expect_equal(res@matches$score, c(0, 0, 0))
    expect_equal(res@matches$ppm_error, c(0, 0, 0))

    expect_error(matchValues(x, cmpds, par, massColname = "other"), "other")

    ## numeric, data.frame
    res <- matchValues(x$mz, cmpds, par)
    expect_equal(query(res), x$mz)
    expect_equal(target(res), cmpds)
    expect_equal(res@matches$query_idx, c(1, 2, 2))
    expect_equal(res@matches$target_idx, c(1, 2, 3))
    expect_equal(res@matches$score, c(0, 0, 0))
    expect_equal(res@matches$ppm_error, c(0, 0, 0))

    expect_error(matchValues(x$mz, cmpds, par, massColname = "other"), "other")

    ## data.frame, numeric
    res <- matchValues(x, cmpds$exactmass, par)
    expect_equal(query(res), x)
    expect_equal(target(res), cmpds$exactmass)
    expect_equal(res@matches$query_idx, c(1, 2, 2))
    expect_equal(res@matches$target_idx, c(1, 2, 3))
    expect_equal(res@matches$score, c(0, 0, 0))
    expect_equal(res@matches$ppm_error, c(0, 0, 0))

    expect_error(matchValues(x, cmpds$exactmass, par, mzColname = "other"), "other")

    adducts <- data.frame(mass_add = c(2, 4), mass_multi = c(2, 1))
    par <- Mass2MzParam(adducts = adducts, tolerance = 0, ppm = 20)
    res <- matchValues(x, cmpds, par)
    expect_true(is(res, "Matched"))
    expect_equal(query(res), x)
    expect_equal(target(res), cmpds)
    expect_true(nrow(res@matches) == 0)
})

test_that(".getMatches works", {
    trgt <- data.frame(index = c(1, 2, 3, 1, 2, 3),
                       mz = c(1, 1.2, 1.21, 1.2, 2, 2.1),
                       adduct = c("A", "A", "A", "B", "B", "B"))
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

test_that("matchValues,Mass2MzRtParam works", {

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

    ## data.frame, data.frame
    res <- matchValues(x, cmpds, par, rtColname = "rt")
    expect_equal(query(res), x)
    expect_error(matchValues(x, cmpds, par, rtColname = c("r", "rt")), "Missin")
    expect_error(matchValues(x, cmpds, par, rtColname = c("rt", "r")), "Missin")
    expect_error(matchValues(x, cmpds, par, mzColname = "r"), "Missin")
    expect_error(matchValues(x, cmpds, par, massColname = "r"), "Missin")

    res <- matchValues(x, cmpds, par)
    expect_equal(query(res), x)
    expect_equal(target(res), cmpds)
    expect_equal(res@matches$query_idx, c(1, 2, 2))
    expect_equal(res@matches$target_idx, c(1, 2, 3))
    expect_equal(res@matches$score, c(0, 0, 0))
    expect_equal(res@matches$ppm_error, c(0, 0, 0))
    expect_equal(res@matches$score_rt, c(0, 0, 0))

    par <- Mass2MzRtParam(adducts = adducts, tolerance = 0, ppm = 20,
                          toleranceRt = 0.2)
    res <- matchValues(x, cmpds, par)
    expect_equal(query(res), x)
    expect_equal(target(res), cmpds)
    expect_equal(res@matches$query_idx, c(1, 2, 2, 3))
    expect_equal(res@matches$target_idx, c(1, 2, 3, 1))
    expect_equal(res@matches$score, c(0, 0, 0, 1e-6))
    expect_equal(res@matches$ppm_error,
                 c(0, 0, 0, 1 / mass2mz(204.089878, "[M+Na]+")))
    expect_equal(res@matches$score_rt, c(0, 0, 0, 0.1))

    par <- Mass2MzRtParam(adducts = adducts, tolerance = 0, ppm = 0,
                          toleranceRt = 0.2)
    res <- matchValues(x, cmpds, par)
    expect_equal(query(res), x)
    expect_equal(target(res), cmpds)
    expect_equal(res@matches$query_idx, c(1, 2, 2))
    expect_equal(res@matches$target_idx, c(1, 2, 3))
    expect_equal(res@matches$score, c(0, 0, 0))
    expect_equal(res@matches$ppm_error, c(0, 0, 0))
    expect_equal(res@matches$score_rt, c(0, 0, 0))

    ## no matches
    adducts <- c("[M+Li]+", "[M+K]+")
    par <- Mass2MzRtParam(adducts = adducts, tolerance = 0, ppm = 20,
                          toleranceRt = 0.2)
    res <- matchValues(x, cmpds, par)
    expect_true(is(res, "Matched"))
    expect_equal(query(res), x)
    expect_equal(target(res), cmpds)
    expect_true(nrow(res@matches) == 0)
})

test_that(".getMatchesMzRt works", {
    trgt <- data.frame(index = rep(1:7, 2),
                       mz = c(11:17, 12:18),
                       adduct = c(rep("A", 7), rep("B", 7)),
                       rt = rep(21:27, 2))
    trgt <- trgt[order(trgt$mz), ]

    res <- .getMatchesMzRt(queryIndex = 3, queryMz = 13, queryRt = 23,
                           target = trgt, tolerance = 0, ppm = 0,
                           toleranceRt = 0)

    expect_true(is.data.frame(res))
    expect_equal(res$query_idx, c(3))
    expect_equal(res$target_idx, c(3))

    res <- .getMatchesMzRt(queryIndex = 3, queryMz = 13, queryRt = 24,
                           target = trgt, tolerance = 0, ppm = 0,
                           toleranceRt = 0)
    expect_true(is.data.frame(res))
    expect_true(nrow(res) == 0)

    res <- .getMatchesMzRt(queryIndex = 3, queryMz = 13, queryRt = 24,
                           target = trgt, tolerance = 0, ppm = 0,
                           toleranceRt = 1)
    expect_true(is.data.frame(res))
    expect_equal(res$query_idx, c(3))
    expect_equal(res$target_idx, c(3))

    res <- .getMatchesMzRt(queryIndex = 3, queryMz = 13, queryRt = 23,
                           target = trgt, tolerance = 1, ppm = 0,
                           toleranceRt = 1)
    expect_true(is.data.frame(res))
    expect_equal(res$query_idx, c(3, 3, 3, 3, 3))
    expect_equal(res$target_idx, c(2, 3, 2, 4, 3))
    expect_equal(res$adduct, c("A", "A", "B", "A", "B"))
    expect_equal(res$score, c(1, 0, 0, -1, -1))
    expect_equal(res$ppm_error, c(1, 0, 0, 1, 1) / c(12, 13, 13, 14, 14) * 10^6)
    expect_equal(res$score_rt, c(1, 0, 1, -1, 0))
})

test_that(".mass_to_mz_df works", {
    mass <- c(4, 2, 5, 10)
    adds <- MetaboCoreUtils::adductNames()
    res <- .mass_to_mz_df(mass, adducts = adds)
    expect_true(is.data.frame(res))
    expect_equal(colnames(res), c("index", "mz", "adduct"))
    expect_equal(res$index, rep(1:4, length(adds)))
    expect_equal(res$adduct, rep(adds, each = 4))

    adds <- MetaboCoreUtils::adducts()
    res1 <- .mass_to_mz_df(mass, adducts = adds)
    expect_true(is.data.frame(res))
    expect_equal(colnames(res), c("index", "mz", "adduct"))
    expect_equal(res$index, rep(1:4, nrow(adds)))
    expect_equal(res$adduct, rep(adds$name, each = 4))
})

test_that("matchValues, MzParam works", {
    qry <- data.frame(mz = c(150, 170, 179))
    trgt <- data.frame(mz = seq(110, 200, 10))
    # qry <- data.frame(mz = c(150, 170, 179))
    # trgt <- seq(110, 200, 10)
    # qry <- c(150, 170, 179)
    # trgt <- data.frame(mz = seq(110, 200, 10))
    # qry <- c(150, 170, 179)
    # trgt <- seq(110, 200, 10)
    par <- MzParam(tolerance = 0)

    ## numeric, data.frame
    expect_error(matchValues(qry$mz, trgt, par, mzColname = "mzz"), "Missing")
    res <- matchValues(qry$mz, trgt, par)
    expect_s4_class(res, "Matched")
    expect_equal(query(res), qry$mz)
    expect_equal(target(res), trgt)
    expect_equal(res@metadata$param, par)

    ## data.frame, numeric
    expect_error(matchValues(qry, trgt$mz, par, mzColname = "mzz"), "Missing")
    res <- matchValues(qry, trgt$mz, par)
    expect_s4_class(res, "Matched")
    expect_equal(query(res), qry)
    expect_equal(target(res), trgt$mz)
    expect_equal(res@metadata$param, par)

    ## data.frame, data.frame
    res <- matchValues(qry, trgt, par, mzColname = "mz")
    expect_s4_class(res, "Matched")
    expect_error(matchValues(qry, trgt, par, mzColname = c("mzz", "mz")),
                 "Missing")
    expect_error(matchValues(qry, trgt, par, mzColname = c("mz", "mzz")),
                 "Missing")

    res <- matchValues(qry, trgt, par)
    expect_equal(query(res), qry)
    expect_equal(target(res), trgt)
    expect_equal(res@matches$query_idx, c(1, 2))
    expect_equal(res@matches$target_idx, c(5, 7))
    expect_equal(res@matches$score, c(0, 0))
    expect_equal(res@matches$ppm_error, c(0, 0))

    ## no matches
    res <- matchValues(qry + 0.1, trgt, par)
    expect_true(is(res, "Matched"))
    expect_equal(query(res), qry + 0.1)
    expect_equal(target(res), trgt)
    expect_true(nrow(res@matches) == 0)

    # positive tolerance
    par <- MzParam(tolerance = 10)
    res <- matchValues(qry, trgt, par)
    expect_equal(query(res), qry)
    expect_equal(target(res), trgt)
    expect_equal(res@matches$query_idx, c(1, 1, 1, 2, 2, 2, 3, 3))
    expect_equal(res@matches$target_idx, c(4, 5, 6, 6, 7, 8, 7, 8))
    expect_equal(res@matches$score, c(10, 0, -10, 10, 0, -10, 9, -1))
    expect_equal(res@matches$ppm_error, c(10, 0, 10, 10, 0, 10, 9, 1) /
                     c(140, 150, 160, 160, 170, 180, 170, 180) * 10^6)
})

test_that("matchValues, MzRtParam works", {

    qry <- data.frame(mz = c(13, 14.1, 17, 18), rt = c(23, 24, 26.8, 23))
    trgt <- data.frame(mz = 11:17, rt = 21:27)

    par <- MzRtParam(tolerance = 0, ppm = 0, toleranceRt = 0)

    ## data.frame, data.frame
    res <- matchValues(qry, trgt, par, rtColname = "rt")
    expect_equal(query(res), qry)
    expect_error(matchValues(qry, trgt, par, mzColname = c("r", "mz")), "Mis")
    expect_error(matchValues(qry, trgt, par, mzColname = c("mz", "m")), "Mis")
    expect_error(matchValues(qry, trgt, par, rtColname = c("r", "rt")), "Mis")
    expect_error(matchValues(qry, trgt, par, rtColname = c("rt", "r")), "Mis")

    res <- matchValues(qry, trgt, par)
    expect_equal(query(res), qry)
    expect_equal(target(res), trgt)
    expect_equal(res@matches$query_idx, c(1))
    expect_equal(res@matches$target_idx, c(3))
    expect_equal(res@matches$score, c(0))
    expect_equal(res@matches$ppm_error, c(0))
    expect_equal(res@matches$score_rt, c(0))

    par <- MzRtParam(tolerance = 0.1, ppm = 0, toleranceRt = 0)
    res <- matchValues(qry, trgt, par)
    expect_equal(res@matches$query_idx, c(1, 2))
    expect_equal(res@matches$target_idx, c(3, 4))
    expect_equal(res@matches$score, c(0, 0.1))
    expect_equal(res@matches$ppm_error, c(0, 0.1 / 14 * 10^6))
    expect_equal(res@matches$score_rt, c(0, 0))

    par <- MzRtParam(tolerance = 0.1, ppm = 0, toleranceRt = 0.2)
    res <- matchValues(qry, trgt, par)
    expect_equal(res@matches$query_idx, c(1, 2, 3))
    expect_equal(res@matches$target_idx, c(3, 4, 7))
    expect_equal(res@matches$score, c(0, 0.1, 0))
    expect_equal(res@matches$ppm_error, c(0, 0.1 / 14 * 10^6, 0))
    expect_equal(res@matches$score_rt, c(0, 0, -0.2))

    ## no matches
    res <- matchValues(qry + 0.5, trgt, par)
    expect_true(is(res, "Matched"))
    expect_equal(query(res), qry + 0.5)
    expect_equal(target(res), trgt)
    expect_true(nrow(res@matches) == 0)

    #### query `SummarizedExperiment`
    library(SummarizedExperiment)
    qry <- SummarizedExperiment(
        assays = data.frame(matrix(NA, 4, 2)),
        colData = data.frame(cD1 = c(NA, NA), cD2 = c(NA, NA)),
        rowData = data.frame(mz = c(13, 14.1, 17, 18),
                             rt = c(23, 24, 26.8, 23)))
    trgt <- data.frame(mz = 11:17, rt = 21:27)

    par <- MzRtParam(tolerance = 0, ppm = 0, toleranceRt = 0)
    res <- matchValues(qry, trgt, par)
    expect_equal(query(res), qry)
    expect_equal(target(res), trgt)
    expect_equal(res@matches$query_idx, c(1))
    expect_equal(res@matches$target_idx, c(3))
    expect_equal(res@matches$score, c(0))
    expect_equal(res@matches$ppm_error, c(0))
    expect_equal(res@matches$score_rt, c(0))

    par <- MzRtParam(tolerance = 0.1, ppm = 0, toleranceRt = 0)
    res <- matchValues(qry, trgt, par)
    expect_equal(res@matches$query_idx, c(1, 2))
    expect_equal(res@matches$target_idx, c(3, 4))
    expect_equal(res@matches$score, c(0, 0.1))
    expect_equal(res@matches$ppm_error, c(0, 0.1 / 14 * 10^6))
    expect_equal(res@matches$score_rt, c(0, 0))

    par <- MzRtParam(tolerance = 0.1, ppm = 0, toleranceRt = 0.2)
    res <- matchValues(qry, trgt, par)
    expect_equal(res@matches$query_idx, c(1, 2, 3))
    expect_equal(res@matches$target_idx, c(3, 4, 7))
    expect_equal(res@matches$score, c(0, 0.1, 0))
    expect_equal(res@matches$ppm_error, c(0, 0.1 / 14 * 10^6, 0))
    expect_equal(res@matches$score_rt, c(0, 0, -0.2))

    ## no matches
    res <- matchValues(qry , trgt + 0.5, par)
    expect_true(is(res, "Matched"))
    expect_equal(query(res), qry)
    expect_equal(target(res), trgt + 0.5)
    expect_true(nrow(res@matches) == 0)
})

test_that("matchValues, Mz2MassParam works", {

    m <- c(200, 300)
    qry <- c(100, as.numeric(mass2mz(m, c("[M+H]+", "[M+K]+"))) + c(0, 0, 0, 5))
    trgt <- c(mass2mz(m, "[M-H]-"), 400, 500)

    par <- Mz2MassParam(queryAdducts = c("[M+H]+", "[M+K]+"),
                        targetAdducts = "[M-H]-")

    ## numeric, data.frame
    expect_error(
        matchMz(qry, data.frame(mz = trgt), par, mzColname = "m"), "Mis")
    res <- matchMz(qry, data.frame(mz = trgt), par)
    expect_s4_class(res, "Matched")
    expect_equal(query(res), qry)
    expect_equal(target(res), data.frame(mz = trgt))
    expect_equal(res@metadata$param, par)
    expect_true(is(matchedData(res), "DataFrame"))

    ## data.frame, numeric
    expect_error(
        matchMz(data.frame(mz = qry), trgt, par, mzColname = "m"), "Mis")
    res <- matchMz(data.frame(mz = qry), trgt, par)
    expect_s4_class(res, "Matched")
    expect_equal(query(res), data.frame(mz = qry))
    expect_equal(target(res), trgt)
    expect_equal(res@metadata$param, par)

    a <- data.frame(mz = qry)
    b <- data.frame(mz = trgt)
    ## data.frame, data.frame
    expect_error(matchMz(a, b, par, mzColname = c("m", "mz")), "Missing")
    expect_error(matchMz(a, b, par, mzColname = c("mz", "m")), "Missing")

    res <- matchValues(qry, trgt, par)
    expect_equal(query(res), qry)
    expect_equal(target(res), trgt)
    expect_equal(res@matches$query_idx, c(2, 3, 4))
    expect_equal(res@matches$target_idx, c(1, 2, 1))
    expect_equal(res@matches$query_adduct, c("[M+H]+", "[M+H]+", "[M+K]+"))
    expect_equal(res@matches$target_adduct, rep("[M-H]-", 3))
    expect_equal(res@matches$score, c(0, 0, 0))
    expect_equal(res@matches$ppm_error, c(0, 0, 0))

    ## no matches
    res <- matchValues(qry + 0.1, trgt, par)
    expect_true(is(res, "Matched"))
    expect_equal(query(res), qry + 0.1)
    expect_equal(target(res), trgt)
    expect_true(nrow(res@matches) == 0)

    # positive tolerance
    par <- Mz2MassParam(queryAdducts = c("[M+H]+", "[M+K]+"),
                        targetAdducts = "[M-H]-",
                        tolerance = 10)
    res <- matchValues(qry, trgt, par)
    expect_equal(query(res), qry)
    expect_equal(target(res), trgt)
    expect_equal(res@matches$query_idx, c(2, 3, 4, 5))
    expect_equal(res@matches$target_idx, c(1, 2, 1, 2))
    expect_equal(res@matches$target_adduct, rep("[M-H]-", 4))
    expect_equal(res@matches$query_adduct, c(rep("[M+H]+", 2),
                                             rep("[M+K]+", 2)))
    expect_equal(res@matches$score, c(0, 0, 0, 5))
    expect_equal(res@matches$ppm_error, c(0, 0, 0, 5 / 300 * 10^6))

    # `query`and `target` data.frames
    qry_df <- data.frame(mzq = qry, other = 1)
    trgt_df <- data.frame(other = "b", mzt = trgt)
    res <- matchValues(qry_df, trgt_df, par, mzColname = c("mzq", "mzt"))
    expect_equal(query(res), qry_df)
    expect_equal(target(res), trgt_df)
    expect_equal(res@matches$query_idx, c(2, 3, 4, 5))
    expect_equal(res@matches$target_adduct, rep("[M-H]-", 4))
    expect_equal(res@matches$query_adduct, c(rep("[M+H]+", 2),
                                             rep("[M+K]+", 2)))
    expect_equal(res@matches$score, c(0, 0, 0, 5))
    expect_equal(res@matches$ppm_error, c(0, 0, 0, 5 / 300 * 10^6))

    par <- Mz2MassParam(queryAdducts = "[M]+", targetAdducts = "[M]+")
    res <- matchValues(qry, qry, par)
    expect_equal(res@matches$query_idx, 1:5)
    expect_equal(res@matches$target_idx, 1:5)

    res_2 <- matchValues(qry, qry, ValueParam())
    expect_equal(res_2@matches$query_idx, 1:5)
    expect_equal(res_2@matches$target_idx, 1:5)
})

test_that("matchValues, Mz2MassRtParam works", {

    m <- c(200, 300)
    qry <- c(100, as.numeric(mass2mz(m, c("[M+H]+", "[M+K]+"))) + c(0, 0, 0, 5))
    trgt <- c(mass2mz(m, "[M-H]-"), 400, 500)
    qry_df <- data.frame(mz = qry, rt = c(50, 100, 151, 100, 150))
    trgt_df <- data.frame(mz = trgt, rt = c(100, 150, 200, 250))

    par <- Mz2MassRtParam(queryAdducts = c("[M+H]+", "[M+K]+"),
                          targetAdducts = "[M-H]-")
    res <- matchValues(qry_df, trgt_df, par, mzColname = "mz")
    expect_s4_class(res, "Matched")
    res <- matchValues(qry_df, trgt_df, par, rtColname = "rt")
    expect_s4_class(res, "Matched")

    expect_error(
        matchValues(qry_df, trgt_df, par, mzColname = c("m", "mz")), "Mis")
    expect_error(
        matchValues(qry_df, trgt_df, par, mzColname = c("mz", "m")), "Mis")
    expect_error(
        matchValues(qry_df, trgt_df, par, rtColname = c("r", "rt")), "Mis")
    expect_error(
        matchValues(qry_df, trgt_df, par, rtColname = c("rt", "r")), "Mis")

    res <- matchValues(qry_df, trgt_df, par)
    expect_equal(query(res), qry_df)
    expect_equal(target(res), trgt_df)
    expect_equal(res@matches$query_idx, c(2, 4))
    expect_equal(res@matches$target_idx, c(1, 1))
    expect_equal(res@matches$query_adduct, c("[M+H]+", "[M+K]+"))
    expect_equal(res@matches$target_adduct, rep("[M-H]-", 2))
    expect_equal(res@matches$score, c(0, 0))
    expect_equal(res@matches$ppm_error, c(0, 0))
    expect_equal(res@matches$score_rt, c(0, 0))
    expect_true(is(matchedData(res), "DataFrame"))

    ## positive toleranceRt
    par <- Mz2MassRtParam(queryAdducts = c("[M+H]+", "[M+K]+"),
                          targetAdducts = "[M-H]-",
                          toleranceRt = 2)
    res <- matchValues(qry_df, trgt_df, par)
    expect_equal(query(res), qry_df)
    expect_equal(target(res), trgt_df)
    expect_equal(res@matches$query_idx, c(2, 3, 4))
    expect_equal(res@matches$target_idx, c(1, 2, 1))
    expect_equal(res@matches$query_adduct, c("[M+H]+", "[M+H]+", "[M+K]+"))
    expect_equal(res@matches$target_adduct, rep("[M-H]-", 3))
    expect_equal(res@matches$score, c(0, 0, 0))
    expect_equal(res@matches$ppm_error, c(0, 0, 0))
    expect_equal(res@matches$score_rt, c(0, 1, 0))

    ## positive tolerance and toleranceRt
    par <- Mz2MassRtParam(queryAdducts = c("[M+H]+", "[M+K]+"),
                          targetAdducts = "[M-H]-", tolerance = 10,
                          toleranceRt = 2)
    res <- matchValues(qry_df, trgt_df, par)
    expect_equal(query(res), qry_df)
    expect_equal(target(res), trgt_df)
    expect_equal(res@matches$query_idx, c(2, 3, 4, 5))
    expect_equal(res@matches$target_idx, c(1, 2, 1, 2))
    expect_equal(res@matches$query_adduct, c(rep("[M+H]+", 2),
                                             rep("[M+K]+", 2)))
    expect_equal(res@matches$target_adduct, rep("[M-H]-", 4))
    expect_equal(res@matches$score, c(0, 0, 0, 5))
    expect_equal(res@matches$ppm_error, c(0, 0, 0, 5 / 300 * 10^6))
    expect_equal(res@matches$score_rt, c(0, 1, 0, 0))

    ## no matches
    trgt_df$rt <- trgt_df$rt + 5
    res <- matchValues(qry_df, trgt_df, par)
    expect_true(is(res, "Matched"))
    expect_equal(query(res), qry_df)
    expect_equal(target(res), trgt_df)
    expect_true(nrow(res@matches) == 0)
})

test_that(".valid_adduct works", {
    expect_equal(.valid_adduct(c("A", "[M+H]+"), name = "`targetAdducts`"),
                 paste("Unknown adducts in `targetAdducts` please check",
                       "MetaboCoreUtils for valid adducts"))
    expect_null(.valid_adduct(c("[M+H]+")))
    expect_equal(.valid_adduct(data.frame(a = 1, b = 2)),
                 paste("Columns \"mass_add\" and \"mass_multi\" must be",
                       "present when `adducts` is a data.frame"))
    expect_null(.valid_adduct(data.frame(mass_add = 1, mass_multi = 2)))
})
