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
    expect_error(matchValues(qry, trgt, par),
                 "`valueColname` has to be provided.")
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

    # no matches
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

    ## data.frame, SummarizedExperiment
    trgt_se <- SummarizedExperiment(
        assays = data.frame(matrix(NA, 10, 2,
                                   dimnames = list(NULL, c("A", "B")))),
        rowData = trgt)
    par <- ValueParam(tolerance = 0)

    expect_error(matchValues(qry, trgt_se, par),
                 "`valueColname` has to be provided.")
    expect_error(matchValues(qry, trgt_se, par,
                             valueColname = c("notcolumn", "value")), "Missing")
    expect_error(matchValues(qry, trgt_se, par,
                             valueColname = c("value", "notcolumn")), "Missing")
    res <- matchValues(qry, trgt_se, par, valueColname = "value")
    expect_equal(query(res), qry)
    expect_equal(target(res), trgt_se)
    expect_equal(res@matches$query_idx, c(1, 2))
    expect_equal(res@matches$target_idx, c(5, 7))
    expect_equal(res@matches$score, c(0, 0))
    expect_equal(res@matches$ppm_error, c(0, 0))

    # no matches
    res <- matchValues(qry + 0.1, trgt_se, par, valueColname = "value")
    expect_true(is(res, "Matched"))
    expect_equal(query(res), qry + 0.1)
    expect_equal(target(res), trgt_se)
    expect_true(nrow(res@matches) == 0)

    # positive tolerance
    par <- ValueParam(tolerance = 10)
    res <- matchValues(qry, trgt_se, par, valueColname = "value")
    expect_equal(query(res), qry)
    expect_equal(target(res), trgt_se)
    expect_equal(res@matches$query_idx, c(1, 1, 1, 2, 2, 2, 3, 3))
    expect_equal(res@matches$target_idx, c(4, 5, 6, 6, 7, 8, 7, 8))
    expect_equal(res@matches$score, c(10, 0, -10, 10, 0, -10, 9, -1))
    expect_equal(res@matches$ppm_error, c(10, 0, 10, 10, 0, 10, 9, 1) /
                     c(140, 150, 160, 160, 170, 180, 170, 180) * 10^6)

    ## data.frame, QFeatures
    trgt_qf <- QFeatures(list(a1 = trgt_se))
    par <- ValueParam(tolerance = 0)

    expect_error(matchValues(qry, trgt_qf, par), "`valueColname` has to be")
    expect_error(matchValues(qry, trgt_qf, par,
                             valueColname = c("notcolumn", "value")), "Missing")
    expect_error(matchValues(qry, trgt_qf, par, targetAssay = "a3",
                             valueColname = "value"), "No assay")
    expect_error(matchValues(qry, trgt_qf, par, targetAssay = "a1",
                             valueColname = c("value", "notcolumn")), "Missing")
    res <- matchValues(qry, trgt_qf, par, targetAssay = "a1",
                       valueColname = "value")
    expect_equal(query(res), qry)
    expect_equal(target(res), trgt_qf)
    expect_equal(res@matches$query_idx, c(1, 2))
    expect_equal(res@matches$target_idx, c(5, 7))
    expect_equal(res@matches$score, c(0, 0))
    expect_equal(res@matches$ppm_error, c(0, 0))

    # no matches
    res <- matchValues(qry + 0.1, trgt_qf, par, valueColname = "value",
                       targetAssay = "a1")
    expect_true(is(res, "Matched"))
    expect_equal(query(res), qry + 0.1)
    expect_equal(target(res), trgt_qf)
    expect_true(nrow(res@matches) == 0)

    # positive tolerance
    par <- ValueParam(tolerance = 10)
    res <- matchValues(qry, trgt_qf, par, valueColname = "value",
                       targetAssay = "a1")
    expect_equal(query(res), qry)
    expect_equal(target(res), trgt_qf)
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

    ## data.frame, data.frame

    # zero tolerance, positive ppm
    expect_error(matchValues(x, cmpds, par, mzColname = "mzz"), "Missing")
    res <- matchValues(x, cmpds, par)
    expect_equal(query(res), x)
    expect_equal(target(res), cmpds)
    expect_equal(res@matches$query_idx, c(1, 2, 2, 3))
    expect_equal(res@matches$target_idx, c(1, 2, 3, 1))
    expect_equal(res@matches$score, c(0, 0, 0, 1e-6))
    expect_equal(res@matches$ppm_error,
                 c(0, 0, 0, 1 / mass2mz(204.089878, "[M+Na]+")))

    # zero tolerance and ppm
    par <- Mass2MzParam(adducts = adducts, tolerance = 0, ppm = 0)
    res <- matchValues(x, cmpds, par)
    expect_equal(query(res), x)
    expect_equal(target(res), cmpds)
    expect_equal(res@matches$query_idx, c(1, 2, 2))
    expect_equal(res@matches$target_idx, c(1, 2, 3))
    expect_equal(res@matches$score, c(0, 0, 0))
    expect_equal(res@matches$ppm_error, c(0, 0, 0))

    # no matches
    adducts <- c("[M+Li]+", "[M+K]+")
    par <- Mass2MzParam(adducts = adducts, tolerance = 0, ppm = 20)
    res <- matchValues(x, cmpds, par)
    expect_true(is(res, "Matched"))
    expect_equal(query(res), x)
    expect_equal(target(res), cmpds)
    expect_true(nrow(res@matches) == 0)

    # with custom adduct definition, zero tolerance, positive ppm
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

    # with custom adduct definition, zero tolerance, positive ppm
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

    expect_error(matchValues(x, cmpds$exactmass, par, mzColname = "other"),
                 "other")

    # different adducts, no match
    adducts2 <- data.frame(mass_add = c(2, 4), mass_multi = c(2, 1))
    par <- Mass2MzParam(adducts = adducts2, tolerance = 0, ppm = 20)
    res <- matchValues(x, cmpds, par)
    expect_true(is(res, "Matched"))
    expect_equal(query(res), x)
    expect_equal(target(res), cmpds)
    expect_true(nrow(res@matches) == 0)

    ## SummarizedExperiment, data.frame
    par <- Mass2MzParam(adducts = adducts, tolerance = 0, ppm = 0)
    cmpds_se <- SummarizedExperiment(
        assays = data.frame(matrix(NA, 3, 2)),
        rowData = cmpds,
        colData = data.frame(cD1 = c(NA, NA),
                             cD2 = c(NA, NA)))
    res <- matchValues(x, cmpds_se, par)
    expect_equal(query(res), x)
    expect_equal(target(res), cmpds_se)
    expect_equal(res@matches$query_idx, c(1, 2, 2))
    expect_equal(res@matches$target_idx, c(1, 2, 3))
    expect_equal(res@matches$score, c(0, 0, 0))
    expect_equal(res@matches$ppm_error, c(0, 0, 0))

    ## data.frame, QFeatures
    cmpds_qf <- QFeatures(list(a1 = cmpds_se))
    res <- matchValues(x, cmpds_qf, par, targetAssay = "a1")
    expect_equal(query(res), x)
    expect_equal(target(res), cmpds_qf)
    expect_equal(res@matches$query_idx, c(1, 2, 2))
    expect_equal(res@matches$target_idx, c(1, 2, 3))
    expect_equal(res@matches$score, c(0, 0, 0))
    expect_equal(res@matches$ppm_error, c(0, 0, 0))
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
    adducts2 <- c("[M+Li]+", "[M+K]+")
    par <- Mass2MzRtParam(adducts = adducts2, tolerance = 0, ppm = 20,
                          toleranceRt = 0.2)
    res <- matchValues(x, cmpds, par)
    expect_true(is(res, "Matched"))
    expect_equal(query(res), x)
    expect_equal(target(res), cmpds)
    expect_true(nrow(res@matches) == 0)

    ## QFeatures, SummarizedExperiment
    par <- Mass2MzRtParam(adducts = adducts, tolerance = 0, ppm = 20,
                          toleranceRt = 0)

    x_qf <- QFeatures(list(a1 = SummarizedExperiment(assays = data.frame(1:3),
                                                     rowData = x,
                                                     colData = data.frame(1))))
    cmpds_se <- SummarizedExperiment(assays = data.frame(2:4), rowData = cmpds,
                                     colData = data.frame(2))
    res <- matchValues(x_qf, cmpds_se, par, queryAssay = "a1")
    expect_equal(query(res), x_qf)
    expect_equal(target(res), cmpds_se)
    expect_equal(res@matches$query_idx, c(1, 2, 2))
    expect_equal(res@matches$target_idx, c(1, 2, 3))
    expect_equal(res@matches$score, c(0, 0, 0))
    expect_equal(res@matches$ppm_error, c(0, 0, 0))
    expect_equal(res@matches$score_rt, c(0, 0, 0))


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

    ## numeric, QFeatures
    qry <- qry$mz
    trgt <- QFeatures(list(a1 = SummarizedExperiment(assays = data.frame(1:10),
                                                     rowData = trgt,
                                                     colData = data.frame(1))))

    par <- MzParam()
    res <- matchValues(qry, trgt, par, mzColname = "mz", targetAssay = "a1")
    expect_s4_class(res, "Matched")
    expect_error(matchValues(qry, trgt, par, mzColname = c("mzz", "mz"),
                             targetAssay = "a1"), "Missing")
    expect_error(matchValues(qry, trgt, par, mzColname = c("mz", "mzz"),
                             targetAssay = "a1"), "Missing")

    res <- matchValues(qry, trgt, par, targetAssay = "a1")
    expect_equal(query(res), qry)
    expect_equal(target(res), trgt)
    expect_equal(res@matches$query_idx, c(1, 2))
    expect_equal(res@matches$target_idx, c(5, 7))
    expect_equal(res@matches$score, c(0, 0))
    expect_equal(res@matches$ppm_error, c(0, 0))

    ## no matches
    res <- matchValues(qry + 0.1, trgt, par, targetAssay = "a1")
    expect_true(is(res, "Matched"))
    expect_equal(query(res), qry + 0.1)
    expect_equal(target(res), trgt)
    expect_true(nrow(res@matches) == 0)

    # positive tolerance
    par <- MzParam(tolerance = 10)
    res <- matchValues(qry, trgt, par, targetAssay = "a1")
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

    res <- matchValues(qry, trgt, par, mzColname = "mz")
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

    # no matches
    res <- matchValues(qry + 0.5, trgt, par)
    expect_true(is(res, "Matched"))
    expect_equal(query(res), qry + 0.5)
    expect_equal(target(res), trgt)
    expect_true(nrow(res@matches) == 0)

    ## `SummarizedExperiment`, data.frame
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

    # no matches
    res <- matchValues(qry , trgt + 0.5, par)
    expect_true(is(res, "Matched"))
    expect_equal(query(res), qry)
    expect_equal(target(res), trgt + 0.5)
    expect_true(nrow(res@matches) == 0)

    ## QFeatures
    ions <- data.frame(id = c("a", "b", "c", "d", "e", "f"),
                       mz = c(3, 4, 3, 5, 5, 6),
                       rt = c(4, 4, 8, 8, 9, 10))
    cmps <- data.frame(ID = c("AB", "CD", "EF"),
                       mz = c(3.5, 4, 5.5),
                       rt = c(4, 8, 9.5))
    trgt <- data.frame(id = c("X", "Y", "Z"),
                       mz = c(3.5, 8, 5.5),
                       rt = c(4.5, 8.7, 10))
    ions <- SummarizedExperiment(matrix(rnorm(12), nrow = 6, ncol = 2,
                                        dimnames = list(NULL, c("a", "b"))),
                                 rowData = ions)
    cmps <- SummarizedExperiment(matrix(rnorm(6), nrow = 3, ncol = 2,
                                        dimnames = list(NULL, c("a", "b"))),
                                 rowData = cmps)
    qf <- QFeatures(list(ions = ions, compounds = cmps))
    prm <- MzRtParam(tolerance = 0.5, toleranceRt = 0.5)
    res <- matchValues(qf, trgt, queryAssay = "ions",
                       param = prm)
    expect_true(all(colnames(rowData(ions)) %in% colnames(res)))
    expect_equal(res$score, c(-0.5, 0.5, NA, NA, NA, 0.5))
    expect_equal(res$score_rt, c(-0.5, -0.5, NA, NA, NA, 0))
    expect_error(matchValues(qf, trgt, param = prm), "queryAssay")
    expect_equal(matchedData(res), matchedData(matchValues(ions, trgt, prm)))

    res <- matchValues(qf, trgt, queryAssay = "compounds",
                       param = prm)
    expect_equal(res$score, c(0, NA, 0))
    expect_equal(res$score_rt, c(-0.5, NA, -0.5))
    expect_equal(matchedData(res), matchedData(matchValues(cmps, trgt, prm)))

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

    res <- matchValues(a, b, par, mzColname = "mz")
    expect_equal(query(res), a)
    expect_equal(target(res), b)
    expect_equal(res@matches$query_idx, c(2, 3, 4))
    expect_equal(res@matches$target_idx, c(1, 2, 1))

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

    # `query` and `target` data.frames
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

    # QFeatures, data.frame
    qry_se <- SummarizedExperiment(assays = data.frame(1:5), rowData = qry_df,
                                   colData = data.frame(a = 1, b = 2))
    qry_qf <- QFeatures(list(a1 = qry_se))
    par <- Mz2MassParam(queryAdducts = c("[M+H]+", "[M+K]+"),
                        targetAdducts = "[M-H]-",
                        tolerance = 10)
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

    ## numeric, numeric
    a <- 1:10
    b <- 1:5
    res <- matchValues(a, b, Mz2MassParam(queryAdducts = "[M]+",
                                          targetAdducts = "[M]-"))
    expect_equal(res@matches$target_idx, 1:5)
    expect_equal(res@matches$query_idx, 1:5)
    expect_true(all(res@matches$query_adduct == "[M]+"))
    expect_true(all(res@matches$target_adduct == "[M]-"))
    expect_s4_class(matchedData(res), "DataFrame")
    expect_equal(res$query, a)
    expect_equal(res$target, c(b, NA, NA, NA, NA, NA))
})

test_that("matchValues, Mz2MassRtParam works", {

    ## `query` and `target` data.frames
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
    md <- matchedData(res)
    expect_s4_class(md, "DataFrame")
    expect_true(all(c("mz", "rt", "target_mz", "target_rt",
                      "query_adduct", "target_adduct") %in% colnames(md)))

    # positive toleranceRt
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

    # positive tolerance and toleranceRt
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

    ## SummarizedExperiment, QFeatures
    qry_se <- SummarizedExperiment(assays = data.frame(1:5), rowData = qry_df,
                                   colData = data.frame(1))
    trgt_se <- SummarizedExperiment(assays = data.frame(1:4), rowData = trgt_df,
                         colData = data.frame(1))
    trgt_qf <- QFeatures(list(a1 = trgt_se))
    par <- Mz2MassRtParam(queryAdducts = c("[M+H]+", "[M+K]+"),
                          targetAdducts = "[M-H]-")

    res <- matchValues(qry_se, trgt_qf, par, mzColname = "mz",
                       targetAssay = "a1")
    expect_s4_class(res, "Matched")
    res <- matchValues(qry_se, trgt_qf, par, rtColname = "rt",
                       targetAssay = "a1")
    expect_s4_class(res, "Matched")

    expect_error(matchValues(qry_se, trgt_qf, par, mzColname = c("m", "mz"),
                             targetAssay = "a1"), "Mis")
    expect_error(matchValues(qry_se, trgt_qf, par, mzColname = c("mz", "m"),
                             targetAssay = "a1"), "Mis")
    expect_error(matchValues(qry_se, trgt_qf, par, rtColname = c("r", "rt"),
                             targetAssay = "a1"), "Mis")
    expect_error(matchValues(qry_se, trgt_qf, par, rtColname = c("rt", "r"),
                             targetAssay = "a1"), "Mis")

    res <- matchValues(qry_se, trgt_qf, par, targetAssay = "a1")
    expect_equal(query(res), qry_se)
    expect_equal(target(res), trgt_qf)
    expect_equal(res@matches$query_idx, c(2, 4))
    expect_equal(res@matches$target_idx, c(1, 1))
    expect_equal(res@matches$query_adduct, c("[M+H]+", "[M+K]+"))
    expect_equal(res@matches$target_adduct, rep("[M-H]-", 2))
    expect_equal(res@matches$score, c(0, 0))
    expect_equal(res@matches$ppm_error, c(0, 0))
    expect_equal(res@matches$score_rt, c(0, 0))
    md <- matchedData(res)
    expect_s4_class(md, "DataFrame")
    expect_true(all(c("mz", "rt", "target_mz", "target_rt",
                      "query_adduct", "target_adduct") %in% colnames(md)))

    # positive toleranceRt
    par <- Mz2MassRtParam(queryAdducts = c("[M+H]+", "[M+K]+"),
                          targetAdducts = "[M-H]-",
                          toleranceRt = 2)
    res <- matchValues(qry_se, trgt_qf, par, targetAssay = "a1")
    expect_equal(query(res), qry_se)
    expect_equal(target(res), trgt_qf)
    expect_equal(res@matches$query_idx, c(2, 3, 4))
    expect_equal(res@matches$target_idx, c(1, 2, 1))
    expect_equal(res@matches$query_adduct, c("[M+H]+", "[M+H]+", "[M+K]+"))
    expect_equal(res@matches$target_adduct, rep("[M-H]-", 3))
    expect_equal(res@matches$score, c(0, 0, 0))
    expect_equal(res@matches$ppm_error, c(0, 0, 0))
    expect_equal(res@matches$score_rt, c(0, 1, 0))

    # positive tolerance and toleranceRt
    par <- Mz2MassRtParam(queryAdducts = c("[M+H]+", "[M+K]+"),
                          targetAdducts = "[M-H]-", tolerance = 10,
                          toleranceRt = 2)
    res <- matchValues(qry_se, trgt_qf, par, targetAssay = "a1")
    expect_equal(query(res), qry_se)
    expect_equal(target(res), trgt_qf)
    expect_equal(res@matches$query_idx, c(2, 3, 4, 5))
    expect_equal(res@matches$target_idx, c(1, 2, 1, 2))
    expect_equal(res@matches$query_adduct, c(rep("[M+H]+", 2),
                                             rep("[M+K]+", 2)))
    expect_equal(res@matches$target_adduct, rep("[M-H]-", 4))
    expect_equal(res@matches$score, c(0, 0, 0, 5))
    expect_equal(res@matches$ppm_error, c(0, 0, 0, 5 / 300 * 10^6))
    expect_equal(res@matches$score_rt, c(0, 1, 0, 0))

    # no matches
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

test_that("matchValues,data.frame,Spectra,MzRtParam works", {
    df <- data.frame(mz = c(279.09, 200, 224.08, 100),
                     rt = c(379, 200, 378, 100))
    expect_error(matchValues(df, pest_ms2, MzRtParam(),
                             mzColname = c("mz", "mz")), "peak variable")
    expect_error(matchValues(df, pest_ms2, MzRtParam(),
                             mzColname = c("mz", "other")), "in target")
    expect_error(matchValues(df, pest_ms2, mzColname = c("mz", "precursorMz"),
                             MzRtParam(tolerance = 0.1,
                                       toleranceRt = 3)), "in target")
    res <- matchValues(df, pest_ms2, mzColname = c("mz", "precursorMz"),
                       MzRtParam(tolerance = 0.1, toleranceRt = 3),
                       rtColname = c("rt", "rtime"))
    expect_true(validObject(res))
    expect_equal(whichQuery(res), c(1L, 3L))
    expect_equal(whichTarget(res), c(4L, 9L, 7L, 11L))
    expect_equal(query(res), df)
    expect_equal(target(res), pest_ms2)

    res <- matchValues(df, pest_ms2, mzColname = c("mz", "precursorMz"),
                       rtColname = c("rt", "rtime"), param = MzRtParam())
    expect_equal(whichQuery(res), integer())
    expect_equal(whichTarget(res), integer())

    colnames(df) <- c("precursorMz", "rtime")
    res <- matchValues(df, pest_ms2, mzColname = c("precursorMz"),
                       MzRtParam(tolerance = 0.1, toleranceRt = 3),
                       rtColname = c("rtime"))
    expect_true(validObject(res))
    expect_equal(whichQuery(res), c(1L, 3L))
    expect_equal(whichTarget(res), c(4L, 9L, 7L, 11L))
    expect_equal(query(res), df)
    expect_equal(target(res), pest_ms2)

})

test_that("matchValues,data.frame,Spectra,MzParam works", {
    df <- data.frame(mz = c(279.09, 200, 224.08, 100),
                     rt = c(379, 200, 378, 100))
    expect_error(matchValues(df, pest_ms2, MzParam(),
                             mzColname = c("mz")), "peak variable")
    expect_error(matchValues(df, pest_ms2, MzParam(),
                             mzColname = c("mz", "other")), "in target")

    res <- matchValues(df, pest_ms2, MzParam(tolerance = 0.01),
                       mzColname = c("mz", "precursorMz"))

    expect_true(validObject(res))
    expect_equal(whichQuery(res), c(1L, 3L))
    expect_equal(whichTarget(res), c(4L, 9L, 7L, 11L))
    expect_equal(query(res), df)
    expect_equal(target(res), pest_ms2)
})

test_that("matchValues,numeric,Spectra,MzParam works", {
    mzs <- c(200, 400, 224.08, 124)

    expect_error(matchValues(mzs, pest_ms2, MzParam(),
                             mzColname = c("mz")), "peak variable")
    expect_error(matchValues(mzs, pest_ms2, MzParam(),
                             mzColname = c("other")), "in target")

    res <- matchValues(mzs, pest_ms2, MzParam(tolerance = 0.01),
                       mzColname = c("precursorMz"))

    expect_true(validObject(res))
    expect_equal(whichQuery(res), c(3L))
    expect_equal(whichTarget(res), c(7L, 11L))
    expect_equal(query(res), mzs)
    expect_equal(target(res), pest_ms2)
})
