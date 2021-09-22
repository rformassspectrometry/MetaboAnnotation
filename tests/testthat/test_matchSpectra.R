test_that("CompareSpectraParam works", {
    res <- CompareSpectraParam()
    expect_true(is(res, "CompareSpectraParam"))

    expect_error(CompareSpectraParam(tolerance = 1:3), "positive number")
    expect_error(CompareSpectraParam(ppm = -4), "positive number")
    expect_error(CompareSpectraParam(threshold = 1:2), "length 1")

    res <- CompareSpectraParam(other_param = 5, b = 3)
    expect_equal(res@dots, list(other_param = 5, b = 3))

    res <- .compare_spectra_parms_list(res)
    expect_true(is.list(res))
    expect_equal(res$other_param, 5L)
    expect_equal(res$b, 3L)
    expect_equal(res$ppm, 5)
})

test_that(".compare_spectra, .compare_spectra_without_precursor work", {
    prm <- CompareSpectraParam(requirePrecursor = FALSE,
                               THRESHFUN = function(x) x > 0.3)

    res <- .match_spectra(pest_ms2, minimb, prm)
    res_2 <- .match_spectra_without_precursor(pest_ms2, minimb, prm)
    res_3 <- .match_spectra_parallel(pest_ms2, minimb, prm, SerialParam())
    expect_equal(res@matches, res_2@matches)
    expect_equal(res@matches, res_3@matches)

    ## returns integer.
    prm <- CompareSpectraParam(requirePrecursor = FALSE,
                               THRESHFUN = function(x) which.max(x))
    res <- .match_spectra(pest_ms2, minimb, prm)
    res_2 <- .match_spectra_without_precursor(pest_ms2, minimb, prm)
    res_3 <- .match_spectra_parallel(pest_ms2, minimb, prm, SerialParam())
    expect_equal(res@matches, res_2@matches)
    expect_equal(res@matches, res_3@matches)

    ## no matches
    prm <- CompareSpectraParam(requirePrecursor = FALSE,
                               THRESHFUN = function (x) which(x > 20))
    res <- .match_spectra(pest_ms2, minimb, prm)
    res_2 <- .match_spectra_without_precursor(pest_ms2, minimb, prm)
    res_3 <- .match_spectra_parallel(pest_ms2, minimb, prm, SerialParam())
    expect_equal(res@matches, res_2@matches)
    expect_equal(res@matches, res_3@matches)
})

test_that(".get_matches_spectra, matchSpectra,CompareSpectraParam works", {
    csp <- CompareSpectraParam(
        requirePrecursor = FALSE,
        THRESHFUN = function(x) x == apply(x, 1, max, na.rm = TRUE))
    res <- matchSpectra(pest_ms2, minimb, csp)
    expect_equal(res@matches$query_idx, 1:13)

    mb2 <- minimb
    spectraNames(mb2) <- seq_along(mb2)
    res <- .get_matches_spectra(
        2, pest_ms2, mb2,
        MetaboAnnotation:::.compare_spectra_parms_list(csp),
        csp@THRESHFUN, precMz = csp@requirePrecursor,
        precMzPeak = csp@requirePrecursorPeak, spectraNames(mb2))
    expect_true(is.data.frame(res))
    expect_equal(res$query_idx, 2L)

    csp <- CompareSpectraParam(requirePrecursor = TRUE,
                               THRESHFUN = function(x) which.max(x))
    res <- matchSpectra(pest_ms2, minimb, csp)
    expect_equal(res@matches$query_idx, c(2, 4, 6, 8, 9))
    expect_true(anyDuplicated(res@matches$query_idx) == 0)

    csp <- CompareSpectraParam(requirePrecursor = TRUE)
    res <- matchSpectra(pest_ms2, minimb, csp)
    expect_equal(unique(res@matches$query_idx), c(2, 4, 6, 8))
    expect_true(anyDuplicated(res@matches$query_idx) == 2)
})

test_that("MatchForwardReverseParam works", {
    res <- MatchForwardReverseParam()
    expect_true(is(res, "MatchForwardReverseParam"))
    expect_equal(res@FUN, MsCoreUtils::ndotproduct)
    expect_equal(res@requirePrecursor, TRUE)

    expect_warning(res <- MatchForwardReverseParam(type = "left"), "supported")
    expect_true(!any(names(res@dots) == "type"))
})

test_that("matchSpectra,MatchForwardReverseParam works", {
    mp <- MatchForwardReverseParam(requirePrecursor = FALSE,
                                   THRESHFUN = function(x) which.max(x))
    res <- matchSpectra(pest_ms2, minimb, mp)
    expect_equal(res@matches$query_idx, 1:13)
    expect_equal(colnames(res@matches), c("query_idx", "target_idx", "score",
                                          "reverse_score", "presence_ratio"))

    mp <- MatchForwardReverseParam(requirePrecursor = TRUE,
                                   THRESHFUN = function(x) which.max(x))
    res <- matchSpectra(pest_ms2, minimb, mp)
    expect_equal(res@matches$query_idx, c(2, 4, 6, 8, 9))
    expect_true(anyDuplicated(res@matches$query_idx) == 0)
    expect_true(all(res$reverse_score > res$score, na.rm = TRUE))

    mp <- MatchForwardReverseParam(requirePrecursor = TRUE)
    res <- matchSpectra(pest_ms2, minimb, mp)
    expect_equal(unique(res@matches$query_idx), c(2, 4, 6, 8))
})
