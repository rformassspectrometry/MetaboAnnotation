test_that("CompareSpectraParam works", {
    res <- CompareSpectraParam()
    expect_true(is(res, "CompareSpectraParam"))

    expect_error(CompareSpectraParam(tolerance = 1:3), "positive number")
    expect_error(CompareSpectraParam(ppm = -4), "positive number")
    expect_error(CompareSpectraParam(toleranceRt = 1:3), "positive number")

    res <- CompareSpectraParam(other_param = 5, b = 3)
    expect_equal(res@dots, list(other_param = 5, b = 3))

    res <- .compare_spectra_parms_list(res)
    expect_true(is.list(res))
    expect_equal(res$other_param, 5L)
    expect_equal(res$b, 3L)
    expect_equal(res$ppm, 5)
    expect_equal(res$toleranceRt, Inf)
})

test_that(".valid_threshfun works", {
    res <- .valid_threshfun(function(x) which(x > 0.4))
    expect_true(length(res) == 0)
    res <- .valid_threshfun(function(x) x > 0.4)
    expect_true(length(res) == 0)
    res <- .valid_threshfun(NULL)
    expect_true(length(res) == 0)

    res <- .valid_threshfun("f")
    expect_match(res, "function")

    res <- .valid_threshfun(function(z) z / 3)
    expect_match(res, "integer")

    res <- .valid_threshfun(function(z) -1)
    expect_match(res, "integer")

    res <- .valid_threshfun(function(z) TRUE)
    expect_match(res, "logical")
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
        THRESHFUN = function(x) which.max(x))
    res <- matchSpectra(pest_ms2, minimb, csp)
    expect_equal(res@matches$query_idx, 1:13)
    expect_equal(length(unique(res$target_spectrum_id)), 11)

    mb2 <- minimb
    spectraNames(mb2) <- seq_along(mb2)
    res <- .get_matches_spectra(
        2, pest_ms2, mb2,
        MetaboAnnotation:::.compare_spectra_parms_list(csp),
        csp@THRESHFUN, precMz = csp@requirePrecursor,
        precMzPeak = csp@requirePrecursorPeak, sn = spectraNames(mb2))
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

    csp <- CompareSpectraParam(toleranceRt = 10)
    res <- matchSpectra(pest_ms2, minimb, csp)
    expect_true(nrow(MetaboAnnotation::matches(res)) == 0)

    mb2 <- minimb
    mb2$rtime <- 350
    csp <- CompareSpectraParam(toleranceRt = Inf)
    res <- matchSpectra(pest_ms2, mb2, csp)
    csp <- CompareSpectraParam(toleranceRt = 100)
    res_2 <- matchSpectra(pest_ms2, mb2, csp)
    expect_equal(MetaboAnnotation::matches(res),
                 MetaboAnnotation::matches(res_2))
    mb2$rtime <- 361
    csp <- CompareSpectraParam(toleranceRt = 1)
    res <- matchSpectra(pest_ms2, mb2, csp)
    expect_true(all(res@matches$query_idx == 2))
    expect_equal(res@matches$target_idx, c(70, 73, 75))
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
                                          "reverse_score", "presence_ratio",
                                          "matched_peaks_count"))

    mp <- MatchForwardReverseParam(requirePrecursor = TRUE,
                                   THRESHFUN = function(x) which.max(x))
    res <- matchSpectra(pest_ms2, minimb, mp)
    expect_equal(res@matches$query_idx, c(2, 4, 6, 8, 9))
    expect_true(anyDuplicated(res@matches$query_idx) == 0)
    expect_true(all(res$reverse_score > res$score, na.rm = TRUE))

    mp <- MatchForwardReverseParam(requirePrecursor = TRUE)
    res <- matchSpectra(pest_ms2, minimb, mp)
    expect_equal(unique(res@matches$query_idx), c(2, 4, 6, 8))

    mp <- MatchForwardReverseParam(requirePrecursor = TRUE,
                                   THRESHFUN_REVERSE = function(z) z > 0.9)
    res_2 <- matchSpectra(pest_ms2, minimb, mp)
    expect_true(all(res_2@matches$reverse_score > 0.9))
    expect_equal(unique(res_2@matches$query_idx), c(2, 4, 8))
})
