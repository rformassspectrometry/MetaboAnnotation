test_that("CompareSpectraParam works", {
    res <- CompareSpectraParam()
    expect_true(is(res, "CompareSpectraParam"))

    expect_error(CompareSpectraParam(tolerance = 1:3), "positive number")
    expect_error(CompareSpectraParam(ppm = -4), "positive number")
    expect_error(CompareSpectraParam(toleranceRt = -1), "positive")
    expect_error(CompareSpectraParam(percentRt = c(-1, 1, 1)), "positive")

    res <- CompareSpectraParam(other_param = 5, b = 3)
    expect_equal(res@dots, list(other_param = 5, b = 3))

    res <- .compare_spectra_parms_list(res)
    expect_true(is.list(res))
    expect_equal(res$other_param, 5L)
    expect_equal(res$b, 3L)
    expect_equal(res$ppm, 5)
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

    expect_match(.valid_threshfun(function(z) -1L), "values between 1 and")
})

test_that(".compare_spectra, .compare_spectra_without_precursor work", {
    prm <- CompareSpectraParam(requirePrecursor = FALSE,
                               THRESHFUN = function(x) x > 0.3)

    res <- .match_spectra(pest_ms2, minimb, prm, BPPARAM = SerialParam())
    res_2 <- .match_spectra_without_precursor(pest_ms2, minimb, prm)
    expect_equal(res@matches, res_2@matches)

    res_3 <- .match_spectra(pest_ms2, minimb, prm, BPPARAM = SerialParam(),
                            rtColname = "rtime")
    expect_equal(res@matches, res_3@matches)

    ## returns integer.
    prm <- CompareSpectraParam(requirePrecursor = FALSE,
                               THRESHFUN = function(x) which.max(x))
    res <- .match_spectra(pest_ms2, minimb, prm, BPPARAM = SerialParam())
    res_2 <- .match_spectra_without_precursor(pest_ms2, minimb, prm)
    expect_equal(res@matches, res_2@matches)

    ## no matches
    prm <- CompareSpectraParam(requirePrecursor = FALSE,
                               THRESHFUN = function (x) which(x > 20))
    res <- .match_spectra(pest_ms2, minimb, prm, BPPARAM = SerialParam())
    res_2 <- .match_spectra_without_precursor(pest_ms2, minimb, prm)
    expect_equal(res@matches, res_2@matches)
})

test_that(".get_matches_spectra, matchSpectra,CompareSpectraParam works", {
    csp <- CompareSpectraParam(
        requirePrecursor = FALSE,
        THRESHFUN = function(x) which.max(x))
    res <- matchSpectra(pest_ms2, minimb, csp)
    expect_equal(res@matches$query_idx, 1:13)
    expect_equal(length(unique(res$target_spectrum_id)), 11)
    expect_true(any(spectraVariables(res) == ".original_query_index"))
    expect_equal(res@query$.original_query_index, seq_along(pest_ms2))
    expect_warning(matchSpectra(res@query, minimb, csp), "Overwriting")

    res <- matchSpectra(pest_ms2, minimb, csp, addOriginalQueryIndex = FALSE)
    expect_false(any(spectraVariables(res) == ".original_query_index"))

    mb2 <- minimb
    spectraNames(mb2) <- seq_along(mb2)
    res <- .get_matches_spectra(
        2, pest_ms2, mb2,
        .compare_spectra_parms_list(csp),
        csp@THRESHFUN, precMz = csp@requirePrecursor,
        precMzPeak = csp@requirePrecursorPeak, sn = spectraNames(mb2))
    expect_true(is.data.frame(res))
    expect_equal(res$query_idx, 2L)

    res <- .get_matches_spectra(
        2, pest_ms2, mb2, .compare_spectra_parms_list(csp),
        csp@THRESHFUN, precMz = FALSE, precMzPeak = TRUE,
        sn = spectraNames(mb2))
    expect_true(is.data.frame(res))
    expect_equal(res$query_idx, 2L)
    expect_equal(res$target_idx, 73L)

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

    csp <- CompareSpectraParam(toleranceRt = c(1, 4))
    expect_error(matchSpectra(pest_ms2, mb2, csp), "equal to the number")
    csp <- CompareSpectraParam(toleranceRt = c(1, 2, 10, 0, 0, 0, 4,
                                               Inf, 4, 0, 0, 0, 0))
    res <- matchSpectra(pest_ms2, mb2, csp)
    expect_equal(res@matches$target_idx, c(70, 73, 75, 47, 51, 53, 59))

    csp <- CompareSpectraParam(percentRt = 3, toleranceRt = 0)
    res <- matchSpectra(pest_ms2, mb2, csp)
    expect_true(all(res@matches$query_idx == 2))
    expect_equal(res@matches$target_idx, c(70, 73, 75))

    csp <- CompareSpectraParam(percentRt = c(5, 3, 0, 0.1, 0.1, 0.1, 0,
                                             10, 0.1, 0, 0.1, 0, 0.1),
                               toleranceRt = 0)
    res <- matchSpectra(pest_ms2, mb2, csp)
    expect_equal(res@matches$target_idx, c(70, 73, 75, 47, 51, 53, 59))

    ## Alternative retention time column.
    csp <- CompareSpectraParam()
    mb2 <- minimb
    mb2$rtime <- 361
    spectraNames(mb2) <- seq_along(mb2)
    res <- .get_matches_spectra(
        2, pest_ms2, mb2,
        .compare_spectra_parms_list(csp),
        csp@THRESHFUN, precMz = csp@requirePrecursor,
        precMzPeak = csp@requirePrecursorPeak, sn = spectraNames(mb2),
        toleranceRt = rep(1, length(pest_ms2)),
        percentRt = rep(0, length(pest_ms2)))
    expect_equal(res$query_idx, c(2, 2, 2))
    expect_equal(res$target_idx, c(70, 73, 75))

    expect_error(.get_matches_spectra(
        2, pest_ms2, mb2,
        .compare_spectra_parms_list(csp),
        csp@THRESHFUN, precMz = csp@requirePrecursor,
        precMzPeak = csp@requirePrecursorPeak, sn = spectraNames(mb2),
        toleranceRt = rep(1, length(pest_ms2)),
        percentRt = rep(0, length(pest_ms2)),
        query_rt_col = "other"), "not available")

    expect_error(.get_matches_spectra(
        2, pest_ms2, mb2,
        .compare_spectra_parms_list(csp),
        csp@THRESHFUN, precMz = csp@requirePrecursor,
        precMzPeak = csp@requirePrecursorPeak, sn = spectraNames(mb2),
        toleranceRt = rep(1, length(pest_ms2)),
        percentRt = rep(0, length(pest_ms2)),
        target_rt_col = "other"), "not available")

    pest_ms2$other <- pest_ms2$rtime
    mb2$other <- mb2$rtime
    mb2$rtime <- 900

    res <- .get_matches_spectra(
        2, pest_ms2, mb2,
        .compare_spectra_parms_list(csp),
        csp@THRESHFUN, precMz = csp@requirePrecursor,
        precMzPeak = csp@requirePrecursorPeak, sn = spectraNames(mb2),
        toleranceRt = rep(1, length(pest_ms2)),
        percentRt = rep(0, length(pest_ms2)),
        query_rt_col = "other", target_rt_col = "other")
    expect_equal(res$query_idx, c(2, 2, 2))
    expect_equal(res$target_idx, c(70, 73, 75))
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
    expect_true(any(spectraVariables(res) == ".original_query_index"))
    expect_equal(query(res)$.original_query_index, seq_along(pest_ms2))

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

    a <- data.frame(
        precursorMz = c(623.1618, 609.1825)
    )
    a$mz <- list(
        c(300.0276, 315.0511),
        c(242.0585, 286.0483, 301.0718)
    )
    a$intensity <- list(
        c(20, 100),
        c(1, 10, 100)
    )
    s <- Spectra(a)
    res <- matchSpectra(
        s[1], s[2],
        param =  MatchForwardReverseParam(
            requirePrecursor = FALSE, tolerance = 0.01,
            THRESHFUN = function(x) which(x >= 0),
            MAPFUN = joinPeaksGnps, FUN = MsCoreUtils::gnps))
    expect_equal(res@matches$matched_peaks_count, 2)
    ## Check that res@matches is not all 0

    res <- matchSpectra(
        s[1], s[2],
        param =  MatchForwardReverseParam(
            requirePrecursor = FALSE, tolerance = 1.1,
            THRESHFUN = function(x) which(x >= 0)))
    expect_equal(res@matches$matched_peaks_count, 1)
    expect_equal(res@matches$presence_ratio, 1/3)
})

test_that("matchSpectra,Spectra,CompDb works", {
    cdb <- new("CompDb")

    res <- matchSpectra(pest_ms2, cdb, CompareSpectraParam())
    expect_s4_class(res, "MatchedSpectra")
    expect_equal(length(res), 13)
    expect_s4_class(target(res), "Spectra")
    expect_equal(length(target(res)), 0)

    fl <- system.file("sql", "CompDb.MassBank.sql", package = "CompoundDb")
    cdb <- CompoundDb::CompDb(fl)
    res <- matchSpectra(pest_ms2, cdb, CompareSpectraParam())
    expect_s4_class(res, "MatchedSpectra")
    expect_equal(length(res), 13)
    expect_s4_class(target(res), "Spectra")
    expect_equal(length(target(res)), 70)

    expect_error(
        matchSpectra(pest_ms2, cdb, CompareSpectraParam(toleranceRt = 1),
                     rtColname = c("other", "other")), "not available")
})

test_that(".check_bpparam works", {
    res <- .check_bpparam(pest_ms2, minimb, SnowParam())
    expect_s4_class(res, "SnowParam")
})
