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

test_that("matchSpectra,CompareSpectraParam works", {
    csp <- CompareSpectraParam(requirePrecursor = FALSE,
                               THRESHFUN = function(x) which.max(x))
    res <- matchSpectra(pest_ms2, minimb, csp)
    expect_equal(res@matches$query_idx, 1:13)

    csp <- CompareSpectraParam(requirePrecursor = TRUE,
                               THRESHFUN = function(x) which.max(x))
    res <- matchSpectra(pest_ms2, minimb, csp)
    expect_equal(res@matches$query_idx, c(2, 4, 6, 8, 9))
    expect_true(anyDuplicated(res@matches$query_idx) == 0)
    
    csp <- CompareSpectraParam(requirePrecursor = TRUE)
    res <- matchSpectra(pest_ms2, minimb, csp)
    expect_equal(unique(res@matches$query_idx), c(2, 4, 6, 8))
    expect_true(anyDuplicated(res@matches$query_idx) == 0)
})
