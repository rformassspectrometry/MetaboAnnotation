csp <- CompareSpectraParam(requirePrecursor = TRUE,
                           THRESHFUN = function(x) x >= 0.7)
norm_int <- function(x) {
    x[, "intensity"] <- x[, "intensity"] / max(x[, "intensity"]) * 100
    x
}
ms <- matchSpectra(addProcessing(pest_ms2, norm_int),
                   addProcessing(minimb, norm_int), csp)
ms2 <- ms[whichQuery(ms)]

test_that(".createChoices works", {
    res <- .createChoices(ms)
    expect_true(is.list(res))
    expect_equal(unname(res), as.list(seq_along(ms)))
})

test_that(".create_dt works", {
    res <- .create_dt(ms)
    expect_true(is.data.frame(res))
    expect_true(all(c("precursorMz", "target_precursorMz", "rtime", "score")
                    %in% colnames(res)))

    ms3 <- ms
    ms3@matches$reverse_score <- 3.2
    ms3@matches$presence_ratio <- 32.1
    res <- .create_dt(ms3)
    expect_true(all(res$reverse_score == 3.2))
    expect_true(all(res$presence_ratio == 32.1))
})

test_that("validateMatchedSpectra works", {
    with_mocked_bindings(
        ".is_shiny_available" = function() FALSE,
        code = expect_error(validateMatchedSpectra(ms), "requires package")
    )
    with_mocked_bindings(
        ".is_shinyjs_available" = function() FALSE,
        code = expect_error(validateMatchedSpectra(ms), "requires package")
    )
    with_mocked_bindings(
        ".is_dt_available" = function() FALSE,
        code = expect_error(validateMatchedSpectra(ms), "requires package")
    )

    expect_error(validateMatchedSpectra(3), "is not TRUE")
    expect_error(validateMatchedSpectra(MatchedSpectra()), "is empty")

    ## validateMatchedSpectra(ms)
})

test_that("packages available", {
    expect_true(.is_shiny_available())
    expect_true(.is_shinyjs_available())
    expect_true(.is_dt_available())
    expect_true(.is_plotly_available())
})

test_that("plotlySpectraMirror works", {
    with_mocked_bindings(
        ".is_plotly_available" = function() FALSE,
        code = expect_error(.plotlySpectraMirror(pest_ms2[1], pest_ms2[1]),
                            "'plotly'")
    )

    expect_error(.plotlySpectraMirror(3, pest_ms2[2]), "not TRUE")
    expect_error(.plotlySpectraMirror(pest_ms2[2], 3), "not TRUE")
    expect_error(.plotlySpectraMirror(pest_ms2, pest_ms2[1]), "length 1")
    expect_error(.plotlySpectraMirror(pest_ms2[1], pest_ms2), "length 1")
    p <- .plotlySpectraMirror(pest_ms2[1], pest_ms2[2])
    expect_true(is(p, "plotly"))
    p <- .plotlySpectraMirror(pest_ms2[6], pest_ms2[8], xLabel = "query",
                              yLabel = "target")

    ## One Spectra empty.
    p <- .plotlySpectraMirror(pest_ms2[2], Spectra())
    expect_true(is(p, "plotly"))
    p <- .plotlySpectraMirror(Spectra(), pest_ms2[2])
    expect_true(is(p, "plotly"))
})
