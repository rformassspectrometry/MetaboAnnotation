csp <- CompareSpectraParam(requirePrecursor = TRUE,
                           THRESHFUN = function(x) x >= 0.7)
norm_int <- function(x) {
    x[, "intensity"] <- x[, "intensity"] / max(x[, "intensity"]) * 100
    x
}
ms <- matchSpectra(addProcessing(pest_ms2, norm_int),
                   addProcessing(minimb, norm_int), csp)

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
})
