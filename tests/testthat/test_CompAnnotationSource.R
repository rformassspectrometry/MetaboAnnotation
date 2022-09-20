test_that("CompAnnotationSource works", {
    setClass("DummySource",
             contains = "CompAnnotationSource")
    d <- new("DummySource")
    expect_true(inherits(d, "CompAnnotationSource"))
    expect_equal(d@version, "0.1")
    expect_output(show(d), "DummySource")
})

test_that("matchSpectra,CompAnnotationSource works", {
    setClass("DummySource",
             contains = "CompAnnotationSource")
    d <- new("DummySource")
    expect_error(matchSpectra(pest_ms2, d, param = CompareSpectraParam()),
                 "Not implemented")
})

test_that("metadata,CompAnnotationSource works", {
    setClass("DummySource",
             contains = "CompAnnotationSource")
    d <- new("DummySource")
    expect_error(metadata(d), "Not implemented")
})
