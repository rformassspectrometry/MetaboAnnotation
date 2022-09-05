test_that("CompDbSource works", {
    expect_error(new("CompDbSource", dbfile = tempfile()), "unable")
    fl <- system.file("extdata", "MS1_example.txt", package = "MetaboAnnotation")
    expect_error(new("CompDbSource", dbfile = fl), "database")

    res <- new("CompDbSource")
    expect_true(validObject(res))

    fl <- system.file("sql", "CompDb.MassBank.sql", package = "CompoundDb")
    res  <- new("CompDbSource", dbfile = fl)
    expect_true(validObject(res))

    expect_true(.validate_dbfile(fl))

    expect_output(show(res), "CompDbSource")
})

test_that("metadata,CompDbSource works", {
    fl <- system.file("sql", "CompDb.MassBank.sql", package = "CompoundDb")
    src  <- new("CompDbSource", dbfile = fl)
    res <- metadata(src)
    expect_true(is.data.frame(res))
})

test_that("matchSpectra,Spectra,CompDbSource works", {
    fl <- system.file("sql", "CompDb.MassBank.sql", package = "CompoundDb")
    src  <- new("CompDbSource", dbfile = fl)

    res <- matchSpectra(pest_ms2, src, param = CompareSpectraParam())
    expect_s4_class(res, "MatchedSpectra")
    expect_equal(query(res), pest_ms2)
    expect_s4_class(target(res)@backend, "MsBackendDataFrame")
    expect_true(length(target(res)) == 0)

    qry <- Spectra(CompoundDb::CompDb(fl))[3]
    res <- matchSpectra(qry, src, param = CompareSpectraParam())
    expect_true(length(target(res)) == 4)
    expect_equal(MetaboAnnotation::matches(res)$target_idx, 1:4)
    expect_s4_class(target(res)@backend, "MsBackendDataFrame")
})
