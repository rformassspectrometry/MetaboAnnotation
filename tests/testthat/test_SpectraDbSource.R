test_that(".get_weizmass_conf works", {
    res <- .get_weizmass_conf()
    expect_equal(names(res), c("user", "pass", "host", "dbname"))
    expect_error(MetaboAnnotation:::.get_weizmass_conf("4"), "No configuration")
})

test_that("WeizMassSource works", {
    expect_error(WeizMassSource(), "license")
    if (requireNamespace("MsBackendWeizMass", quietly = TRUE) &&
        requireNamespace("RSQLite", quietly = TRUE)) {
        res <- WeizMassSource(TRUE, sqlite = TRUE, dbname = tempfile())
        expect_s4_class(res@backend, "MsBackendWeizMass")
        expect_s4_class(res@drv, "SQLiteDriver")
        expect_true(length(res@dbname) == 1L)
        expect_output(show(res), "MsBackendWeizMass")
        expect_output(show(res), "Spectra source: WeizMass")
    }
    if (requireNamespace("MsBackendWeizMass", quietly = TRUE) &&
        requireNamespace("RMariaDB", quietly = TRUE)) {
        res <- WeizMassSource(TRUE, host = "localhost")
        expect_s4_class(res@backend, "MsBackendWeizMass")
        expect_s4_class(res@drv, "MariaDBDriver")
        expect_equal(res@host, "localhost")
        expect_equal(res@user, character())
        expect_output(show(res), "MsBackendWeizMass")
        expect_output(show(res), "Spectra source: WeizMass")
    }
})

test_that(".sqlite_weiz_mass works", {
    if (requireNamespace("MsBackendWeizMass", quietly = TRUE) &&
        requireNamespace("RMariaDB", quietly = TRUE)) {
        res <- .sqlite_weiz_mass("dbname")
        expect_equal(res@dbname, "dbname")
        expect_s4_class(res@drv, "SQLiteDriver")
        expect_s4_class(res@backend, "MsBackendWeizMass")
        expect_error(.sqlite_weiz_mass(), "'dbname' is mandatory")
    }
})

test_that(".sql_weiz_mass works", {
    if (requireNamespace("MsBackendWeizMass", quietly = TRUE) &&
        requireNamespace("RMariaDB", quietly = TRUE)) {
        res <- .sql_weiz_mass(host = "localhost", version = "4")
        expect_equal(res@version, "4")
        expect_equal(res@user, character())
        expect_equal(res@host, "localhost")
        expect_s4_class(res@drv, "MariaDBDriver")
        expect_s4_class(res@backend, "MsBackendWeizMass")
    }
})

test_that("matchSpectra,SpectraDbSource SQLite WeizMass works", {
    ## only evaluate if MsBackendWeizMass is installed
    if (requireNamespace("MsBackendWeizMass", quietly = TRUE) &&
        requireNamespace("RSQLite", quietly = TRUE)) {
        fl <- system.file("sqlite", "weizmassv2.sqlite",
                          package = "MsBackendWeizMass")
        src <- WeizMassSource(TRUE, sqlite = TRUE, dbname = fl)
        res <- matchSpectra(pest_ms2, src, param = CompareSpectraParam())
        expect_s4_class(res, "MatchedSpectra")
        expect_equal(query(res), pest_ms2)
        expect_s4_class(target(res)@backend, "MsBackendDataFrame")
        expect_true(length(target(res)) == 0)

        be <- backendInitialize(MsBackendWeizMass::MsBackendWeizMass(),
                                DBI::dbConnect(RSQLite::SQLite(), fl))
        qry <- Spectra(be)[1L]
        res <- matchSpectra(qry, src, param = CompareSpectraParam())
        expect_s4_class(res, "MatchedSpectra")
        expect_equal(query(res), qry)
        expect_s4_class(target(res)@backend, "MsBackendDataFrame")
        expect_true(length(target(res)) == 1)
        expect_equal(mz(target(res)), mz(qry))
    }
})

test_that("matchSpectra,SpectraDbSource MySQL WeizMass works", {
    ## Only perform test if connection to database works, i.e. on a development
    ## machine.
    do_test <- FALSE
    if (requireNamespace("RMariaDB", quietly = TRUE)) {
        up <- MetaboAnnotation:::.get_weizmass_conf()
        con <- try(DBI::dbConnect(RMariaDB::MariaDB(), host = up["host"],
                                  dbname = up["dbname"], user = up["user"],
                                  pass = up["pass"]))
        if (!is(con, "try-error")) {
            DBI::dbDisconnect(con)
            do_test <- TRUE
        }
    }
    if (do_test && requireNamespace("MsBackendWeizMass", quietly = TRUE)) {
        src <- WeizMassSource(TRUE)
        expect_equal(src@host, character())
        expect_equal(src@dbname, character())
        expect_equal(src@user, character())
        expect_equal(src@pass, character())
        expect_equal(src@version, character())
        res <- matchSpectra(pest_ms2, src, param = CompareSpectraParam())
        expect_s4_class(res, "MatchedSpectra")
        expect_equal(query(res), pest_ms2)
        expect_s4_class(target(res)@backend, "MsBackendDataFrame")
        expect_true(length(target(res)) == 0)

        fl <- system.file("sqlite", "weizmassv2.sqlite",
                          package = "MsBackendWeizMass")
        be <- backendInitialize(MsBackendWeizMass::MsBackendWeizMass(),
                                DBI::dbConnect(RSQLite::SQLite(), fl))
        qry <- Spectra(be)[1L]
        res <- matchSpectra(qry, src, param = CompareSpectraParam())
        expect_s4_class(res, "MatchedSpectra")
        expect_equal(query(res), qry)
        expect_s4_class(target(res)@backend, "MsBackendDataFrame")
        expect_true(length(target(res)) == 1L)
        expect_equal(peaksData(qry), peaksData(target(res)))

        src <- WeizMassSource(TRUE, dbname = "weiz_mass", host = "localhost",
                              user = "other", pass = "na")
        expect_equal(src@user, "other")
        expect_equal(src@pass, "na")

        expect_error(matchSpectra(pest_ms2, src, param = CompareSpectraParam()),
                     "Access denied")
    }
})
