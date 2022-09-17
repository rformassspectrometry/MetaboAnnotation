test_that(".get_weizmass_conf works", {
    res <- .get_weizmass_conf()
    expect_equal(names(res), c("user", "pass", "host", "dbname"))
})

test_that("RemoteWeizMassSource works", {
    expect_error(RemoteWeizMassSource(), "license")
    res <- RemoteWeizMassSource(TRUE)
    expect_s4_class(res@backend, "MsBackendWeizMass")
    expect_s4_class(res@drv, "MariaDBDriver")

    expect_output(show(res), "MsBackendWeizMass")
    expect_output(show(res), "Spectra source: WeizMass")
})

test_that("LocalWeizMassSource works", {
    expect_error(LocalWeizMassSource(), "license")
    res <- LocalWeizMassSource(TRUE)
    expect_s4_class(res@backend, "MsBackendWeizMass")
    expect_s4_class(res@drv, "SQLiteDriver")
})

test_that("matchSpectra,SpectraDbSource LocalWeizMass works", {
    ## only evaluate if MsBackendWeizMass is installed
    if (requireNamespace("MsBackendWeizMass", quietly = TRUE)) {
        fl <- system.file("sqlite", "weizmassv2.sqlite",
                          package = "MsBackendWeizMass")
        src <- LocalWeizMassSource(TRUE, dbfile = fl)
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

test_that("matchSpectra,SpectraDbSource RemoteWeizMass works", {
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
        src <- RemoteWeizMassSource(TRUE, dbname = "weiz_mass",
                                    host = "localhost")
        expect_equal(src@host, "localhost")
        expect_equal(src@dbname, "weiz_mass")
        expect_equal(src@user, character())
        expect_equal(src@pass, character())
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

        src <- RemoteWeizMassSource(TRUE, dbname = "weiz_mass",
                                    host = "localhost", user = "other",
                                    pass = "na")
        expect_equal(src@user, "other")
        expect_equal(src@pass, "na")

        expect_error(matchSpectra(pest_ms2, src, param = CompareSpectraParam()),
                     "Access denied")
    }
})
