q1 <- data.frame(col1 = 1:5, col2 = 6:10)
t1 <- data.frame(col1 = 11:16, col2 = 17:22)
q2 <- list("A", c(1,2), 1L, list(), data.frame())
t2 <- list("A", c(1,2), 1L, list(), data.frame(), matrix(1:3))
library(SummarizedExperiment)
q3 <- SummarizedExperiment(assays = data.frame(matrix(NA, 5, 2)),
                           rowData = q1,
                           colData = data.frame(cD1 = c(NA, NA),
                                                cD2 = c(NA, NA)))
library(QFeatures)
q4 <- QFeatures(list(a1 = q3, a2 = q3[4:5]))
t3 <- t4 <- t1

m <- data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                target_idx = c(2L, 2L, 3L, 4L, 5L),
                score = seq(0.5, 0.9, by = 0.1))

test_that("validator work", {
    ## .validate_matches_format
    res <- .validate_matches_format(4)
    expect_match(res, "data.frame")
    res <- .validate_matches_format(data.frame())
    expect_match(res, "Not all ")
    tmp <- data.frame(query_idx = 1:3, target_idx = 1:6, score = 12)
    res <- .validate_matches_format(tmp)
    expect_equal(res, NULL)
    tmp$query_idx <- "a"
    res <- .validate_matches_format(tmp)
    expect_match(res, "query_idx")
    tmp$query_idx <- 1:6
    tmp$target_idx <- "a"
    res <- .validate_matches_format(tmp)
    expect_match(res, "target_idx")

    ## .validate_matches_content
    tmp$target_idx <- c(1:3, 1:3)
    res <- .validate_matches_content(tmp, 3, 3)
    expect_match(res, "query_idx")
    res <- .validate_matches_content(tmp, 9, 2)
    expect_match(res, "target_idx")
    res <- .validate_matches_content(tmp, 10, 10)
    expect_equal(res, NULL)

    ## .validate_qt
    q <-  matrix(rep(1,4),2,2)
    res <- .validate_qt(q)
    expect_match(res, "no column names")
    colnames(q) <- c("col1", "col2")
    res <- .validate_qt(q)
    expect_equal(res, NULL)
    q <- array(rep(1,8), dim = c(2,2,2))
    res <- .validate_qt(q)
    expect_match(res, "unsupported dimensions ")

})

test_that("Matched works", {
    mo <- Matched()
    expect_true(validObject(mo))
    expect_output(show(mo), "Matched")
    expect_equal(whichQuery(mo), integer())
    expect_equal(whichTarget(mo), integer())

    ## With data. query and target of type data.frame
    mo <- Matched(query = q1, target = t1,
                  matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                       target_idx = c(2L, 2L, 3L, 4L, 5L),
                                       score = seq(0.5, 0.9, by = 0.1)))
    expect_true(validObject(mo))
    expect_true(length(mo) == 5)
    expect_output(show(mo), "4 matched")
    expect_true(is(query(mo), "data.frame"))
    expect_true(is(target(mo), "data.frame"))
    expect_equal(whichQuery(mo), c(1L, 2L, 5L))
    expect_equal(whichTarget(mo), c(2L, 3L, 4L, 5L))

    ## With data. query and target of type list
    mo <- Matched(query = q2, target = t2,
                  matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                       target_idx = c(2L, 2L, 3L, 4L, 5L),
                                       score = seq(0.5, 0.9, by = 0.1)))
    expect_true(validObject(mo))
    expect_true(length(mo) == 5)
    expect_output(show(mo), "4 matched")
    expect_true(is(query(mo), "list"))
    expect_true(is(target(mo), "list"))
    expect_equal(whichQuery(mo), c(1L, 2L, 5L))
    expect_equal(whichTarget(mo), c(2L, 3L, 4L, 5L))

    ## With data. query SummarizedExperiment and target data.frame
    mo <- Matched(query = q3, target = t3,
                  matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                       target_idx = c(2L, 2L, 3L, 4L, 5L),
                                       score = seq(0.5, 0.9, by = 0.1)))
    expect_true(validObject(mo))
    expect_true(length(mo) == 5)
    expect_output(show(mo), "4 matched")
    expect_true(is(query(mo), "SummarizedExperiment"))
    expect_true(is(target(mo), "data.frame"))
    expect_equal(whichQuery(mo), c(1L, 2L, 5L))
    expect_equal(whichTarget(mo), c(2L, 3L, 4L, 5L))

})

test_that(".subset_matches_nodim and [ works", {
    mo <- Matched()
    expect_error(.subset_matches_nodim(mo, 1), "out-of-bounds")

    #### query and target: data.frames
    mo <- Matched(
        q1, t1, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                     target_idx = c(2L, 2L, 3L, 4L, 5L),
                                     score = seq(0.5, 0.9, by = 0.1)))

    res <- .subset_matches_nodim(mo, 2)
    expect_true(length(res) == 1)
    expect_equal(res@matches$query_idx, c(1L, 1L, 1L))
    expect_equal(res@matches$target_idx, mo@matches$target_idx[c(2, 3, 4)])
    expect_equal(res@matches$score, mo@matches$score[c(2, 3, 4)])
    expect_equal(res@query, mo@query[2, ])
    expect_equal(res@target, mo@target)

    ## No matching
    res <- .subset_matches_nodim(mo, c(3, 4))
    expect_true(validObject(res))
    expect_true(nrow(res@matches) == 0)
    expect_equal(res@query, mo@query[c(3, 4), ])

    expect_error(.subset_matches_nodim(mo, 12), "out-of-bounds")

    res <- .subset_matches_nodim(mo, c(1, 5))
    expect_equal(res@query, mo@query[c(1, 5), ])
    expect_true(length(res) == 2)
    expect_equal(res@matches$score, mo@matches$score[c(1, 5)])

    ## arbitrary order
    res <- .subset_matches_nodim(mo, c(3, 2, 1, 1, 4))
    expect_equal(query(res), query(mo)[c(3, 2, 1, 1, 4), ])
    expect_equal(target(res), target(mo))
    expect_equal(res@matches$query_idx, c(2L, 2L, 2L, 3L, 4L))
    expect_equal(res@matches$score, c(0.6, 0.7, 0.8, 0.5, 0.5))

    res <- mo[]
    expect_equal(res, mo)
    res <- mo[c(TRUE, FALSE)]
    expect_equal(res@matches$query_idx, 1L)
    expect_equal(res@matches$target_idx, 2L)
    expect_equal(res@matches$score, 0.5)

    ## All works even after subsetting and pruning.
    mo <- Matched(
        q1, t1, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                     target_idx = c(2L, 2L, 3L, 4L, 5L),
                                     score = seq(0.5, 0.9, by = 0.1)))

    res <- mo[whichQuery(mo)]
    expect_equal(query(res), q1[whichQuery(mo), ])
    res <- pruneTarget(res)
    expect_equal(target(res), t1[c(2, 3, 4, 5), ])
    expect_equal(res@matches$target_idx, c(1L, 1L, 2L, 3L, 4L))

    #### query and target: lists
    mo <- Matched(
        q2, t2, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                     target_idx = c(2L, 2L, 3L, 4L, 5L),
                                     score = seq(0.5, 0.9, by = 0.1)))
    res <- .subset_matches_nodim(mo, 2)
    expect_true(length(res) == 1)
    expect_equal(res@matches$query_idx, c(1L, 1L, 1L))
    expect_equal(res@matches$target_idx, mo@matches$target_idx[c(2, 3, 4)])
    expect_equal(res@matches$score, mo@matches$score[c(2, 3, 4)])
    expect_equal(res@query, mo@query[2])
    expect_equal(res@target, mo@target)

    expect_error(.subset_matches_nodim(mo, 12), "out-of-bounds")

    res <- .subset_matches_nodim(mo, c(1, 5))
    expect_equal(res@query, mo@query[c(1, 5)])
    expect_true(length(res) == 2)
    expect_equal(res@matches$score, mo@matches$score[c(1, 5)])

    ## arbitrary order
    res <- .subset_matches_nodim(mo, c(3, 2, 1, 1, 4))
    expect_equal(query(res), query(mo)[c(3, 2, 1, 1, 4)])
    expect_equal(target(res), target(mo))
    expect_equal(res@matches$query_idx, c(2L, 2L, 2L, 3L, 4L))
    expect_equal(res@matches$score, c(0.6, 0.7, 0.8, 0.5, 0.5))

    res <- mo[]
    expect_equal(res, mo)
    res <- mo[c(TRUE, FALSE)]
    expect_equal(res@matches$query_idx, 1L)
    expect_equal(res@matches$target_idx, 2L)
    expect_equal(res@matches$score, 0.5)

    ## All works even after subsetting and pruning.
    mo <- Matched(
        q2, t2, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                     target_idx = c(2L, 2L, 3L, 4L, 5L),
                                     score = seq(0.5, 0.9, by = 0.1)))

    res <- mo[whichQuery(mo)]
    expect_equal(query(res), q2[whichQuery(mo)])
    res <- pruneTarget(res)
    expect_equal(target(res), t2[c(2, 3, 4, 5)])
    expect_equal(res@matches$target_idx, c(1L, 1L, 2L, 3L, 4L))

    #### query SummarizedExperiment and target data.frame
    mo <- Matched(
        q3, t3, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                     target_idx = c(2L, 2L, 3L, 4L, 5L),
                                     score = seq(0.5, 0.9, by = 0.1)))

    res <- .subset_matches_nodim(mo, 2)
    expect_true(length(res) == 1)
    expect_equal(res@matches$query_idx, c(1L, 1L, 1L))
    expect_equal(res@matches$target_idx, mo@matches$target_idx[c(2, 3, 4)])
    expect_equal(res@matches$score, mo@matches$score[c(2, 3, 4)])
    expect_equal(res@query, mo@query[2, ])
    expect_equal(res@target, mo@target)


    expect_error(.subset_matches_nodim(mo, 12), "out-of-bounds")

    res <- .subset_matches_nodim(mo, c(1, 5))
    expect_equal(res@query, mo@query[c(1, 5), ])
    expect_true(length(res) == 2)
    expect_equal(res@matches$score, mo@matches$score[c(1, 5)])

    ## arbitrary order
    res <- .subset_matches_nodim(mo, c(3, 2, 1, 1, 4))
    expect_equal(query(res), query(mo)[c(3, 2, 1, 1, 4), ])
    expect_equal(target(res), target(mo))
    expect_equal(res@matches$query_idx, c(2L, 2L, 2L, 3L, 4L))
    expect_equal(res@matches$score, c(0.6, 0.7, 0.8, 0.5, 0.5))

    res <- mo[]
    expect_equal(res, mo)
    res <- mo[c(TRUE, FALSE)]
    expect_equal(res@matches$query_idx, 1L)
    expect_equal(res@matches$target_idx, 2L)
    expect_equal(res@matches$score, 0.5)

    ## All works even after subsetting and pruning.
    mo <- Matched(
        q3, t3, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                     target_idx = c(2L, 2L, 3L, 4L, 5L),
                                     score = seq(0.5, 0.9, by = 0.1)))

    res <- mo[whichQuery(mo)]
    expect_equal(query(res), q3[whichQuery(mo), ])
    res <- pruneTarget(res)
    expect_equal(target(res), t3[c(2, 3, 4, 5), ])
    expect_equal(res@matches$target_idx, c(1L, 1L, 2L, 3L, 4L))
})

test_that(".fill_index works", {
    res <- .fill_index(1:5, integer())
    expect_equal(res, 1:5)

    res <- .fill_index(1:20, c(2, 4, 4, 4, 6, 9, 10))
    expect_equal(res, c(1, 2, 3, 4, 4, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
                        14, 15, 16, 17, 18, 19, 20))
})

test_that("$,Matched and .dollar works", {
    mo <- Matched(
        q1, t1, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                     target_idx = c(2L, 2L, 3L, 4L, 5L),
                                     score = seq(0.5, 0.9, by = 0.1)))
    expect_equal(mo$score, c(0.5, 0.6, 0.7, 0.8, NA, NA, 0.9))
    expect_equal(mo$score, .dollar(mo@query, mo@target, mo@matches, "score"))
    expect_equal(mo$col1, q1[c(1, 2, 2, 2, 3, 4, 5), "col1"])
    expect_equal(mo$col1, .dollar(mo@query, mo@target, mo@matches, "col1"))
    expect_equal(mo$target_col1, t1[c(2,2,3,4, NA, NA, 5), "col1"])
    expect_equal(mo$target_col1,
                 .dollar(mo@query, mo@target, mo@matches, "target_col1"))

    mo <- Matched(
        q2, t2, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                     target_idx = c(2L, 2L, 3L, 4L, 5L),
                                     score = seq(0.5, 0.9, by = 0.1)))
    expect_equal(mo$score, c(0.5, 0.6, 0.7, 0.8, NA, NA, 0.9))
    expect_equal(mo$query, q2[c(1, 2, 2, 2, 3, 4, 5)])
    tmp <- t2[c(2,2,3,4, NA, NA, 5)]
    tmp[5:6] <- NA
    expect_equal(mo$target, tmp)

    mo <- Matched(
        q2, t1, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                     target_idx = c(2L, 2L, 3L, 4L, 5L),
                                     score = seq(0.5, 0.9, by = 0.1)))
    expect_equal(mo$score, c(0.5, 0.6, 0.7, 0.8, NA, NA, 0.9))
    expect_equal(mo$query, q2[c(1, 2, 2, 2, 3, 4, 5)])
    expect_equal(mo$target_col1, t1[c(2,2,3,4, NA, NA, 5), "col1"])

    mo <- Matched(
        q1, t2, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                     target_idx = c(2L, 2L, 3L, 4L, 5L),
                                     score = seq(0.5, 0.9, by = 0.1)))
    expect_equal(mo$score, c(0.5, 0.6, 0.7, 0.8, NA, NA, 0.9))
    expect_equal(mo$col1, q1[c(1, 2, 2, 2, 3, 4, 5), "col1"])
    tmp <- t2[c(2,2,3,4, NA, NA, 5)]
    tmp[5:6] <- NA
    expect_equal(mo$target, tmp)

    mo <- Matched(
        q3, t3, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                     target_idx = c(2L, 2L, 3L, 4L, 5L),
                                     score = seq(0.5, 0.9, by = 0.1)))
    expect_equal(mo$score, c(0.5, 0.6, 0.7, 0.8, NA, NA, 0.9))
    expect_equal(mo$col1, rowData(q3)[c(1, 2, 2, 2, 3, 4, 5), "col1"])
    expect_equal(mo$target_col1, t3[c(2,2,3,4, NA, NA, 5), "col1"])
})

test_that(".dollar works with unordered @matches", {
    mo <- Matched(
        q1, t1, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                     target_idx = c(2L, 2L, 3L, 4L, 5L),
                                     score = seq(0.5, 0.9, by = 0.1)))
    ## Both query and target idx are unsorted.
    mo2 <- mo
    mo2@matches <- mo@matches[c(2, 1, 3, 5, 4), ]
    ## query
    expect_equal(mo$col1, mo2$col1)
    expect_equal(mo$col2, mo2$col2)
    ## matches
    expect_equal(mo$score, mo2$score)
    ## target
    expect_equal(mo$target_col1, mo2$target_col1)
    expect_equal(mo$target_col2, mo2$target_col2)

    ## query idx is sorted, target idx not.
    mo2 <- mo
    mo2@target <- mo@target[c(4, 5, 1, 2, 3, 6), ]
    mo2@matches$target_idx <- c(4L, 4L, 5L, 1L, 2L)

    ## matches
    expect_equal(mo$score, mo2$score)
    ## target
    expect_equal(mo$target_col1, mo2$target_col1)
    expect_equal(mo$target_col2, mo2$target_col2)
})

test_that("colnames,Matched works", {
    mo <- Matched(
        q1, t1, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                     target_idx = c(2L, 2L, 3L, 4L, 5L),
                                     score = seq(0.5, 0.9, by = 0.1)))
    expect_equal(colnames(mo), c(colnames(q1),
                                 paste0("target_", colnames(t1)), "score"))
    mo <- Matched(
        q2, t2, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                     target_idx = c(2L, 2L, 3L, 4L, 5L),
                                     score = seq(0.5, 0.9, by = 0.1)))
    expect_equal(colnames(mo), c("query", "target", "score"))

    mo <- Matched(
        q1, t2, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                     target_idx = c(2L, 2L, 3L, 4L, 5L),
                                     score = seq(0.5, 0.9, by = 0.1)))
    expect_equal(colnames(mo), c(colnames(q1), "target", "score"))

    mo <- Matched(
        q2, t1, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                     target_idx = c(2L, 2L, 3L, 4L, 5L),
                                     score = seq(0.5, 0.9, by = 0.1)))
    expect_equal(colnames(mo), c("query",
                                 paste0("target_", colnames(t1)), "score"))

    mo <- Matched(
        q3, t3, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                     target_idx = c(2L, 2L, 3L, 4L, 5L),
                                     score = seq(0.5, 0.9, by = 0.1)))
    expect_equal(colnames(mo), c(colnames(rowData(q3)),
                                 paste0("target_", colnames(t3)), "score"))
})

test_that("matchedData,Matched works", {

    mo <- Matched()
    res <- matchedData(mo)
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) == 0)
    expect_equal(colnames(res), colnames(mo))

    #### query and target data.frames
    mo <- Matched(
        q1, t1, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                     target_idx = c(2L, 2L, 3L, 4L, 5L),
                                     score = seq(0.5, 0.9, by = 0.1)))
    res <- matchedData(mo)
    expect_equal(res$col2, mo$col2)
    expect_equal(res$target_col2, mo$target_col2)
    expect_equal(res$score, mo$score)

    expect_error(matchedData(mo, columns = "other"), "other not available")

    ## Only query object variables
    res <- matchedData(mo, columns = c("col1", "col2"))
    expect_equal(colnames(res), c("col1", "col2"))
    expect_equal(res$col1, mo$col1)

    ## Only target object variables
    res <- matchedData(mo, columns = c("target_col1", "target_col2"))
    expect_equal(colnames(res), c("target_col1", "target_col2"))
    expect_equal(res$target_col1, mo$target_col1)

    ## Only matches
    res <- matchedData(mo, columns = c("score"))
    expect_equal(colnames(res), c("score"))
    expect_equal(res$score, mo$score)

    ## A Matched with no matching target elements (query and target data.frames)
    mo <- Matched(q1, t1, matches = data.frame(query_idx = integer(),
                                               target_idx = integer(),
                                               score = numeric()))
    res <- matchedData(mo)
    expect_equal(res$col2, mo$col2)
    expect_equal(res$target_col2, mo$target_col2)
    expect_equal(res$score, mo$score)

    ## Only query object variables
    res <- matchedData(mo, columns = c("col1", "col2"))
    expect_equal(colnames(res), c("col1", "col2"))
    expect_equal(res$col1, mo$col1)

    ## Only target object variables
    res <- matchedData(mo, columns = c("target_col1", "target_col2"))
    expect_equal(colnames(res), c("target_col1", "target_col2"))
    expect_equal(res$target_col1, mo$target_col1)

    ## Only matches
    res <- matchedData(mo, columns = c("score"))
    expect_equal(colnames(res), c("score"))
    expect_equal(res$score, mo$score)

    #### query data.frame and target list
    mo <- Matched(
        q1, t2, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                     target_idx = c(2L, 2L, 3L, 4L, 5L),
                                     score = seq(0.5, 0.9, by = 0.1)))
    res <- matchedData(mo)
    expect_equal(res$col2, mo$col2)
    expect_equal(res$target, I(mo$target))
    expect_equal(res$score, mo$score)

    ## Only query object variables
    res <- matchedData(mo, columns = c("col1", "col2"))
    expect_equal(colnames(res), c("col1", "col2"))
    expect_equal(res$col1, mo$col1)

    ## Only target object variables
    res <- matchedData(mo, columns = c("target"))
    expect_equal(colnames(res), c("target"))
    expect_equal(res$target, mo$target)

    ## Only matches
    res <- matchedData(mo, columns = c("score"))
    expect_equal(colnames(res), c("score"))
    expect_equal(res$score, mo$score)

    #### A Matched with no matching target elements: query list and target
    #### data.frame
    mo <- Matched(q2, t1, matches = data.frame(query_idx = integer(),
                                               target_idx = integer(),
                                               score = numeric()))
    res <- matchedData(mo)
    expect_equal(res$query, I(mo$query))
    expect_equal(res$target_col1, mo$target_col1)
    expect_equal(res$score, mo$score)

    ## Only query object variables
    res <- matchedData(mo, columns = c("query"))
    expect_equal(colnames(res), c("query"))
    expect_equal(res$query, mo$query)
    expect_equal(res$query, q2)

    ## Only target object variables
    res <- matchedData(mo, columns = c("target_col1", "target_col2"))
    expect_equal(colnames(res), c("target_col1", "target_col2"))
    expect_equal(res$target_col1, mo$target_col1)

    ## Only matches
    res <- matchedData(mo, columns = c("score"))
    expect_equal(colnames(res), c("score"))
    expect_equal(res$score, mo$score)

    #### query SummarizedExoeriment and target data.frame
    mo <- Matched(
        q3, t3, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                     target_idx = c(2L, 2L, 3L, 4L, 5L),
                                     score = seq(0.5, 0.9, by = 0.1)))
    res <- matchedData(mo)
    expect_equal(res$col2, mo$col2)
    expect_equal(res$target_col2, mo$target_col2)
    expect_equal(res$score, mo$score)

    expect_error(matchedData(mo, columns = "other"), "other not available")

    ## Only query object variables
    res <- matchedData(mo, columns = c("col1", "col2"))
    expect_equal(colnames(res), c("col1", "col2"))
    expect_equal(res$col1, mo$col1)

    ## Only target object variables
    res <- matchedData(mo, columns = c("target_col1", "target_col2"))
    expect_equal(colnames(res), c("target_col1", "target_col2"))
    expect_equal(res$target_col1, mo$target_col1)

    ## Only matches
    res <- matchedData(mo, columns = c("score"))
    expect_equal(colnames(res), c("score"))
    expect_equal(res$score, mo$score)

    ## A MatchedSummarizedExperiment with no matching target elements
    mo <- Matched(q3, t3, matches = data.frame(query_idx = integer(),
                                               target_idx = integer(),
                                               score = numeric()))
    res <- matchedData(mo)
    expect_equal(res$col2, mo$col2)
    expect_equal(res$target_col2, mo$target_col2)
    expect_equal(res$score, mo$score)

    ## Only query object variables
    res <- matchedData(mo, columns = c("col1", "col2"))
    expect_equal(colnames(res), c("col1", "col2"))
    expect_equal(res$col1, mo$col1)

    ## Only target object variables
    res <- matchedData(mo, columns = c("target_col1", "target_col2"))
    expect_equal(colnames(res), c("target_col1", "target_col2"))
    expect_equal(res$target_col1, mo$target_col1)

    ## Only matches
    res <- matchedData(mo, columns = c("score"))
    expect_equal(colnames(res), c("score"))
    expect_equal(res$score, mo$score)
})

test_that("pruneTarget,Matched works", {
    mo <- Matched()
    res <- pruneTarget(mo)
    expect_true(is(res, "Matched"))
    expect_equal(res, mo)

    #### query and target data.frames
    mo <- Matched(
        q1, t1, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                     target_idx = c(2L, 2L, 3L, 4L, 5L),
                                     score = seq(0.5, 0.9, by = 0.1)))
    res <- pruneTarget(mo)
    expect_equal(matchedData(res), matchedData(mo))
    expect_true(nrow(res@target) < nrow(mo@target))

    mo <- Matched(q1, t1, matches = data.frame(query_idx = integer(),
                                               target_idx = integer(),
                                               score = numeric()))
    res <- pruneTarget(mo)
    expect_equal(matchedData(res), matchedData(mo))
    expect_true(nrow(res@target) < nrow(mo@target))

    #### query data.frame and target list
    mo <- Matched(
        q1, t2, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                     target_idx = c(2L, 2L, 3L, 4L, 5L),
                                     score = seq(0.5, 0.9, by = 0.1)))
    res <- pruneTarget(mo)
    expect_equal(matchedData(res), matchedData(mo))
    expect_true(length(res@target) < length(mo@target))

    mo <- Matched(q1, t2, matches = data.frame(query_idx = integer(),
                                               target_idx = integer(),
                                               score = numeric()))
    res <- pruneTarget(mo)
    expect_equal(matchedData(res), matchedData(mo))
    expect_true(length(res@target) < length(mo@target))

    #### query SummarizedExperiment and target data.frame
    mo <- Matched(
        q3, t3, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                     target_idx = c(2L, 2L, 3L, 4L, 5L),
                                     score = seq(0.5, 0.9, by = 0.1)))
    res <- pruneTarget(mo)
    expect_equal(matchedData(res), matchedData(mo))
    expect_true(nrow(res@target) < nrow(mo@target))

    mo <- Matched(q3, t3, matches = data.frame(query_idx = integer(),
                                               target_idx = integer(),
                                               score = numeric()))
    res <- pruneTarget(mo)
    expect_equal(matchedData(res), matchedData(mo))
    expect_true(nrow(res@target) < nrow(mo@target))
})

test_that("filterMatches,Matched works", {
    #### query and target data.frames
    mo <- Matched(
        q1, t1, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                     target_idx = c(2L, 2L, 3L, 4L, 5L),
                                     score = seq(0.5, 0.9, by = 0.1)))
    ## out of bounds indexes
    expect_error(filterMatches(mo, index = c(1, 10)), "out-of-bounds")
    ## no index : every match is removed
    idxs <- integer(0)
    mosub <- filterMatches(mo, index = idxs)
    expect_equal(mosub@matches, mo@matches[idxs, ])
    expect_equal(query(mosub), query(mo))
    expect_equal(target(mosub), target(mo))
    ## in range indexes
    idxs <- c(1, 3, 5)
    mosub <- filterMatches(mo, index = idxs)
    expect_equal(mosub@matches, mo@matches[idxs, ])
    expect_equal(query(mosub), query(mo))
    expect_equal(target(mosub), target(mo))
    mosub <- filterMatches(mo, index = idxs, keep = FALSE)
    expect_equal(mosub@matches, mo@matches[-idxs, ])
    ## keep matches based on query and target input values
    queryValue <- c(q1[mo@matches[idxs, "query_idx"], "col1"], -1, - 2)
    targetValue <- c(t1[mo@matches[idxs, "target_idx"], "col2"], -2, -1)
    mosub <- filterMatches(mo, queryValue = queryValue,
                         targetValue = targetValue,
                         queryColname = "col1", targetColname = "col2")
    expect_equal(mosub@matches, mo@matches[idxs, ])
    expect_equal(query(mosub), query(mo))
    expect_equal(target(mosub), target(mo))
    mosub <- filterMatches(mo, queryValue = queryValue,
                           targetValue = targetValue, queryColname = "col1",
                           targetColname = "col2", keep = FALSE)
    expect_equal(mosub@matches, mo@matches[-idxs, ])
    ## no matches corresponding to the input values
    mosub <- filterMatches(mo, queryValue = queryValue + 100,
                         targetValue = targetValue + 100,
                         queryColname = "col1", targetColname = "col2")
    expect_equal(mosub@matches, data.frame(query_idx = integer(),
                                           target_idx = integer(),
                                           score = numeric()))

    #### query and target vectors (multiple matches corresponding to a given
    #### couple of query and target input values)
    mo <- Matched(
        c(1, 1, 2, 3, 4), c("A", "A", "A", "C", "D"),
        matches = data.frame(query_idx = c(1L, 2L, 4L, 5L),
                             target_idx = c(1L, 2L, 4L, 5L),
                             score = seq(0.5, 0.8, by = 0.1)))
    queryValue <- c(1, 3)
    targetValue <- c("A", "C")
    mosub <- filterMatches(mo, queryValue = queryValue,
                         targetValue = targetValue)
    expect_equal(mosub@matches, mo@matches[c(1, 2, 3), ])

    #### query SummarizedExperiment and target data.frame
    mo <- Matched(
        q3, t3, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                     target_idx = c(2L, 2L, 3L, 4L, 5L),
                                     score = seq(0.5, 0.9, by = 0.1)))
    ## out of bounds indexes
    expect_error(filterMatches(mo, index = c(1, 10)), "out-of-bounds")
    ## no index : every match is removed
    idxs <- integer(0)
    mosub <- filterMatches(mo, index = idxs)
    expect_equal(mosub@matches, mo@matches[idxs, ])
    expect_equal(query(mosub), query(mo))
    expect_equal(target(mosub), target(mo))
    ## in range indexes
    idxs <- c(1, 3, 5)
    mosub <- filterMatches(mo, index = idxs)
    expect_equal(mosub@matches, mo@matches[idxs, ])
    expect_equal(query(mosub), query(mo))
    expect_equal(target(mosub), target(mo))
    ## keep matches based on query and target input values
    queryValue <- c(rowData(q3)[mo@matches[idxs, "query_idx"], "col1"], -1, - 2)
    targetValue <- c(t3[mo@matches[idxs, "target_idx"], "col2"], -2, -1)
    mosub <- filterMatches(mo, queryValue = queryValue,
                           targetValue = targetValue,
                           queryColname = "col1", targetColname = "col2")
    expect_equal(mosub@matches, mo@matches[idxs, ])
    expect_equal(query(mosub), query(mo))
    expect_equal(target(mosub), target(mo))
    ## no matches corresponding to the input values
    mosub <- filterMatches(mo, queryValue = queryValue + 100,
                           targetValue = targetValue + 100,
                           queryColname = "col1", targetColname = "col2")
    expect_equal(mosub@matches, data.frame(query_idx = integer(),
                                           target_idx = integer(),
                                           score = numeric()))

})

test_that("filterMatches,Matched,SelectMatchesParam works", {
    #### query and target data.frames
    mo <- Matched(
        q1, t1, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                     target_idx = c(2L, 2L, 3L, 4L, 5L),
                                     score = seq(0.5, 0.9, by = 0.1)))
    ## out of bounds indexes
    expect_error(filterMatches(mo, SelectMatchesParam(index = c(1L, 10L))),
                 "out-of-bounds")
    ## no index : every match is removed
    idxs <- integer(0)
    mosub <- filterMatches(mo, SelectMatchesParam(index = idxs))
    expect_equal(mosub@matches, mo@matches[idxs, ])
    expect_equal(query(mosub), query(mo))
    expect_equal(target(mosub), target(mo))
    ## in range indexes
    idxs <- c(1L, 3L, 5L)
    mosub <- filterMatches(mo, SelectMatchesParam(index = idxs))
    expect_equal(mosub@matches, mo@matches[idxs, ])
    expect_equal(query(mosub), query(mo))
    expect_equal(target(mosub), target(mo))
    mosub_old <- filterMatches(mo, index = idxs)
    expect_equal(mosub@matches, mosub_old@matches)
    mosub <- filterMatches(mo, SelectMatchesParam(index = idxs, keep = FALSE))
    expect_equal(mosub@matches, mo@matches[-idxs, ])
    ## keep matches based on query and target input values
    queryValue <- c(q1[mo@matches[idxs, "query_idx"], "col1"], -1, - 2)
    targetValue <- c(t1[mo@matches[idxs, "target_idx"], "col2"], -2, -1)
    mosub <- filterMatches(mo, SelectMatchesParam(queryValue = queryValue,
                                                  targetValue = targetValue,
                                                  queryColname = "col1",
                                                  targetColname = "col2"))
    expect_equal(mosub@matches, mo@matches[idxs, ])
    expect_equal(query(mosub), query(mo))
    expect_equal(target(mosub), target(mo))
    mosub <- filterMatches(mo, SelectMatchesParam(queryValue = queryValue,
                                                  targetValue = targetValue,
                                                  queryColname = "col1",
                                                  targetColname = "col2",
                                                  keep = FALSE))
    expect_equal(mosub@matches, mo@matches[-idxs, ])
    ## no matches corresponding to the input values
    mosub <- filterMatches(mo,
                           SelectMatchesParam(queryValue = queryValue + 100,
                                              targetValue = targetValue + 100,
                                              queryColname = "col1",
                                              targetColname = "col2"))
    expect_equal(mosub@matches, data.frame(query_idx = integer(),
                                           target_idx = integer(),
                                           score = numeric()))

    #### query and target vectors (multiple matches corresponding to a given
    #### couple of query and target input values)
    mo <- Matched(
        c(1, 1, 2, 3, 4), c("A", "A", "A", "C", "D"),
        matches = data.frame(query_idx = c(1L, 2L, 4L, 5L),
                             target_idx = c(1L, 2L, 4L, 5L),
                             score = seq(0.5, 0.8, by = 0.1)))
    queryValue <- c(1, 3)
    targetValue <- c("A", "C")
    mosub <- filterMatches(mo, SelectMatchesParam(queryValue = queryValue,
                                                  targetValue = targetValue))
    expect_equal(mosub@matches, mo@matches[c(1, 2, 3), ])

    #### query SummarizedExperiment and target data.frame
    mo <- Matched(
        q3, t3, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                     target_idx = c(2L, 2L, 3L, 4L, 5L),
                                     score = seq(0.5, 0.9, by = 0.1)))
    ## out of bounds indexes
    expect_error(filterMatches(mo, SelectMatchesParam(index = c(1L, 10L))),
                 "out-of-bounds")
    ## no index : every match is removed
    idxs <- integer(0)
    mosub <- filterMatches(mo, SelectMatchesParam(index = idxs))
    expect_equal(mosub@matches, mo@matches[idxs, ])
    expect_equal(query(mosub), query(mo))
    expect_equal(target(mosub), target(mo))
    ## in range indexes
    idxs <- c(1L, 3L, 5L)
    mosub <- filterMatches(mo, SelectMatchesParam(index = idxs))
    expect_equal(mosub@matches, mo@matches[idxs, ])
    expect_equal(query(mosub), query(mo))
    expect_equal(target(mosub), target(mo))
    ## keep matches based on query and target input values
    queryValue <- c(rowData(q3)[mo@matches[idxs, "query_idx"], "col1"], -1, - 2)
    targetValue <- c(t3[mo@matches[idxs, "target_idx"], "col2"], -2, -1)
    mosub <- filterMatches(mo, SelectMatchesParam(queryValue = queryValue,
                                                  targetValue = targetValue,
                                                  queryColname = "col1",
                                                  targetColname = "col2"))
    expect_equal(mosub@matches, mo@matches[idxs, ])
    expect_equal(query(mosub), query(mo))
    expect_equal(target(mosub), target(mo))
    ## no matches corresponding to the input values
    mosub <- filterMatches(mo,
                           SelectMatchesParam(queryValue = queryValue + 100,
                                              targetValue = targetValue + 100,
                                              queryColname = "col1",
                                              targetColname = "col2"))
    expect_equal(mosub@matches, data.frame(query_idx = integer(),
                                           target_idx = integer(),
                                           score = numeric()))
})

test_that("filterMatches,Matched,TopRankedMatchesParam works", {
    mo <- Matched(
        q1, t1, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                     target_idx = c(2L, 2L, 3L, 4L, 5L),
                                     score = c(4, 4, 1, 3, 1)))
    mosub <- filterMatches(mo, TopRankedMatchesParam())
    expect_equal(mosub@matches, mo@matches[c(1L, 3L, 5L), ])
    expect_equal(query(mosub), query(mo))
    expect_equal(target(mosub), target(mo))

    mosub <- filterMatches(mo, TopRankedMatchesParam(n = 2L))
    expect_equal(mosub@matches, mo@matches[c(1L, 3L, 4L, 5L), ])
    expect_equal(query(mosub), query(mo))
    expect_equal(target(mosub), target(mo))

    mosub <- filterMatches(mo, TopRankedMatchesParam(n = 10L))
    expect_equal(mosub@matches, mo@matches)
    expect_equal(query(mosub), query(mo))
    expect_equal(target(mosub), target(mo))

    mo <- Matched(
        q1, t1, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                     target_idx = c(2L, 2L, 3L, 4L, 5L),
                                     score = c(4, 4, 1, 3, 1),
                                     score_rt = c(1, -1, -1, -2, 1)))

    mosub <- filterMatches(mo, TopRankedMatchesParam())
    expect_equal(mosub@matches, mo@matches[c(1L, 3L, 5L), ])
    expect_equal(query(mosub), query(mo))
    expect_equal(target(mosub), target(mo))

    mosub <- filterMatches(mo, TopRankedMatchesParam(n = 2L))
    expect_equal(mosub@matches, mo@matches[c(1L, 2L, 3L, 5L), ])
    expect_equal(query(mosub), query(mo))
    expect_equal(target(mosub), target(mo))

    mosub <- filterMatches(mo, TopRankedMatchesParam(n = 10L))
    expect_equal(mosub@matches, mo@matches)
    expect_equal(query(mosub), query(mo))
    expect_equal(target(mosub), target(mo))

    q1 <- data.frame(col1 = 1:5, col2 = 6:10)
    t2 <- 11:16
    mo <- Matched(q1, t2,
                  matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                       target_idx = c(2L, 2L, 3L, 4L, 5L),
                                       score = seq(0.5, 0.9, by = 0.1)))
    res <- filterMatches(mo, TopRankedMatchesParam(n = 1L))
    expect_equal(res@matches$query_idx, c(1, 2, 5))
    expect_equal(res@matches$target_idx, c(2, 2, 5))
    expect_equal(res@matches$score, c(0.5, 0.6, 0.9))

    res <- filterMatches(mo, TopRankedMatchesParam(n = 1L, decreasing = TRUE))
    expect_equal(res@matches$query_idx, c(1, 2, 5))
    expect_equal(res@matches$target_idx, c(2, 4, 5))
    expect_equal(res@matches$score, c(0.5, 0.8, 0.9))
})

test_that("addMatches,Matched works", {
    #### query and target data.frames
    mo <- Matched(
        q1, t1, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                     target_idx = c(2L, 2L, 3L, 4L, 5L),
                                     score = seq(0.5, 0.9, by = 0.1)))
    moadd <- addMatches(mo, queryValue = q1[c(3, 5), "col1"],
                        targetValue = t1[c(1, 6), "col2"],
                        queryColname = "col1", targetColname = "col2",
                        score = data.frame(score = c(1, 1.1)), isIndex = FALSE)
    expect_equal(moadd@matches$query_idx, c(mo@matches$query_idx, 3, 5))
    expect_equal(moadd@matches$target_idx, c(mo@matches$target_idx, 1, 6))
    expect_equal(moadd@matches$score, c(mo@matches$score, 1, 1.1))
    expect_equal(query(moadd), query(mo))
    expect_equal(target(moadd), target(mo))

    moadd <- addMatches(mo, queryValue = c(3L, 5L),
                        targetValue = c(1L, 6L),
                        queryColname = "col1", targetColname = "col2",
                        score = data.frame(score = c(1, 1.1)), isIndex = TRUE)
    expect_equal(moadd@matches$query_idx, c(mo@matches$query_idx, 3, 5))
    expect_equal(moadd@matches$target_idx, c(mo@matches$target_idx, 1, 6))
    expect_equal(moadd@matches$score, c(mo@matches$score, 1, 1.1))
    expect_equal(query(moadd), query(mo))
    expect_equal(target(moadd), target(mo))

    #### query vector and target data.frame
    mo <- Matched(q1[, "col1"], t1,
                  matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                       target_idx = c(2L, 2L, 3L, 4L, 5L),
                                       score = seq(0.5, 0.9, by = 0.1)))
    moadd <- addMatches(mo, queryValue = c(q1[c(3, 5), "col1"], 50),
                        targetValue = c(t1[c(1, 6), "col2"], 50),
                        targetColname = "col2",
                        score = data.frame(score = c(1, 1.1, 2)))
    expect_equal(moadd@matches$query_idx, c(mo@matches$query_idx, 3, 5))
    expect_equal(moadd@matches$target_idx, c(mo@matches$target_idx, 1, 6))
    expect_equal(moadd@matches$score, c(mo@matches$score, 1, 1.1))
    expect_equal(query(moadd), query(mo))
    expect_equal(target(moadd), target(mo))

    moadd <- addMatches(mo, queryValue = c(3L, 5L),
                        targetValue = c(1L, 6L), targetColname = "col2",
                        score = data.frame(score = c(1, 1.1)), isIndex = TRUE)
    expect_equal(moadd@matches$query_idx, c(mo@matches$query_idx, 3, 5))
    expect_equal(moadd@matches$target_idx, c(mo@matches$target_idx, 1, 6))
    expect_equal(moadd@matches$score, c(mo@matches$score, 1, 1.1))
    expect_equal(query(moadd), query(mo))
    expect_equal(target(moadd), target(mo))

    #### query SummarizedExperiment and target data.frame
    mo <- Matched(
        query = q3, target = t3,
        matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                             target_idx = c(2L, 2L, 3L, 4L, 5L),
                             score = seq(0.5, 0.9, by = 0.1)))

    moadd <- addMatches(mo, queryValue = 2L, targetValue = 1L, isIndex = TRUE,
                        score = data.frame(score = 100))
    expect_true(matches(moadd)$query_idx[6] == 2)
    expect_true(matches(moadd)$target_idx[6] == 1)
    expect_true(matches(moadd)$score[6] == 100)
})

test_that(".extract_elements works", {
    a <- data.frame(a = 1:5)
    expect_equal(.extract_elements(a, 3)[, 1], 3)
    expect_equal(.extract_elements(a, c(3, 6))[, 1], c(3, NA))
    expect_equal(.extract_elements(a, c(NA, 1))[, 1], c(NA, 1))

    a <- 1:5
    expect_equal(.extract_elements(a, 3), 3)
    expect_equal(.extract_elements(a, c(3, 6)), c(3, NA))
    expect_equal(.extract_elements(a, c(NA, 1)), c(NA, 1))

    a <- list(1, 2, 3, 4, 5)
    expect_equal(.extract_elements(a, 3), list(3))
    expect_equal(.extract_elements(a, c(3, 6)), list(3, NULL))
    expect_equal(.extract_elements(a, c(NA, 1)), list(NA, 1))
})

test_that(".cnt works", {
    t <- data.frame(a = 1:5, b = 2:6)
    expect_equal(.cnt(t), c("target_a", "target_b"))

    t <- 1:5
    expect_equal(.cnt(t), "target")

    t <- array(data = 1, dim = c(1, 1, 1))
    expect_error(.cnt(t), "unsupported")
})

test_that(".validate_assay works", {
    expect_identical(.validate_assay(q4, c("a1", "a2")),
                     "`queryAssay` must be `character(1)`")
    expect_identical(.validate_assay(q4, "a3", "target"),
                     "No assay in `target` with name \"a3\"")
    expect_null(.validate_assay(q4, "a1"))
})

test_that(".objectToMatch works", {
    expect_identical(.objectToMatch(q1), q1)
    expect_identical(.objectToMatch(q2), q2)
    expect_identical(.objectToMatch(q3), rowData(q3))
    expect_error(.objectToMatch(q4), "has to be provided")
    expect_error(.objectToMatch(q4, c("a1", "a2")), "must be `character")
    expect_error(.objectToMatch(q4, c("a3")), "with name")
    expect_identical(.objectToMatch(q4, c("a1")), rowData(q4[["a1"]]))
    expect_identical(.objectToMatch(q1, "assayname", "col2"), q1[, "col2"])
    expect_identical(.objectToMatch(q3, colnames = "col2"),
                     rowData(q3)[, "col2"])
    expect_identical(.objectToMatch(q4, "a1", "col2"),
                     rowData(q4[["a1"]])[, "col2"])
    expect_error(.objectToMatch(q1, "a1", c("col1", "col3")), "Missing column")
    expect_error(.objectToMatch(q3, "a1", c("col1", "col3")), "Missing column")
    expect_error(.objectToMatch(q4, "a1", c("col1", "col3")), "Missing column")
})

test_that(".subset_qt works", {
    i <- c(1, 3)
    expect_identical(.subset_qt(q1, i = i), q1[i, ])
    expect_identical(.subset_qt(q2, i = i), q2[i])
    expect_identical(.subset_qt(q3, i = i), (q3)[i, ])
    res <- .subset_qt(q4, "a1", i = i)
    expect_is(res, "QFeatures")
    expect_equal(res[["a1"]], q4[["a1"]][i, ])
})

test_that("SelectMatchesParam works", {
    res <- SelectMatchesParam()
    expect_true(is(res, "SelectMatchesParam"))

    expect_error(SelectMatchesParam(queryValue = "A", targetValue = 1:2),
                 "must have the same length")
    expect_error(SelectMatchesParam(queryColname = c("A", "b")),
                 "cannot be of length greater than 1")
    expect_error(SelectMatchesParam(targetColname = c("A", "b")),
                 "cannot be of length greater than 1")
    expect_error(SelectMatchesParam(index = c(1L, -2L)),
                 "must contain positive integers")
    expect_error(SelectMatchesParam(keep = rep(TRUE, 2)),
                 "must be a logical of length 1")
})

test_that("TopRankedMatchesParam works", {
    res <- TopRankedMatchesParam()
    expect_true(is(res, "TopRankedMatchesParam"))

    expect_error(TopRankedMatchesParam(n = c(2L, 3L)), "length 1")
    expect_error(TopRankedMatchesParam(n = -4L), "positive integer")
})
