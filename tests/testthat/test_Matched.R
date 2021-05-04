q1 <- data.frame(col1 = 1:5, col2 = 6:10)
t1 <- data.frame(col1 = 11:16, col2 = 17:22)
q2 <- list("A", c(1,2), 1L, list(), data.frame())
t2 <- list("A", c(1,2), 1L, list(), data.frame(), matrix(1:3))

m <- data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                target_idx = c(2L, 2L, 3L, 4L, 5L),
                score = seq(0.5, 0.9, by = 0.1))

#### maybe to remove in test_MatchedSpectra.R?
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
    
    ## query and target: lists
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
})

#### maybe to remove in test_MatchedSpectra.R?
test_that(".fill_index works", {
    res <- .fill_index(1:5, integer())
    expect_equal(res, 1:5)
    
    res <- .fill_index(1:20, c(2, 4, 4, 4, 6, 9, 10))
    expect_equal(res, c(1, 2, 3, 4, 4, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
                        14, 15, 16, 17, 18, 19, 20))
})

test_that("$,Matched works", {
    mo <- Matched(
        q1, t1, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                     target_idx = c(2L, 2L, 3L, 4L, 5L),
                                     score = seq(0.5, 0.9, by = 0.1)))
    expect_equal(mo$score, c(0.5, 0.6, 0.7, 0.8, NA, NA, 0.9))
    expect_equal(mo$col1, q1[c(1, 2, 2, 2, 3, 4, 5), "col1"])
    expect_equal(mo$target_col1, t1[c(2,2,3,4, NA, NA, 5), "col1"])
    
    mo <- Matched(
        q2, t2, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                     target_idx = c(2L, 2L, 3L, 4L, 5L),
                                     score = seq(0.5, 0.9, by = 0.1)))
    expect_equal(mo$score, c(0.5, 0.6, 0.7, 0.8, NA, NA, 0.9))
    expect_equal(mo$query, q2[c(1, 2, 2, 2, 3, 4, 5)])
    expect_equal(mo$target, t2[c(2,2,3,4, NA, NA, 5)])
    
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
    expect_equal(mo$target, t2[c(2,2,3,4, NA, NA, 5)])
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
    
})

library(S4Vectors)

test_that("data,Matched works", {
    
    mo <- Matched()
    res <- data(mo)
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) == 0)
    expect_equal(colnames(res), colnames(mo))
    
    #### query and target data.frames 
    mo <- Matched(
        q1, t1, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                     target_idx = c(2L, 2L, 3L, 4L, 5L),
                                     score = seq(0.5, 0.9, by = 0.1)))
    res <- data(mo)
    expect_equal(res$col2, mo$col2)
    expect_equal(res$target_col2, mo$target_col2)
    expect_equal(res$score, mo$score)
    
    expect_error(data(mo, columns = "other"), "other not available")
    
    ## Only query object variables
    res <- data(mo, columns = c("col1", "col2"))
    expect_equal(colnames(res), c("col1", "col2"))
    expect_equal(res$col1, mo$col1)

    ## Only target object variables
    res <- data(mo, columns = c("target_col1", "target_col2"))
    expect_equal(colnames(res), c("target_col1", "target_col2"))
    expect_equal(res$target_col1, mo$target_col1)
    
    ## Only matches
    res <- data(mo, columns = c("score"))
    expect_equal(colnames(res), c("score"))
    expect_equal(res$score, mo$score)
    
    #### A Matched with no matching target elements (query and target data.frames)
    mo <- Matched(q1, t1, matches = data.frame(query_idx = integer(),
                                               target_idx = integer(),
                                               score = numeric()))
    res <- data(mo)
    expect_equal(res$col2, mo$col2)
    expect_equal(res$target_col2, mo$target_col2)
    expect_equal(res$score, mo$score)
    
    ## Only query object variables
    res <- data(mo, columns = c("col1", "col2"))
    expect_equal(colnames(res), c("col1", "col2"))
    expect_equal(res$col1, mo$col1)
    
    ## Only target object variables
    res <- data(mo, columns = c("target_col1", "target_col2"))
    expect_equal(colnames(res), c("target_col1", "target_col2"))
    expect_equal(res$target_col1, mo$target_col1)
    
    ## Only matches
    res <- data(mo, columns = c("score"))
    expect_equal(colnames(res), c("score"))
    expect_equal(res$score, mo$score)
    
    #### query data.frame and target list 
    mo <- Matched(
        q1, t2, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                     target_idx = c(2L, 2L, 3L, 4L, 5L),
                                     score = seq(0.5, 0.9, by = 0.1)))
    res <- data(mo)
    expect_equal(res$col2, mo$col2)
    expect_equal(res$target, I(mo$target))
    expect_equal(res$score, mo$score)
    
    ## Only query object variables
    res <- data(mo, columns = c("col1", "col2"))
    expect_equal(colnames(res), c("col1", "col2"))
    expect_equal(res$col1, mo$col1)
    
    ## Only target object variables
    res <- data(mo, columns = c("target"))
    expect_equal(colnames(res), c("target"))
    expect_equal(res$target, mo$target)
    
    ## Only matches
    res <- data(mo, columns = c("score"))
    expect_equal(colnames(res), c("score"))
    expect_equal(res$score, mo$score)
    
    #### A Matched with no matching target elements: query list and target 
    #### data.frame
    mo <- Matched(q2, t1, matches = data.frame(query_idx = integer(),
                                               target_idx = integer(),
                                               score = numeric()))
    res <- data(mo)
    expect_equal(res$query, I(mo$query))
    expect_equal(res$target_col1, mo$target_col1)
    expect_equal(res$score, mo$score)
    
    ## Only query object variables
    res <- data(mo, columns = c("query"))
    expect_equal(colnames(res), c("query"))
    expect_equal(res$col1, mo$col1)
    
    ## Only target object variables
    res <- data(mo, columns = c("target_col1", "target_col2"))
    expect_equal(colnames(res), c("target_col1", "target_col2"))
    expect_equal(res$target_col1, mo$target_col1)
    
    ## Only matches
    res <- data(mo, columns = c("score"))
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
    expect_equal(data(res), data(mo))
    expect_true(nrow(res@target) < nrow(mo@target))
    
    mo <- Matched(q1, t1, matches = data.frame(query_idx = integer(),
                                               target_idx = integer(),
                                               score = numeric()))
    res <- pruneTarget(mo)
    expect_equal(data(res), data(mo))
    expect_true(nrow(res@target) < nrow(mo@target))
    
    #### query data.frame and target list
    mo <- Matched(
        q1, t2, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                     target_idx = c(2L, 2L, 3L, 4L, 5L),
                                     score = seq(0.5, 0.9, by = 0.1)))
    res <- pruneTarget(mo)
    expect_equal(data(res), data(mo))
    expect_true(length(res@target) < length(mo@target))
    
    mo <- Matched(q1, t2, matches = data.frame(query_idx = integer(),
                                               target_idx = integer(),
                                               score = numeric()))
    res <- pruneTarget(mo)
    expect_equal(data(res), data(mo))
    expect_true(length(res@target) < length(mo@target))
    
    
})

