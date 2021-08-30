library(SummarizedExperiment)
q1 <- SummarizedExperiment(
  assays = data.frame(matrix(NA, 5, 2)), 
  rowData = data.frame(col1 = 1:5, col2 = 6:10),
  colData = data.frame(cD1 = c(NA, NA), cD2 = c(NA, NA)))
t1 <- data.frame(col1 = 11:16, col2 = 17:22)

m <- data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                target_idx = c(2L, 2L, 3L, 4L, 5L),
                score = seq(0.5, 0.9, by = 0.1))

test_that("MatchedSummarizedExperiment works", {
  mo <- MatchedSummarizedExperiment()
  expect_true(validObject(mo))
  expect_output(show(mo), "MatchedSummarizedExperiment")
  expect_equal(whichQuery(mo), integer())
  expect_equal(whichTarget(mo), integer())
  
  mo <- MatchedSummarizedExperiment(query = q1, target = t1,
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
  mo <- MatchedSummarizedExperiment()
  expect_error(.subset_matches_nodim(mo, 1), "out-of-bounds")
  
  mo <- MatchedSummarizedExperiment(
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
  mo <- MatchedSummarizedExperiment(
    q1, t1, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                 target_idx = c(2L, 2L, 3L, 4L, 5L),
                                 score = seq(0.5, 0.9, by = 0.1)))
  
  res <- mo[whichQuery(mo)]
  expect_equal(query(res), q1[whichQuery(mo), ])
  res <- pruneTarget(res)
  expect_equal(target(res), t1[c(2, 3, 4, 5), ])
  expect_equal(res@matches$target_idx, c(1L, 1L, 2L, 3L, 4L))
})

test_that("$,MatchedSummarizedExperiment works", {
  mo <- MatchedSummarizedExperiment(
    q1, t1, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                 target_idx = c(2L, 2L, 3L, 4L, 5L),
                                 score = seq(0.5, 0.9, by = 0.1)))
  expect_equal(mo$score, c(0.5, 0.6, 0.7, 0.8, NA, NA, 0.9))
  expect_equal(mo$col1, rowData(q1)[c(1, 2, 2, 2, 3, 4, 5), "col1"])
  expect_equal(mo$target_col1, t1[c(2,2,3,4, NA, NA, 5), "col1"])
})

test_that("colnames,MatchedSummarizedExperiment works", {
  mo <- MatchedSummarizedExperiment(
    q1, t1, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                 target_idx = c(2L, 2L, 3L, 4L, 5L),
                                 score = seq(0.5, 0.9, by = 0.1)))
  expect_equal(colnames(mo), c(colnames(rowData(q1)),
                               paste0("target_", colnames(t1)), "score"))
})

test_that("matchedData,MatchedSummarizedExperiment works", {
  
  mo <- MatchedSummarizedExperiment()
  res <- matchedData(mo)
  expect_true(is(res, "DataFrame"))
  expect_true(nrow(res) == 0)
  expect_equal(colnames(res), colnames(mo))
  
  mo <- MatchedSummarizedExperiment(
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
  
  ## A MatchedSummarizedExperiment with no matching target elements
  mo <- MatchedSummarizedExperiment(q1, t1, matches = data.frame(query_idx = integer(),
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

test_that("pruneTarget,MatchedSummarizedExperiment works", {
  mo <- MatchedSummarizedExperiment()
  res <- pruneTarget(mo)
  expect_true(is(res, "MatchedSummarizedExperiment"))
  expect_equal(res, mo)
  
  mo <- MatchedSummarizedExperiment(
    q1, t1, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                 target_idx = c(2L, 2L, 3L, 4L, 5L),
                                 score = seq(0.5, 0.9, by = 0.1)))
  res <- pruneTarget(mo)
  expect_equal(matchedData(res), matchedData(mo))
  expect_true(nrow(res@target) < nrow(mo@target))
  
  mo <- MatchedSummarizedExperiment(q1, t1, matches = data.frame(query_idx = integer(),
                                             target_idx = integer(),
                                             score = numeric()))
  res <- pruneTarget(mo)
  expect_equal(matchedData(res), matchedData(mo))
  expect_true(nrow(res@target) < nrow(mo@target))
})

test_that("filterMatches,Matched works", {
  #### query and target data.frames
  mo <- MatchedSummarizedExperiment(
    q1, t1, matches = data.frame(query_idx = c(1L, 2L, 2L, 2L, 5L),
                                 target_idx = c(2L, 2L, 3L, 4L, 5L),
                                 score = seq(0.5, 0.9, by = 0.1)))
  ## out of bounds indexes
  expect_error(filterMatches(mo, idxs = c(1, 10)), "out of bounds")
  ## no index : every match is removed
  idxs <- integer(0)
  mosub <- filterMatches(mo, idxs = idxs)
  expect_equal(mosub@matches, mo@matches[idxs, ])
  expect_equal(query(mosub), query(mo))
  expect_equal(target(mosub), target(mo))
  ## in range indexes
  idxs <- c(1, 3, 5)
  mosub <- filterMatches(mo, idxs = idxs)
  expect_equal(mosub@matches, mo@matches[idxs, ])
  expect_equal(query(mosub), query(mo))
  expect_equal(target(mosub), target(mo))
  ## keep matches based on query and target input values
  queryValues <- c(rowData(q1)[mo@matches[idxs, "query_idx"], "col1"], -1, - 2)
  targetValues <- c(t1[mo@matches[idxs, "target_idx"], "col2"], -2, -1)
  mosub <- filterMatches(mo, queryValues = queryValues, 
                       targetValues = targetValues, 
                       queryColname = "col1", targetColname = "col2")
  expect_equal(mosub@matches, mo@matches[idxs, ])
  expect_equal(query(mosub), query(mo))
  expect_equal(target(mosub), target(mo))
  ## no matches corresponding to the input values
  mosub <- filterMatches(mo, queryValues = queryValues + 100, 
                       targetValues = targetValues + 100, 
                       queryColname = "col1", targetColname = "col2")
  expect_equal(mosub@matches, data.frame(query_idx = integer(), 
                                         target_idx = integer(), 
                                         score = numeric()))
})