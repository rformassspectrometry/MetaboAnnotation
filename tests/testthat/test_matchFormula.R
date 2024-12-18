test_that("matchFormula works", {
  # input formula
  x <- c("H12C6O6", "C11H12O2", "HN3", "ClH", "HCl")
  y <- c("HCl", "C2H4O", "C6H12O6", "N3H")

  # character vs character
  res <- matchFormula(x, y)
  expect_s4_class(res, "Matched")
  expect_equal(query(res), x)
  expect_equal(target(res), y)
  expect_equal(res@matches$query_idx, c(1L, 3L, 4L, 5L))
  expect_equal(res@matches$target_idx, c(3L, 4L, 1L, 1L))
  expect_equal(res@matches$score, rep(1L, 4))

  x_df <- data.frame(chem_formula = x,
                     name = c("A", "B", "C", "D", "E"))

  y_df <- data.frame(formula = y,
                     name = c("D", "E", "F", "G"))

  ## data.frame vs data.frame
  expect_error(matchFormula(x_df, y_df), "Missing column")
  expect_error(matchFormula(x_df, y_df, formulaColname = c("chem_formula")),
               "Missing column")
  res <- matchFormula(x_df, y_df, formulaColname = c("chem_formula", "formula"))
  expect_equal(query(res), x_df)
  expect_equal(target(res), y_df)
  expect_equal(res@matches$query_idx, c(1L, 3L, 4L, 5L))
  expect_equal(res@matches$target_idx, c(3L, 4L, 1L, 1L))
  expect_equal(res@matches$score, rep(1L, 4))

  # data.frame vs character
  res <- matchFormula(x_df, y, formulaColname = "chem_formula")
  expect_equal(query(res), x_df)
  expect_equal(target(res), y)
  expect_equal(res@matches$query_idx, c(1L, 3L, 4L, 5L))
  expect_equal(res@matches$target_idx, c(3L, 4L, 1L, 1L))
  expect_equal(res@matches$score, rep(1L, 4))
  expect_error(matchFormula(x_df, y, formulaColname = "other"),
               "Missing column")

  # character vs data.frame
  res <- matchFormula(x, y_df)
  expect_equal(query(res), x)
  expect_equal(target(res), y_df)
  expect_equal(res@matches$query_idx, c(1L, 3L, 4L, 5L))
  expect_equal(res@matches$target_idx, c(3L, 4L, 1L, 1L))
  expect_equal(res@matches$score, rep(1L, 4))
  expect_error(matchFormula(x, y_df, formulaColname = "other"),
               "Missing column")
})

test_that(".getFormulaMatches works", {
    res <- MetaboAnnotation:::.getFormulaMatches(3L, "A", c("A", "B", "A"))
    expect_true(is.data.frame(res))
    expect_equal(colnames(res), c("query_idx", "target_idx", "score"))
    expect_equal(res$query_idx, c(3L, 3L))
    expect_equal(res$target_idx, c(1L, 3L))
})
