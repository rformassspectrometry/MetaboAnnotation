test_that("matchFormula works", {

  # input formula
  x <- c("H12C6O6", "C11H12O2", "HN3")
  y <- c("HCl", "C2H4O", "C6H12O6")

  x_df <- data.frame(formula = c("H12C6O6", "C11H12O2", "HN3"),
                         name = c("A", "B", "C"))

  y_df <- data.frame(formula = c("HCl", "C2H4O", "C6H12O6"),
                          name = c("D", "E", "F"))

  # character vs character
  res <- matchFormula(x, y)
  expect_equal(query(res), x)
  expect_equal(target(res), y)
  expect_equal(res@matches$query_idx, c(1))
  expect_equal(res@matches$target_idx, c(3))
  expect_equal(res@matches$score, c(1))

  # data.frame vs data.frame
  res <- matchFormula(x_df, y_df)
  expect_equal(query(res), x_df)
  expect_equal(target(res), y_df)
  expect_equal(res@matches$query_idx, c(1))
  expect_equal(res@matches$target_idx, c(3))
  expect_equal(res@matches$score, c(1))

  # data.frame vs character
  res <- matchFormula(x_df, y)
  expect_equal(query(res), x_df)
  expect_equal(target(res), y)
  expect_equal(res@matches$query_idx, c(1))
  expect_equal(res@matches$target_idx, c(3))
  expect_equal(res@matches$score, c(1))

  # character vs data.frame
  res <- matchFormula(x, y_df)
  expect_equal(query(res), x)
  expect_equal(target(res), y_df)
  expect_equal(res@matches$query_idx, c(1))
  expect_equal(res@matches$target_idx, c(3))
  expect_equal(res@matches$score, c(1))
})
