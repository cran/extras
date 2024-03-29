test_that("translations", {
  expect_identical(pow(5, 2), 25)
  expect_identical(phi(0), 0.5)
  expect_equal(phi(2), 0.9772499, tolerance = 0.0000001)
  x <- NA
  log(x) <- log(5)
  expect_equal(x, 5)
  expect_equal(logit(0.5), 0)
  expect_equal(logit(1), Inf)
  x <- NA
  logit(x) <- logit(0.75)
  expect_equal(x, 0.75)
  expect_equal(ilogit(logit(0.67)), 0.67)
  expect_equal(invlogit(logit(0.67)), 0.67)
  expect_equal(inv_logit(logit(0.67)), 0.67)
})

test_that("translations2", {
  x <- seq(0, 1, by = 0.25)
  expect_identical(logit(x), qlogis(x))
  expect_identical(ilogit(logit(x)), x)
  expect_identical(invlogit(logit(x)), x)
  expect_identical(inv_logit(logit(x)), x)

  logit(x) <- c(0.5, 1)
  expect_identical(x, ilogit(c(0.5, 1)))

  log(x) <- c(0.5, 1)
  expect_identical(x, exp(c(0.5, 1)))

  expect_identical(pow(3, 4), 3^4)
  expect_equal(phi(0:2), c(0.5, 0.8413447, 0.9772499), tolerance = 0.0000001)
})
