skip_on_cran()

test_that("a lasso over two identical columns chooses what it chooses without one", {
  # the effect can be shared between the two copies, the information over
  # them is singular, and the points of the path where both were free read
  # as unsettled and were dropped: the path chose the empty fit, lambda 41.8
  # on an effect of 2
  set.seed(2)
  x1 <- rnorm(100)
  d <- data.frame(x1 = x1, x2 = x1, x3 = rnorm(100))
  d$y <- 1 + 2 * x1 + rnorm(100)
  f <- suppressWarnings(statmod(y ~ lasso(~ x1 + x2 + x3),
                                distrib = distributions7::gaussian1_distrib(),
                                data = d))
  f0 <- statmod(y ~ lasso(~ x1 + x3),
                distrib = distributions7::gaussian1_distrib(), data = d)
  expect_equal(hyper(f)$estimate, hyper(f0)$estimate, tolerance = 1e-8)
  cf <- coef(f)$mu
  # the two copies carry the effect of one between them
  expect_equal(sum(f@coefficients$mu[2:3]), coef(f0)$mu[["lasso.x1"]],
               tolerance = 1e-6)
  expect_gt(cf[["lasso.x1"]], 1)
  # the split between the copies is not identified, and the summary says so
  expect_identical(f@aliased, "mu:lasso.x2")
  out <- capture.output(print(summary(f)))
  expect_true(any(grepl("1 not identified", out, fixed = TRUE)))
})
