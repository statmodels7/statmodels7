test_that("the raw interval table reads each standard error by name beside a filter", {
  # The raw variance matrix carries a structural term's free parameters right
  # after the coefficients of the equation the term sits in. Read by
  # position, a filter in the mean gave the scale's intercept the standard
  # error of the filter's log-loading, and summary() printed it.
  set.seed(3)
  n <- 120L
  f <- numeric(n)
  y <- numeric(n)
  e <- stats::rnorm(n)
  for (t in seq_len(n)) {
    if (t > 1L) f[t] <- 0.4 * e[t - 1L] + 0.7 * f[t - 1L]
    y[t] <- 2 + f[t] + e[t]
  }
  dd <- data.frame(y = y, x = stats::rnorm(n))
  fit <- statmod(y ~ gas(1, 1) | sigma ~ x,
                 distributions7::gaussian1_distrib(), dd)
  V <- vcov(fit, readable = FALSE)
  ci <- confint(fit, readable = FALSE)
  keys <- c("mu:(Intercept)", "sigma:(Intercept)", "sigma:x")
  expect_equal(ci[keys, "se"], unname(sqrt(diag(V)[keys])), tolerance = 1e-12)
  # the readable table agrees on the coefficients both report
  expect_equal(ci[keys, "se"], unname(confint(fit)[keys, "se"]),
               tolerance = 1e-12)
  # negative control: the matrix's positions 2 and 3 are the filter's, so a
  # positional read would have put them on the scale's rows
  expect_false(isTRUE(all.equal(ci["sigma:(Intercept)", "se"],
                                unname(sqrt(diag(V))[2L]))))
})

test_that("a loading far up its log link is an estimate, not an edge", {
  # the score of a gaussian mean is in units of 1/variance, so on data in the
  # hundreds the loading is in the thousands and its free value passes 8;
  # a log link has a bound towards zero only, so the loading keeps its
  # standard error and interval
  set.seed(3)
  n <- 120L
  f <- numeric(n)
  e <- stats::rnorm(n, sd = 200)
  y <- numeric(n)
  for (t in seq_len(n)) {
    if (t > 1L) f[t] <- 0.4 * e[t - 1L] + 0.7 * f[t - 1L]
    y[t] <- 900 + f[t] + e[t]
  }
  fit <- statmod(y ~ gas(1, 1), distributions7::gaussian1_distrib(),
                 data.frame(y = y))
  z <- fit@structural[[1]]$unconstrained
  expect_gt(z[["alpha1"]], 8)
  ci <- confint(fit)
  expect_true(is.finite(ci["mu:gas.alpha1", "se"]))
  expect_true(ci["mu:gas.alpha1", "lower"] > 0)
})
