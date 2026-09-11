# A prediction-error criterion over a structural term's own parameters is
# refused by name, and the certificate of a fit with no outer gradient says
# which case it is.

sim_filter_panel <- function(seed, m_grp, ti, a = 0.25, b = 0.55, beta = 0.5) {
  set.seed(seed)
  id <- factor(rep(seq_len(m_grp), each = ti))
  n <- length(id)
  x <- stats::rnorm(n)
  om <- stats::rnorm(m_grp, 0.2, 0.35)
  y <- numeric(n)
  for (g in seq_len(m_grp)) {
    rows <- which(as.integer(id) == g)
    f <- om[g] / (1 - b)
    s <- 0
    for (k in seq_along(rows)) {
      f <- om[g] + a * s + b * f
      eta <- beta * x[rows[k]] + f
      y[rows[k]] <- eta + stats::rnorm(1)
      s <- y[rows[k]] - eta
    }
  }
  data.frame(id = id, t = rep(seq_len(ti), m_grp), x = x, y = y)
}

test_that("a prediction-error criterion refuses a penalty over a filter's own parameters", {
  # Measured before the refusal existed: the exact gradient of aic() there was
  # exactly 0, so the search stopped at its first evaluation reporting
  # convergence, and the traces it priced were 2 and 1 where the traces over
  # the joint vector are 7.433 and 22.410.
  skip_on_cran()
  dp <- sim_filter_panel(7, 3L, 35L)
  g1 <- distributions7::gaussian1_distrib()
  form <- y ~ x + gas(p = 1, q = 1, omega ~ random(~1 | id), by = id, time = t)
  expect_error(statmod(form, g1, dp, outer_criterion = aic()),
               "aic\\(\\) cannot select the hyperparameter 'sigma'.*Use reml\\(\\) or ml\\(\\)")
  expect_error(statmod(form, g1, dp, outer_criterion = bic()),
               "bic\\(\\) cannot select the hyperparameter 'sigma'.*Use reml\\(\\) or ml\\(\\)")
  expect_error(statmod(form, g1, dp, outer_criterion = cv()),
               "cv\\(\\) cannot select the hyperparameter 'sigma'.*cross-validation sweeps a path")
  # the route the message names is a fit
  f <- statmod(form, g1, dp)
  expect_true(is.finite(f@criterion))
})

test_that("a lasso over a filter's own parameters is refused by name rather than by a crash", {
  # It used to die building the path, on "'from' must be a finite number". No
  # remedy is named: holding lambda does not fit either, which a measurement
  # of the published release says as well.
  skip_on_cran()
  dp <- sim_filter_panel(7, 3L, 35L)
  form <- y ~ x + gas(p = 1, q = 1, omega ~ 0 + lasso(~id), by = id, time = t)
  expect_error(statmod(form, distributions7::gaussian1_distrib(), dp),
               "bic\\(\\) cannot select the hyperparameter 'lambda'.*read a curvature the kink does not have")
})

test_that("the certificate says which case a fit without an outer gradient is", {
  set.seed(3)
  ds <- data.frame(x = runif(200))
  ds$y <- sin(6 * ds$x) + stats::rnorm(200, sd = 0.3)
  g1 <- distributions7::gaussian1_distrib()
  form <- y ~ s(x, bspline_smooth(k = 8))

  # a smooth penalty chosen by aic(): the reason used to be the kinked one, of
  # a model carrying no kink
  fa <- statmod(form, g1, ds, outer_criterion = aic())
  ca <- statmod_certificate(fa)
  expect_equal(ca$state, "unknown")
  expect_match(ca$reason, "chosen by aic\\(\\), a prediction-error criterion")
  expect_false(grepl("kink", ca$reason))

  fb <- statmod(form, g1, ds, outer_criterion = bic())
  expect_match(statmod_certificate(fb)$reason, "chosen by bic\\(\\)")

  fh <- statmod(y ~ s(x, bspline_smooth(k = 8), hyper = c(lambda = 2)), g1, ds)
  expect_match(statmod_certificate(fh)$reason, "every hyperparameter here was held")

  fn <- statmod(form, g1, ds, outer_criterion = NULL)
  expect_match(statmod_certificate(fn)$reason, "no criterion was asked for")
  expect_false(grepl("kink", statmod_certificate(fn)$reason))
})
