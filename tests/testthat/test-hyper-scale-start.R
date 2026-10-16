test_that("a matrix prior starts at half the working variance", {
  set.seed(31)
  m <- 15; k <- 12
  g <- factor(rep(seq_len(m), each = k)); x <- runif(m * k)
  y <- 1000 * (1 + 2 * x + rnorm(m, 0, 0.8)[g] + rnorm(m, 0, 0.5)[g] * x +
                 rnorm(m * k, 0, 0.6))
  dd <- data.frame(y = y, x = x, g = g)
  spec <- statmod_spec(y ~ x + random(~ x | g), distributions7::gaussian1_distrib(),
                       dd, NULL, NULL, linpar = linpar_options())
  des <- statmod_design(spec)
  h0 <- statmod_hyper_start(spec, des)
  h1 <- statmod_hyper_scale_start(spec, des, h0)
  th <- h1$mu[[1L]]
  # log_cholesky: the diagonal of the factor is the square root of the variance
  expect_equal(unname(th[c("sigma_log_L1", "sigma_log_L2")]),
               rep(0.5 * log(stats::var(y) * (length(y) - 1) / length(y) / 2), 2))
  expect_equal(unname(th[["sigma_L2.1"]]), 0)
  # a univariate prior is left at its neutral start
  sp2 <- statmod_spec(y ~ random(~ 1 | g), distributions7::gaussian1_distrib(),
                      dd, NULL, NULL, linpar = linpar_options())
  d2 <- statmod_design(sp2)
  expect_identical(statmod_hyper_scale_start(sp2, d2, statmod_hyper_start(sp2, d2)),
                   statmod_hyper_start(sp2, d2))
})

test_that("a correlated random effect is fitted equivariantly in the response scale", {
  set.seed(31)
  m <- 15; k <- 12
  g <- factor(rep(seq_len(m), each = k)); x <- runif(m * k)
  y <- 1 + 2 * x + rnorm(m, 0, 0.8)[g] + rnorm(m, 0, 0.5)[g] * x +
    rnorm(m * k, 0, 0.6)
  f1 <- statmod(y ~ x + random(~ x | g), data = data.frame(y = y, x = x, g = g),
                distrib = distributions7::gaussian1_distrib())
  f3 <- statmod(y ~ x + random(~ x | g),
                data = data.frame(y = 1000 * y, x = x, g = g),
                distrib = distributions7::gaussian1_distrib())
  # REML over n - p = 178 residual dimensions: the criterion moves by
  # -178 log(1000) and the factor of the covariance by 1000
  expect_equal(f3@criterion - f1@criterion, -178 * log(1000), tolerance = 1e-6)
  h1 <- unlist(f1@hyper); h3 <- unlist(f3@hyper)
  expect_equal(unname(h3[1:2] - h1[1:2]), rep(log(1000), 2), tolerance = 1e-4)
  expect_equal(unname(h3[3] / h1[3]), 1000, tolerance = 1e-3)
})
