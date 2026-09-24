# A zero-inflated model's mixing weight starts at the proportion of zeros,
# on its own link.
#
# The intercept-only fit puts the weight at the edge of its domain wherever the
# parent's overdispersion absorbs the excess zeros, and a model with
# covariates started there never leaves: the link is flat, the score on the
# unconstrained scale vanishes, and the fit reports convergence at a weight of
# zero. The sample below is the one that showed it: 7.6 log-likelihood units
# below the interior maximum.

zi_sample <- function() {
  n <- 1500
  set.seed(20260924)
  x <- runif(n, -1, 1)
  mu <- exp(0.6 + 0.9 * x)
  set.seed(5)
  y <- rnbinom(n, size = 2, mu = mu)
  y[runif(n) < 0.25] <- 0
  data.frame(y = y, x = x)
}

test_that("the start is the link of the proportion of zeros, whatever the link", {
  dat <- zi_sample()
  p0 <- mean(dat$y == 0)
  for (lk in list(linkfunctions7::logit_link(), linkfunctions7::probit_link(),
                  linkfunctions7::cloglog_link())) {
    d <- distributions7::zero_inflated(distributions7::negbin2_distrib(),
                                       link_zi = lk)
    spec <- statmod_spec(y ~ x, d, dat)
    expect_equal(statmod_intercepts(spec)$zi, linkfunctions7::linkfun(lk, p0),
                 label = lk@link_name)
  }
  # and nothing else moves on a family that asks for nothing
  spec <- statmod_spec(y ~ x, distributions7::negbin2_distrib(), dat)
  expect_null(statmod_intercepts(spec)$zi)
})

test_that("from there the fit reaches the interior maximum", {
  skip_on_cran()
  dat <- zi_sample()
  d <- distributions7::zero_inflated(distributions7::negbin2_distrib())
  fit <- statmod(y ~ x, d, dat)
  expect_true(fit@converged)
  expect_gt(plogis(fit@coefficients$zi), 0.1)
  # the profile peaks between 0.2 and 0.25 at -2396.02; the edge sits at
  # -2403.62, where the intercept-only start used to leave the fit
  expect_gt(as.numeric(logLik(fit)), -2397)
})
