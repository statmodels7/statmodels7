# The least-squares start of a location on the identity link.

test_that("a location on the identity link starts from least squares", {
  set.seed(8)
  dat <- data.frame(x = runif(200) - 0.5)
  dat$y <- 1 + 6 * dat$x + rnorm(200, sd = 0.5)
  spec <- statmod_spec(y ~ x, distrib = distributions7::gaussian1_distrib(),
                       data = dat)
  design <- statmod_design(spec)
  st <- location_ols_start(spec, design)
  expect_false(is.null(st))
  ols <- stats::coef(stats::lm(y ~ x, data = dat))
  expect_equal(st$coef, unname(ols), tolerance = 1e-6)
  # sigma starts at the spread of the residuals, not of the response
  expect_equal(exp(st$eta0$sigma), sd(residuals(stats::lm(y ~ x, data = dat))),
               tolerance = 0.02)
  start <- start_at(start_intercepts(), spec, design, NULL)
  expect_equal(start$mu, unname(ols), tolerance = 1e-6)
})

test_that("the least-squares start leaves other equations alone", {
  set.seed(9)
  dat <- data.frame(x = runif(100))
  dat$y <- rpois(100, exp(1 + dat$x))
  # a log link: not a location on the identity link
  spec <- statmod_spec(y ~ x, distrib = distributions7::poisson_distrib(),
                       data = dat)
  expect_null(location_ols_start(spec, statmod_design(spec)))
  # an intercept-only location equation has nothing to fit
  dat$z <- rnorm(100)
  spec <- statmod_spec(z ~ 1, distrib = distributions7::gaussian1_distrib(),
                       data = dat)
  expect_null(location_ols_start(spec, statmod_design(spec)))
})

test_that("an offset is removed before the least-squares fit", {
  set.seed(10)
  dat <- data.frame(x = runif(150) - 0.5, o = rnorm(150))
  dat$y <- 2 + 3 * dat$x + dat$o + rnorm(150, sd = 0.3)
  spec <- statmod_spec(y ~ x + offset(o), distrib = distributions7::gaussian1_distrib(),
                       data = dat)
  st <- location_ols_start(spec, statmod_design(spec))
  expect_equal(st$coef, unname(stats::coef(stats::lm(I(y - o) ~ x, data = dat))),
               tolerance = 1e-6)
})

test_that("the direct skew normal reaches its maximum from the least-squares start", {
  set.seed(7)
  d <- distributions7::skewnormal1_distrib()
  dat <- data.frame(x = runif(300) - 0.5)
  dat$y <- distributions7::distrib_rng(d, 300, list(mu = 1 + 10 * dat$x,
                                                    sigma = 2, alpha = 3))
  f <- statmod(y ~ x, distrib = d, data = dat)
  ref <- statmod(y ~ x, distrib = d, data = dat,
                 start = list(mu = c(1, 10), sigma = log(2), alpha = 3))
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(ref)), tolerance = 1e-6)
})
