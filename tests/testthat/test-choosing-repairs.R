# The repairs of statmodels7 0.166.0, found while measuring the claims of the
# book's chapter on choosing a model. Each block states the defect it pins.

test_that("prior weights naming a column of data are read there", {
  set.seed(1)
  d <- data.frame(x = runif(120))
  d$y <- 1 + d$x + rnorm(120)
  d$wt <- runif(120, 0.5, 2)
  wv <- d$wt
  a <- statmod(y ~ x, distrib = distributions7::gaussian1_distrib(), data = d,
               weights = wt)
  b <- statmod(y ~ x, distrib = distributions7::gaussian1_distrib(), data = d,
               weights = wv)
  expect_identical(a@coefficients, b@coefficients)
  expect_equal(a@coefficients$mu,
               unname(stats::coef(stats::lm(y ~ x, d, weights = wt))),
               tolerance = 1e-8)
})

test_that("a multivariate family is rejected by name", {
  set.seed(2)
  d <- data.frame(x = runif(50))
  d$Y <- cbind(rnorm(50), rnorm(50))
  expect_error(statmod(Y ~ x, distrib = distributions7::mvgaussian1_distrib(2),
                       data = d),
               "does not fit a multivariate response")
})

test_that("a Laplace response with no smooth hyperparameter fits by default", {
  set.seed(3)
  d <- data.frame(x = runif(150))
  d$y <- 1 + 2 * d$x + stats::rexp(150) * sample(c(-1, 1), 150, TRUE)
  f <- statmod(y ~ x, distrib = distributions7::laplace_distrib(), data = d)
  expect_null(f@methods$outer)
  expect_true(all(is.finite(f@coefficients$mu)))
  # a criterion named by the caller is still refused, by name
  expect_error(statmod(y ~ x, distrib = distributions7::laplace_distrib(),
                       data = d, outer_criterion = reml()),
               "needs 2 derivatives")
})

test_that("a count far in the upper tail has a finite quantile residual", {
  set.seed(4)
  d <- data.frame(x = rep(0, 60))
  d$y <- c(stats::rpois(58, 15), 67, 69)
  f <- statmod(y ~ 1, distrib = distributions7::poisson_distrib(), data = d)
  r <- residuals(f, seed = 7)
  expect_true(all(is.finite(r)))
  # where F <= 1/2 the residual is the lower-tail one, unchanged
  mu <- stats::predict(f, "mu", d)
  fy <- stats::ppois(d$y, mu)
  lo <- pmax(fy - stats::dpois(d$y, mu), 0)
  set.seed(7)
  v <- stats::runif(nrow(d))
  low <- fy <= 0.5
  expect_identical(r[low], stats::qnorm(lo + v * (fy - lo))[low])
  # and in the upper tail it is read from the survival function
  hi <- !low
  s <- stats::ppois(d$y, mu, lower.tail = FALSE)
  w <- s + (1 - v) * pmin(stats::dpois(d$y, mu), fy)
  expect_equal(r[hi], stats::qnorm(w, lower.tail = FALSE)[hi], tolerance = 1e-12)
  expect_gt(max(r), 8)
})

test_that("beside a kinked penalty the default marginal names nothing", {
  d <- uscrime
  d$lrate <- log(d$y)
  xs <- ~ M + So + Ed + Po1 + Po2 + LF + M.F + Pop + NW + U1 + U2 + GDP + Ineq +
    Prob + Time
  a <- statmod(lrate ~ scad(xs, standardize = TRUE),
               distrib = distributions7::gaussian1_distrib(), data = d)
  b <- statmod(lrate ~ scad(xs, standardize = TRUE),
               distrib = distributions7::gaussian1_distrib(), data = d,
               outer_criterion = reml(marginal = "none"))
  expect_identical(a@coefficients, b@coefficients)
  expect_true(a@converged)
  # the path scores every point: no point was lost to a search on sigma
  expect_false(anyNA(a@history$outer$criterion))
})

test_that("a kinked block's step is halved until the objective falls", {
  d <- uscrime
  d$lrate <- log(d$y)
  xs <- ~ M + So + Ed + Po1 + Po2 + LF + M.F + Pop + NW + U1 + U2 + GDP + Ineq +
    Prob + Time
  f <- statmod(lrate ~ Po1 + Ed + Ineq + Prob |
                 sigma ~ lasso(xs, standardize = TRUE, lambda = 0.4543283),
               distrib = distributions7::gaussian1_distrib(), data = d,
               outer_criterion = NULL)
  expect_true(f@converged)
  ob <- f@history$blocks$objective
  expect_true(all(diff(ob) <= 1e-8 * max(1, abs(ob))))
})

test_that("a Newton search that did not settle is run again by lbfgs", {
  skip_on_cran()
  form <- GAG ~ s(Age, bspline_smooth(k = 15)) |
    sigma ~ s(Age, bspline_smooth(k = 8)) |
    nu ~ s(Age, bspline_smooth(k = 6))
  # with an inner budget of 100 the Newton search stopped at -845.09 and
  # lbfgs, run again from the same start, reaches -816.59. Since 0.197.0 the
  # default inner route repairs an indefinite observed information and
  # settles within that budget, so the unsettled search is reproduced on
  # Fisher scoring named explicitly
  f <- statmod(form, distrib = distributions7::student_t1_distrib(),
               data = gagurine,
               inner_optimizer = iwls(maxit = 100, hessian = "expected"))
  expect_true(inherits(f@methods$search, "optimizers7::Lbfgs"))
  expect_gt(as.numeric(logLik(f, type = "marginal")), -817)
  expect_identical(statmod_certificate(f)$state, "boundary")
  # since 0.182.0 the default budget of 1000 lets Newton settle there itself
  g <- statmod(form, distrib = distributions7::student_t1_distrib(),
               data = gagurine)
  expect_gt(as.numeric(logLik(g, type = "marginal")), -817)
  expect_identical(statmod_certificate(g)$state, "boundary")
})
