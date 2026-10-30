# The starts of the score-driven terms: two by default, tried one term at a
# time.

gas_series <- function(seed, n = 800) {
  set.seed(seed)
  y <- numeric(n)
  f <- -0.15 / (1 - 0.95)
  for (t in seq_len(n)) {
    y[t] <- rnorm(1, 0, exp(f))
    f <- -0.15 + 0.95 * f + 0.15 * ((y[t] / exp(f))^2 - 1)
  }
  data.frame(t = seq_len(n), y = y)
}

test_that("the start at a high persistence reaches the maximum that the first misses", {
  d <- gas_series(1)
  fam <- distributions7::gaussian1_distrib()
  f1 <- statmod(y ~ 1 | sigma ~ gas(1, 1, time = t, n_start = 1), fam, d,
                outer_criterion = ml())
  f2 <- statmod(y ~ 1 | sigma ~ gas(1, 1, time = t), fam, d, outer_criterion = ml())
  expect_lt(coef(f1)$sigma[["gas.xi1"]], 0.1)
  expect_gt(coef(f2)$sigma[["gas.xi1"]], 0.9)
  expect_gt(as.numeric(logLik(f2)) - as.numeric(logLik(f1)), 100)
})

test_that("the joint filter carries the starts of its terms in turn", {
  set.seed(2)
  d <- data.frame(t = 1:120, y = rnorm(120))
  sp <- statmod_spec(y ~ gas(1, 1, time = t) | sigma ~ gas(1, 1, time = t, n_start = 3),
                     distributions7::gaussian1_distrib(), d)
  ms <- structural_multistarts(sp, statmod_design(sp))
  expect_length(ms, 1L)
  expect_identical(length(ms[[1L]]$zeta), 1L + 1L + 2L)
  sp1 <- statmod_spec(y ~ gas(1, 1, time = t, n_start = 1), distributions7::gaussian1_distrib(), d)
  expect_length(structural_multistarts(sp1, statmod_design(sp1)), 0L)
})

test_that("the intercept-only fit is computed once per call and gives the same fit", {
  d <- gas_series(3, n = 300)
  fam <- distributions7::gaussian1_distrib()
  fit_fun <- statmod_intercepts_fit
  calls <- 0L
  local_mocked_bindings(statmod_intercepts_fit = function(spec) {
    calls <<- calls + 1L
    fit_fun(spec)
  })
  f_memo <- statmod(y ~ 1 | sigma ~ gas(1, 1, time = t), fam, d,
                    outer_criterion = ml())
  expect_equal(calls, 1L)
  expect_null(statmod_memo$store)
  # without the memory: every reading fits again
  local_mocked_bindings(statmod_intercepts = function(spec) {
    calls <<- calls + 1L
    fit_fun(spec)
  })
  calls <- 0L
  f_plain <- statmod(y ~ 1 | sigma ~ gas(1, 1, time = t), fam, d,
                     outer_criterion = ml())
  expect_gt(calls, 1L)
  expect_identical(coef(f_memo), coef(f_plain))
  expect_identical(as.numeric(logLik(f_memo)), as.numeric(logLik(f_plain)))
})

test_that("the gradient check runs once per site within a fit, and the fit is unchanged", {
  d <- gas_series(3, n = 300)
  fam <- distributions7::gaussian1_distrib()
  sites <- character(0)
  wgc <- with_gradient_check
  local_mocked_bindings(with_gradient_check = function(site, expr) {
    sites <<- c(sites, site)
    wgc(site, expr)
  })
  checks <- logical(0)
  local_mocked_bindings(check_gradient_consistency = function(fn, gr, par) {
    checks <<- c(checks, isTRUE(getOption("optimizers7.check_gradient", TRUE)))
    invisible(TRUE)
  }, .package = "optimizers7")
  f_once <- statmod(y ~ 1 | sigma ~ gas(1, 1, time = t), fam, d, outer_criterion = ml())
  expect_gt(length(sites), length(unique(sites)))
  expect_true(any(checks) && any(!checks))
  expect_null(statmod_memo$store)
  expect_true(isTRUE(getOption("optimizers7.check_gradient", TRUE)))
  # every check, as before
  local_mocked_bindings(with_gradient_check = function(site, expr) expr)
  f_all <- statmod(y ~ 1 | sigma ~ gas(1, 1, time = t), fam, d, outer_criterion = ml())
  expect_identical(coef(f_once), coef(f_all))
  expect_identical(as.numeric(logLik(f_once)), as.numeric(logLik(f_all)))
})

test_that("outside a call of statmod() every gradient check runs", {
  expect_null(statmod_memo$store)
  seen <- logical(0)
  for (i in 1:2) with_gradient_check("x", seen <- c(seen, isTRUE(getOption("optimizers7.check_gradient", TRUE))))
  expect_identical(seen, c(TRUE, TRUE))
})
