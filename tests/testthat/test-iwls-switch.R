# An auto_info() settled on the expected information continues on the observed
# one after the `switch_after` scoring steps of iwls().

test_that("an auto iwls() on an exact family switches, a named one does not", {
  g <- iwls_resolve(iwls(), distributions7::gaussian1_distrib())
  expect_identical(info_label(g@hessian), "expected")
  expect_true(g@switching)
  expect_identical(g@switch_after, 10)
  p <- iwls_resolve(iwls(), distributions7::pig1_distrib())
  expect_false(p@switching)
  e <- iwls_resolve(iwls(), distributions7::gaussian1_distrib(),
                    expected_info())
  expect_false(e@switching)
  expect_false(iwls()@switching)
  # the constants are arguments
  expect_identical(iwls(switch_after = 4)@switch_after, 4)
  expect_identical(iwls(switch_rate = 0.7)@switch_rate, 0.7)
  expect_true(is.na(iwls(switch_rate = NA)@switch_rate))
  expect_error(iwls(switch_after = -1), "switch_after")
  expect_error(iwls(switch_rate = "a"), "switch_rate")
})

test_that("a scoring run that converges slowly moves to the observed pieces", {
  # the objective's curvature is Q, the 'expected' pieces a diagonal matrix
  # that differs from it, so scoring on them converges only linearly
  Q <- matrix(c(4, 3.5, 3.5, 4), 2)
  m <- c(1, -2)
  obj <- list(fn = function(b) 0.5 * sum((b - m) * as.numeric(Q %*% (b - m))),
              gr = function(b) as.numeric(Q %*% (b - m)))
  expd <- function(b) list(R = NULL, C = NULL, A = diag(diag(Q)))
  obsd <- function(b) list(R = NULL, C = NULL, A = Q)
  meth <- iwls_resolve(iwls(tol = 1e-10), NULL, expected_info())
  slow <- iwls_fit(obj, c(0, 0), meth, 1, expd)
  # the count alone: the rate rule off, the switch comes after three steps
  local({
    meth_na <- iwls_resolve(iwls(tol = 1e-10, switch_rate = NA), NULL,
                            expected_info())
    fast <- iwls_fit(obj, c(0, 0), meth_na, 1, expd, switch_at = obsd,
                     switch_after = 3)
    expect_true(fast$converged)
    expect_equal(fast$par, m, tolerance = 1e-8)
    expect_identical(fast$switched, 4L)
    expect_lt(fast$iterations, slow$iterations)
  })
  expect_identical(slow$switched, 0L)
  # the rate alone: scoring contracts by 3.5/4 a step here, above the rate,
  # so the run switches at the third step with no count to reach
  rate <- iwls_fit(obj, c(0, 0), meth, 1, expd, switch_at = obsd,
                   switch_after = Inf)
  expect_true(rate$converged)
  expect_equal(rate$par, m, tolerance = 1e-8)
  expect_identical(rate$switched, 3L)
  # and where scoring contracts fast, the rate does not fire
  Q2 <- matrix(c(4, 0.5, 0.5, 4), 2)
  obj2 <- list(fn = function(b) 0.5 * sum((b - m) * as.numeric(Q2 %*% (b - m))),
               gr = function(b) as.numeric(Q2 %*% (b - m)))
  quick <- iwls_fit(obj2, c(0, 0), meth, 1,
                    function(b) list(R = NULL, C = NULL, A = diag(diag(Q2))),
                    switch_at = function(b) list(R = NULL, C = NULL, A = Q2),
                    switch_after = Inf)
  # (it lands on the minimum exactly, where the stall guard rather than the
  # rule ends the run, so the flag is not what is read here)
  expect_equal(quick$par, m, tolerance = 1e-8)
  expect_identical(quick$switched, 0L)
})

test_that("the rate rule reads the last two contractions", {
  expect_false(iwls_slow(c(1, 0.9), 0.5))
  expect_true(iwls_slow(c(1, 0.9, 0.8), 0.5))
  expect_false(iwls_slow(c(1, 0.9, 0.1), 0.5))
  expect_false(iwls_slow(c(1, 0.9, 0.8), NA_real_))
})

test_that("a gaussian with smooths in the mean and in sigma is fitted", {
  cw <- as.data.frame(ChickWeight)
  cw$Diet <- factor(cw$Diet, ordered = FALSE)
  f <- weight ~ Diet + s(Time, bspline_smooth(k = 8)) |
    sigma ~ s(Time, bspline_smooth(k = 8))
  fit <- statmod(f, distrib = gaussian1_distrib(), data = cw)
  ref <- statmod(f, distrib = gaussian1_distrib(), data = cw,
                 inner_information = observed_info())
  expect_identical(statmod_certificate(fit)$state, "converged")
  expect_equal(as.numeric(logLik(fit, loglik_type = "marginal")),
               as.numeric(logLik(ref, loglik_type = "marginal")), tolerance = 1e-7)
  expect_equal(hyper(fit)$estimate, hyper(ref)$estimate, tolerance = 1e-4)
  # without either switch the first inner fit exhausts a budget of 100
  # scoring steps and the criterion is unavailable at its start
  expect_error(statmod(f, distrib = gaussian1_distrib(), data = cw,
                       inner_optimizer = iwls(maxit = 100, switch_after = Inf,
                                              switch_rate = NA)),
               "unavailable at the starting hyperparameters")
  # the rate rule alone is enough at that budget (an infinite count would
  # turn the switch off altogether, as a named expected_info() does)
  byrate <- statmod(f, distrib = gaussian1_distrib(), data = cw,
                    inner_optimizer = iwls(maxit = 100, switch_after = 1e6))
  expect_equal(as.numeric(logLik(byrate, loglik_type = "marginal")),
               as.numeric(logLik(ref, loglik_type = "marginal")), tolerance = 1e-7)
})

test_that("a run that converges within the switch is unchanged", {
  cw <- as.data.frame(ChickWeight)
  f <- weight ~ s(Time, bspline_smooth(k = 8)) + random(~ 1 | Chick)
  a <- statmod(f, distrib = gaussian1_distrib(), data = cw)
  b <- statmod(f, distrib = gaussian1_distrib(), data = cw,
               inner_optimizer = iwls(switch_after = Inf, switch_rate = NA))
  expect_identical(a@coefficients, b@coefficients)
  expect_identical(hyper(a)$estimate, hyper(b)$estimate)
})
