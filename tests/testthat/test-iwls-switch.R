# iwls(hessian = "auto") settled on the expected information continues on the
# observed one after iwls_switch_after() scoring steps.

test_that("an auto iwls() on an exact family switches, a named one does not", {
  g <- iwls_resolve(iwls(), distributions7::gaussian1_distrib())
  expect_identical(g@hessian, "expected")
  expect_identical(g@switch_after, iwls_switch_after())
  p <- iwls_resolve(iwls(), distributions7::pig1_distrib())
  expect_identical(p@switch_after, Inf)
  e <- iwls_resolve(iwls(hessian = "expected"), distributions7::gaussian1_distrib())
  expect_identical(e@switch_after, Inf)
  expect_identical(iwls()@switch_after, Inf)
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
  meth <- iwls(hessian = "expected", tol = 1e-10)
  slow <- iwls_fit(obj, c(0, 0), meth, 1, expd)
  fast <- iwls_fit(obj, c(0, 0), meth, 1, expd, switch_at = obsd,
                   switch_after = 3)
  expect_true(fast$converged)
  expect_equal(fast$par, m, tolerance = 1e-8)
  expect_identical(fast$switched, 4L)
  expect_identical(slow$switched, 0L)
  expect_lt(fast$iterations, slow$iterations)
})

test_that("a gaussian with smooths in the mean and in sigma is fitted", {
  cw <- as.data.frame(ChickWeight)
  cw$Diet <- factor(cw$Diet, ordered = FALSE)
  f <- weight ~ Diet + s(Time, bspline_smooth(k = 8)) |
    sigma ~ s(Time, bspline_smooth(k = 8))
  fit <- statmod(f, distrib = gaussian1_distrib(), data = cw)
  ref <- statmod(f, distrib = gaussian1_distrib(), data = cw,
                 inner_optimizer = iwls(hessian = "observed"))
  expect_identical(statmod_certificate(fit)$state, "converged")
  expect_equal(as.numeric(logLik(fit, type = "marginal")),
               as.numeric(logLik(ref, type = "marginal")), tolerance = 1e-7)
  expect_equal(hyper(fit)$estimate, hyper(ref)$estimate, tolerance = 1e-4)
  # without the switch the first inner fit exhausts its budget of 100
  # scoring steps and the criterion is unavailable at its start
  local_mocked_bindings(iwls_switch_after = function() Inf)
  expect_error(statmod(f, distrib = gaussian1_distrib(), data = cw),
               "unavailable at the starting hyperparameters")
})

test_that("a run that converges within the switch is unchanged", {
  cw <- as.data.frame(ChickWeight)
  f <- weight ~ s(Time, bspline_smooth(k = 8)) + random(~ 1 | Chick)
  a <- statmod(f, distrib = gaussian1_distrib(), data = cw)
  local_mocked_bindings(iwls_switch_after = function() Inf)
  b <- statmod(f, distrib = gaussian1_distrib(), data = cw)
  expect_identical(a@coefficients, b@coefficients)
  expect_identical(hyper(a)$estimate, hyper(b)$estimate)
})
