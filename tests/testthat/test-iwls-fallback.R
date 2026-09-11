# iwls(hessian = "auto"), settled against the family, and the expected
# information taking the step where the observed one cannot.

test_that("an auto iwls() is settled by the family and a named one is not", {
  a <- iwls()
  expect_identical(a@hessian, "auto")
  g <- iwls_resolve(a, distributions7::gaussian1_distrib())
  expect_identical(g@hessian, "expected")
  expect_false(g@fallback)
  for (d in list(distributions7::pig1_distrib(), distributions7::skewnormal1_distrib(),
                 distributions7::truncated(distributions7::gaussian1_distrib(), lower = 0))) {
    r <- iwls_resolve(a, d)
    expect_identical(r@hessian, "observed")
    expect_true(r@fallback)
  }
  e <- iwls_resolve(iwls(hessian = "expected"), distributions7::pig1_distrib())
  expect_identical(e@hessian, "expected")
  expect_false(e@fallback)
  o <- iwls_resolve(iwls(hessian = "observed"), distributions7::pig1_distrib())
  expect_identical(o@hessian, "observed")
  expect_false(o@fallback)
  # a method nobody settled is refused where the curvature is read
  expect_error(inner_settings(iwls()), "settled by the family")
  expect_true(inner_settings(iwls(), distributions7::gaussian1_distrib())$expected)
  expect_false(inner_settings(iwls(), distributions7::pig1_distrib())$expected)
})

test_that("a fit records the curvature its family settled", {
  set.seed(11)
  n <- 150
  d <- data.frame(x = stats::runif(n))
  mu <- exp(1 + 0.5 * d$x)
  lam <- distributions7::distrib_rng(distributions7::invgauss1_distrib(), n,
                                     list(mu = mu, phi = 0.8 / mu))
  d$y <- stats::rpois(n, lam)
  fp <- statmod(y ~ x, distributions7::pig1_distrib(), d)
  expect_identical(fp@methods$smooth@hessian, "observed")
  expect_true(fp@methods$smooth@fallback)
  fe <- statmod(y ~ x, distributions7::pig1_distrib(), d,
                inner = iwls(hessian = "expected"))
  expect_identical(fe@methods$smooth@hessian, "expected")
  expect_false(fe@methods$smooth@fallback)
  # the same maximum, reached through either curvature
  expect_equal(as.numeric(logLik(fp)), as.numeric(logLik(fe)), tolerance = 1e-6)

  dg <- data.frame(x = d$x, y = 1 + 2 * d$x + stats::rnorm(n))
  fg <- statmod(y ~ x, distributions7::gaussian1_distrib(), dg)
  expect_identical(fg@methods$smooth@hessian, "expected")
  expect_false(fg@methods$smooth@fallback)
})

test_that("the expected pieces take the step where the observed ones cannot", {
  # A quadratic, so a step on the right curvature lands on the minimum and
  # every other route has to be taken on purpose.
  Q <- matrix(c(4, 1, 1, 3), 2)
  m <- c(1, -2)
  toy <- function(bad_outside = Inf) {
    list(fn = function(b) {
           if (max(abs(b)) > bad_outside) stop("outside the region the objective is defined on")
           0.5 * sum((b - m) * as.numeric(Q %*% (b - m)))
         },
         gr = function(b) as.numeric(Q %*% (b - m)))
  }
  pd <- function(b) list(R = NULL, C = NULL, A = Q)

  # an observed curvature that is not positive definite is not used
  neg <- function(b) list(R = NULL, C = NULL, A = -Q)
  r1 <- iwls_fit(toy(), c(0, 0), iwls(hessian = "observed"), 1, neg, backup_at = pd)
  expect_true(r1$converged)
  expect_equal(r1$par, m, tolerance = 1e-6)
  expect_gt(r1$fallback[["indefinite"]], 0L)
  expect_identical(r1$fallback[["search"]], 0L)

  # a positive definite one whose step is 1e12 times too long, beyond what
  # five halvings can shorten: the search fails and the expected pieces step
  far <- function(b) list(R = NULL, C = NULL, A = Q * 1e-12)
  short <- iwls(hessian = "observed", step_halving = 5)
  r2 <- iwls_fit(toy(), c(0, 0), short, 1, far, backup_at = pd)
  expect_true(r2$converged)
  expect_equal(r2$par, m, tolerance = 1e-6)
  expect_gt(r2$fallback[["search"]], 0L)
  expect_identical(r2$fallback[["error"]], 0L)

  # a trial point whose objective raises is read as rejected; with no
  # fallback the error ends the fit, as it always has
  r3 <- iwls_fit(toy(bad_outside = 10), c(0, 0), short, 1, far, backup_at = pd)
  expect_true(r3$converged)
  expect_equal(r3$par, m, tolerance = 1e-6)
  expect_gt(r3$fallback[["error"]], 0L)
  expect_error(iwls_fit(toy(bad_outside = 10), c(0, 0), short, 1, far),
               "outside the region")

  # and a run that never needed it reports none
  r4 <- iwls_fit(toy(), c(0, 0), iwls(hessian = "observed"), 1, pd, backup_at = pd)
  expect_identical(r4$fallback, c(indefinite = 0L, search = 0L, error = 0L))
})
