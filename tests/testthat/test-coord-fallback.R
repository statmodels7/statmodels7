# A kinked block fitted on the OBSERVED information falls back to the expected
# one where the observed curvature is not positive, rather than leaving the
# compiled coordinate descent for the proximal route.

set.seed(3)
nf <- 400L
xf <- matrix(stats::rnorm(nf * 8L), nf, 8L)
zf <- matrix(stats::rnorm(nf * 4L), nf, 4L)
df <- data.frame(y = stats::rnbinom(nf, mu = exp(1 + xf[, 1] / 2 - xf[, 2] / 3),
                                     size = exp(1 + 0.5 * zf[, 1])))
df$X <- xf
df$Z <- zf
ff <- y ~ lasso(X, lambda = 5) | theta ~ lasso(Z, lambda = 5)

test_that("the theta block of a negbin2 stays on the coordinate descent", {
  skip_on_cran()
  ns <- asNamespace("statmodels7")
  cw <- get("coord_working", ns)
  cf <- get("coord_fit", ns)
  seen <- new.env()
  seen$observed_null <- 0L
  seen$cd_null <- 0L
  seen$cd_theta <- 0L
  local_mocked_bindings(
    coord_working = function(spec, ep, coef, design, p, expected, approx) {
      r <- cw(spec, ep, coef, design, p, expected, approx)
      if (is.null(r) && !expected) seen$observed_null <- seen$observed_null + 1L
      r
    },
    coord_fit = function(obj, beta, block, ...) {
      r <- cf(obj, beta, block, ...)
      if (identical(block$param, "theta")) {
        seen$cd_theta <- seen$cd_theta + 1L
        if (is.null(r)) seen$cd_null <- seen$cd_null + 1L
      }
      r
    }
  )
  fo <- statmod(ff, distributions7::negbin2_distrib(), df,
                inner_optimizer = iwls(hessian = "observed"))
  # the premise: the observed curvature really is not positive somewhere,
  # so the fallback is exercised and the test cannot pass vacuously
  expect_gt(seen$observed_null, 0L)
  expect_gt(seen$cd_theta, 0L)
  expect_identical(seen$cd_null, 0L)

  # the mode does not depend on the curvature the steps are weighted by; the
  # two runs stop at the inner rule's own resolution, about 2e-6 relative
  fe <- statmod(ff, distributions7::negbin2_distrib(), df,
                inner_optimizer = iwls(hessian = "expected"))
  expect_equal(fo@coefficients, fe@coefficients, tolerance = 1e-5)
  expect_identical(unlist(fo@coefficients) != 0, unlist(fe@coefficients) != 0)
})
