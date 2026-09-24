# The exact outer Hessian on the EXPECTED information.
#
# Order 2 reads the second derivative of the expected information in the
# predictors, which distributions7 writes out family by family. The assembly
# reads two matrices: the criterion's K = H_exp + S for the determinant, and
# J = H_obs + S for how the mode moves. Reading one for both is the defect this
# file exists to catch.

set.seed(47)
n_eh <- 400
d_eh <- data.frame(x = runif(n_eh, -2, 2), z = runif(n_eh))
f_eh <- exp(0.4 * sin(1.4 * d_eh$x) + 1)
d_eh$yg <- rgamma(n_eh, shape = 4, rate = 4 / f_eh)
d_eh$yn <- rnbinom(n_eh, size = 3, mu = f_eh)
d_eh$yr <- sin(1.4 * d_eh$x) + rnorm(n_eh, 0, 0.5)

# the criterion's exact gradient and Hessian as functions of the free
# hyperparameters, the mode refitted from a FIXED start and then POLISHED by
# Newton steps: the inner fit stops on its stall guard with a residual score
# that is a smooth function of the hyperparameter, and differenced that is a
# bias flat in the step (see helper-hyper.R)
eh_harness <- function(formula, data, distrib, method) {
  fit0 <- statmod(formula, distrib, data, outer_criterion = NULL)
  spec <- fit0@spec
  design <- statmod_design(spec)
  idx <- outer_hyper_index(spec, statmod_blocks(spec, design))
  basis <- integrated_basis(spec, design, method@kind)
  at <- function(eta) {
    hy <- eta_to_hyper(eta, idx, fit0@hyper)
    list(hy = hy, cf = fit_at_hyper(formula, distrib, data, hy,
                                    iwls(tol = 1e-10, maxit = 500),
                                    polish = TRUE)$coefficients)
  }
  list(spec = spec, design = design, idx = idx, basis = basis, at = at,
       gr = function(eta) {
         a <- at(eta)
         statmod_marginal_grad(spec, design, a$cf, a$hy, method, idx, basis)
       },
       he = function(eta, m = method) {
         a <- at(eta)
         statmod_marginal_hess(spec, design, a$cf, a$hy, m, idx, basis)
       },
       eta0 = hyper_to_eta(fit0@hyper, idx))
}

eh_fd <- function(h, eta, step) {
  nh <- length(eta)
  H <- matrix(0, nh, nh)
  for (m in seq_len(nh)) {
    ep <- eta; ep[m] <- ep[m] + step
    em <- eta; em[m] <- em[m] - step
    H[, m] <- (h$gr(ep) - h$gr(em)) / (2 * step)
  }
  (H + t(H)) / 2
}

test_that("the expected route's Hessian is the derivative of its gradient", {
  skip_on_cran()
  # Two penalized equations, so the off-diagonal entry and the keying of
  # d2E[l_ab]/deta_c deta_d by two separate pairs are both exercised. With the
  # mode polished the gap to a difference of the exact gradient is 8.2e-08 at
  # h = 1e-3 and to numDeriv 2.0e-10, where the assembly that traced the
  # OBSERVED third and fourth derivatives against the expected matrix was out
  # by 1e-4 to 3e-4. Unpolished the reference sat at a flat 4e-06, which is
  # what the earlier bound of 2e-5 was set against.
  h <- eh_harness(yg ~ s(x, bspline_smooth(k = 8)) |
                    phi ~ s(z, bspline_smooth(k = 6)),
                  d_eh, distributions7::gamma1_distrib(), reml("expected"))
  expect_true(outer_gradient_ok(h$spec, h$design, h$idx, reml("expected"), 2L))
  eta <- h$eta0 + c(0.3, -0.4)
  Ha <- h$he(eta)
  Hf <- eh_fd(h, eta, 1e-3)
  expect_lt(max(abs(Ha - Hf)) / max(abs(Hf)), 1e-6)
  # and it is not the observed route's Hessian under another name: the two
  # criteria differ at this sample size by far more than the gap above
  Ho <- h$he(eta, reml("observed"))
  expect_gt(max(abs(Ha - Ho)) / max(abs(Hf)), 1e-3)
})

test_that("it holds under ml and on a sum over the support", {
  skip_on_cran()
  # negbin2's d2E[l_theta_theta] is a sum over the support whose mass moves
  # with the parameters; ml() projects the determinant onto the range space
  # while the mode moves in the full one
  h <- eh_harness(yn ~ s(x, bspline_smooth(k = 8)), d_eh,
                  distributions7::negbin2_distrib(), ml("expected"))
  expect_true(outer_gradient_ok(h$spec, h$design, h$idx, ml("expected"), 2L))
  eta <- h$eta0 + 0.4
  Ha <- h$he(eta)
  Hf <- eh_fd(h, eta, 1e-3)
  expect_lt(max(abs(Ha - Hf)) / max(abs(Hf)), 1e-6)
})

test_that("with the dispersion modelled it converges O(h^2) onto the assembly", {
  skip_on_cran()
  # A negative binomial with a smooth in each equation is where the
  # UNPOLISHED reference stopped hardest: the inner fit leaves a residual score
  # of 1.6e-04 there, and differenced it read a gap of 9.5e-06 FLAT in the
  # step, the same on the observed route, which looks like a missing term and
  # is not. Polished, the gap falls by the square of the step (5.4e-07 and
  # 6.0e-08 at 3e-3 and 1e-3 on 600 observations). A floor could never pass
  # the ratio asserted here; a missing term would not either.
  h <- eh_harness(yn ~ s(x, bspline_smooth(k = 8)) |
                    theta ~ s(z, bspline_smooth(k = 6)),
                  d_eh, distributions7::negbin2_distrib(), reml("expected"))
  expect_true(outer_gradient_ok(h$spec, h$design, h$idx, reml("expected"), 2L))
  eta <- h$eta0 + c(0.3, -0.4)
  Ha <- h$he(eta)
  gap <- vapply(c(3e-3, 1e-3), function(s) {
    Hf <- eh_fd(h, eta, s)
    max(abs(Ha - Hf)) / max(abs(Hf))
  }, 0)
  expect_lt(gap[[2L]], 1e-6)
  expect_gt(gap[[1L]] / gap[[2L]], 5)
  expect_lt(gap[[1L]] / gap[[2L]], 15)
})

test_that("a family that does not write the second derivative is differenced", {
  skip_on_cran()
  # The three callers that reach the assembly without outer_gradient_ok()'s
  # gate -- statmod_hyper_vcov(), hyper_correction(), statmod_edf_correction()
  # -- used to get the mixed assembly on the expected route; they get the
  # stencil of the exact gradient now, by identity. Every family whose
  # expected information is written out answers the second derivative, so
  # the refusal is injected rather than taken from a family.
  local_mocked_bindings(
    distrib_d2expected_hessian = function(...) stop("no analytic second derivative"),
    .package = "distributions7")
  h <- eh_harness(yr ~ s(x, bspline_smooth(k = 8)), d_eh,
                  distributions7::gaussian2_distrib(), reml("expected"))
  expect_false(outer_gradient_ok(h$spec, h$design, h$idx, reml("expected"), 2L))
  a <- h$at(h$eta0 + 0.3)
  got <- statmod_marginal_hess(h$spec, h$design, a$cf, a$hy, reml("expected"),
                               h$idx, h$basis)
  want <- statmod_hess_stencil(h$spec, h$design, a$cf, a$hy, reml("expected"),
                               h$idx, h$basis)
  expect_identical(got, want)
})

test_that("the search on the expected route keeps lbfgs", {
  skip_on_cran()
  # D5: newton() on the exact Hessian reaches the criterion lbfgs() reaches and
  # is slower or level in wall time on the six models measured
  fit <- statmod(yg ~ s(x, bspline_smooth(k = 8)),
                 distributions7::gamma1_distrib(), d_eh,
                 outer_criterion = reml("expected"))
  expect_match(class(fit@methods$search)[1L], "Lbfgs", fixed = TRUE)
  spec <- fit@spec
  design <- statmod_design(spec)
  expect_false(outer_newton_ok(spec, design, TRUE, reml("expected")))
  expect_true(outer_newton_ok(spec, design, TRUE, reml("observed")))
  cert <- statmod_certificate(fit)
  expect_identical(cert$curvature, "analytic")
})

test_that("the families of the elementary kernels converge O(h^2) onto the assembly", {
  skip_on_cran()
  # One family per shape of kernel: one parameter (a rational term), a gamma
  # by its variance (the polygamma remainders), and a beta by its two shapes
  # (the polygamma differences), the second parameter modelled where there is
  # one. Measured, the gap falls by exactly 9 between the steps 3e-3 and 1e-3
  # on all sixteen families of distributions7 0.62.0, onto 5e-9 .. 2e-7.
  set.seed(11)
  n <- 400; x <- runif(n); z <- runif(n); f <- sin(2 * pi * x)
  cases <- list(
    list(y ~ s(x, bspline_smooth(k = 8)), distributions7::bernoulli_distrib(),
         data.frame(x = x, z = z, y = rbinom(n, 1, plogis(0.3 + 0.5 * f)))),
    list(y ~ s(x, bspline_smooth(k = 8)) | sigma2 ~ s(z, bspline_smooth(k = 6)),
         distributions7::gamma2_distrib(),
         data.frame(x = x, z = z, y = rgamma(n, 4, 4 / exp(1 + 0.5 * f)))),
    list(y ~ s(x, bspline_smooth(k = 8)) | beta ~ s(z, bspline_smooth(k = 6)),
         distributions7::beta2_distrib(),
         data.frame(x = x, z = z, y = rbeta(n, exp(1 + 0.5 * f), 4))))
  for (cs in cases) {
    h <- eh_harness(cs[[1]], cs[[3]], cs[[2]], reml("expected"))
    expect_true(outer_gradient_ok(h$spec, h$design, h$idx, reml("expected"), 2L))
    eta <- h$eta0 + seq(0.3, -0.4, length.out = length(h$eta0))
    Ha <- h$he(eta)
    gap <- vapply(c(3e-3, 1e-3), function(s) {
      Hf <- eh_fd(h, eta, s)
      max(abs(Ha - Hf)) / max(abs(Hf))
    }, 0)
    expect_lt(gap[[2L]], 1e-6)
    expect_gt(gap[[1L]] / gap[[2L]], 5)
    expect_lt(gap[[1L]] / gap[[2L]], 15)
  }
})
