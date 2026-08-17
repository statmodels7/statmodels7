# The exact gradient with the EXPECTED information.
#
# reml("expected") had no gradient at all, so its search was derivative-free
# and its evaluation count grew with the number of hyperparameters: measured on
# three smooths and a random effect, 965 evaluations and 481 s against the
# observed route's 13 and 5.9, for the same hyperparameters. What was missing
# is dK/dbeta, which is -l''' with the observed information and
# -dE[l'']/deta with the expected one -- and those are different objects,
# because differentiating an expectation moves the measure as well as the
# integrand.

set.seed(41)
ng <- 300
dge <- data.frame(x = runif(ng, -2, 2), z = runif(ng))
dge$y <- rgamma(ng, shape = 4,
                rate = 4 / exp(0.4 * sin(1.4 * dge$x) + 1))

# the criterion as a function of the free hyperparameters, the mode refitted
# from a FIXED start at each -- the same shape as test-outer-gradient.R's
# harness, with the family and the information as parameters
crit_of_eta_e <- function(formula, data, distrib, method) {
  fit0 <- statmod(formula, distrib, data, outer_criterion = NULL)
  spec <- fit0@spec
  design <- statmod_design(spec)
  blocks <- statmod_blocks(spec, design)
  idx <- outer_hyper_index(spec, blocks)
  basis <- integrated_basis(spec, design, method@kind)
  at <- function(eta) {
    hy <- eta_to_hyper(eta, idx, fit0@hyper)
    list(hy = hy, cf = fit_at_hyper(formula, distrib, data, hy,
                                    iwls(tol = 1e-9))$coefficients)
  }
  list(idx = idx,
       fn = function(eta) {
         a <- at(eta)
         statmod_marginal(spec, design, a$cf, a$hy, method, basis = basis)$value
       },
       gr = function(eta) {
         a <- at(eta)
         statmod_marginal_grad(spec, design, a$cf, a$hy, method, idx, basis)
       },
       eta0 = hyper_to_eta(fit0@hyper, idx),
       spec = spec, design = design, idx = idx)
}

test_that("the expected route's gradient matches numDeriv on a gamma", {
  skip_if_not_installed("numDeriv")
  # a gamma, not a poisson: at a CANONICAL link the expected and observed
  # informations are the same matrix, so a poisson would confirm a tautology.
  #
  # The tolerance is the observed route's own, deliberately. A first version of
  # this measured 1.9e-03 and the honest reading was not a looser tolerance:
  # the gap was FLAT in the inner tolerance while the mode's score fell by two
  # decades, which is what says a reference is not the weak side. It was the
  # mode's movement being read off the criterion's matrix instead of the
  # penalized likelihood's; corrected, the same comparison is 1.8e-08.
  h <- crit_of_eta_e(y ~ s(x, k = 10), dge, distributions7::gamma1_distrib(),
                     reml(hessian = "expected"))
  eta <- h$eta0 + 0.4
  expect_equal(h$gr(eta), numDeriv::grad(h$fn, eta), tolerance = 1e-6)
})

test_that("it matches with a second penalized equation", {
  skip_if_not_installed("numDeriv")
  # two equations make dE[l_ab]/deta genuinely three-index, and the component
  # is symmetric in (a, b) and NOT in the parameter differentiated in, so a
  # version keyed like the third derivative would pass the previous test and
  # fail this one
  h <- crit_of_eta_e(y ~ s(x, k = 8) | phi ~ s(z, k = 6), dge,
                     distributions7::gamma1_distrib(),
                     reml(hessian = "expected"))
  eta <- h$eta0 + c(0.3, -0.4)
  expect_equal(h$gr(eta), numDeriv::grad(h$fn, eta), tolerance = 1e-5)
})

test_that("it holds under ml, on the range space", {
  skip_if_not_installed("numDeriv")
  # ml() projects the curvature onto the coefficients' range basis while the
  # MODE goes on moving in the full space, so after the correction two
  # different matrices are read for two different reasons -- the determinant's
  # is the criterion's and projected, the mode's is the penalized likelihood's
  # and full. This is where a version that confused them shows it worst, and
  # every other test in this file is reml.
  h <- crit_of_eta_e(y ~ s(x, k = 10), dge, distributions7::gamma1_distrib(),
                     ml(hessian = "expected"))
  eta <- h$eta0 + 0.4
  expect_equal(h$gr(eta), numDeriv::grad(h$fn, eta), tolerance = 1e-6)

  h2 <- crit_of_eta_e(y ~ s(x, k = 8) | phi ~ s(z, k = 6), dge,
                      distributions7::gamma1_distrib(),
                      ml(hessian = "expected"))
  eta2 <- h2$eta0 + c(0.3, -0.4)
  expect_equal(h2$gr(eta2), numDeriv::grad(h2$fn, eta2), tolerance = 1e-5)
})

test_that("the two routes are one route where the link is canonical", {
  # The control, and the one place where this phase is a tautology: with a log
  # link a poisson's expected information IS its observed one. Stated at a
  # FIXED point it is an exact identity rather than a comparison of two
  # searches, so it is asked at machine precision -- and it is the check that
  # would catch the expected array being read with the wrong key, which no
  # tolerance on a fitted hyperparameter is sharp enough to see.
  set.seed(42)
  dp <- data.frame(x = runif(400, -2, 2))
  dp$y <- rpois(400, exp(0.5 * sin(1.4 * dp$x) + 1))
  f <- y ~ s(x, k = 8)
  d <- distributions7::poisson_distrib()
  h <- crit_of_eta_e(f, dp, d, reml(hessian = "expected"))
  eta <- h$eta0 + 0.35
  hy <- eta_to_hyper(eta, h$idx, statmod_hyper_start(h$spec, h$design))
  cf <- fit_at_hyper(f, d, dp, hy, iwls(tol = 1e-9))$coefficients
  basis <- integrated_basis(h$spec, h$design, "reml")
  go <- statmod_marginal_grad(h$spec, h$design, cf, hy,
                              reml(hessian = "observed"), h$idx, basis)
  ge <- statmod_marginal_grad(h$spec, h$design, cf, hy,
                              reml(hessian = "expected"), h$idx, basis)
  expect_equal(go, ge, tolerance = 1e-10)
  # and the two whole fits therefore land on the same optimum, to whatever
  # their own searches resolve
  fo <- statmod(f, d, dp, outer_criterion = reml(hessian = "observed"))
  fe <- statmod(f, d, dp, outer_criterion = reml(hessian = "expected"))
  expect_equal(unlist(lapply(fo@hyper, unlist)),
               unlist(lapply(fe@hyper, unlist)), tolerance = 1e-4)
  expect_equal(as.numeric(logLik(fo)), as.numeric(logLik(fe)),
               tolerance = 1e-6)
})

test_that("the mode's movement is read off the penalized likelihood", {
  # The defect the first version of this file measured: v = db/dt solves
  # (H_obs + S) v = -d2rho/dbeta dt whatever matrix the determinant is of, and
  # reading it off the criterion's K is wrong by the gap between the two
  # informations -- systematic, and shrinking with n as they converge.
  skip_if_not_installed("numDeriv")
  h300 <- crit_of_eta_e(y ~ s(x, k = 10), dge,
                        distributions7::gamma1_distrib(),
                        reml(hessian = "expected"))
  eta <- h300$eta0 + 0.4
  g <- h300$gr(eta)
  fd <- numDeriv::grad(h300$fn, eta)
  # the wrong reading gave 1.9e-03 here; anything of that size is the
  # conflation coming back
  expect_lt(max(abs(g - fd) / abs(fd)), 1e-5)
})

test_that("the expected route is admitted only where the family answers", {
  h <- crit_of_eta_e(y ~ s(x, k = 8), dge, distributions7::gamma1_distrib(),
                     reml(hessian = "expected"))
  expect_true(outer_gradient_ok(h$spec, h$design, h$idx,
                                reml(hessian = "expected"), 1L))
  # order 2 is not extended: the criterion's own second derivative would want
  # the next order of the same object, and lbfgs on the exact gradient buys
  # most of what newton would
  expect_false(outer_gradient_ok(h$spec, h$design, h$idx,
                                 reml(hessian = "expected"), 2L))
  # and a family that approximates its expected information leaves the search
  # derivative-free rather than reporting a gradient it cannot compute
  expect_false(expected_deriv_ok(distributions7::skewt_distrib()))
  expect_true(expected_deriv_ok(distributions7::gamma1_distrib()))
})

test_that("the criterion's own information is the one that is factorized", {
  # ctx_penalized() hard-coded the OBSERVED information, which was right while
  # the exact gradient ran on no other route and became a gradient of the
  # wrong function the moment the expected route was admitted.
  h <- crit_of_eta_e(y ~ s(x, k = 8), dge, distributions7::gamma1_distrib(),
                     reml(hessian = "expected"))
  eta <- h$eta0 + 0.2
  hy <- eta_to_hyper(eta, h$idx, statmod_hyper_start(h$spec, h$design))
  cf <- fit_at_hyper(y ~ s(x, k = 8), distributions7::gamma1_distrib(), dge,
                     hy)$coefficients
  ctx <- outer_context(h$spec, h$design, cf, hy)
  ko <- ctx_penalized(ctx, h$spec, h$design, cf, hy, FALSE)
  ke <- ctx_penalized(ctx, h$spec, h$design, cf, hy, TRUE)
  expect_false(isTRUE(all.equal(as.matrix(ko$K), as.matrix(ke$K))))
  expect_false(isTRUE(all.equal(ko$logdet, ke$logdet)))
})
