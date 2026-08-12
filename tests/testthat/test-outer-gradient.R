# The exact gradient of the marginal criterion.

set.seed(31)
n <- 250
dg <- data.frame(x = runif(n, -2, 2), z = runif(n))
dg$y <- sin(1.4 * dg$x) + stats::rnorm(n, sd = 0.3)

# The criterion as a function of the free hyperparameters, refitting the
# coefficients at each: this is what the outer search minimizes, and
# differentiating it numerically shares no arithmetic with the assembly the
# gradient uses.
crit_of_eta <- function(formula, data, method, hyper0 = NULL) {
  # the base point is the PROBE value, not the estimated optimum: this
  # differentiates the criterion at a point where its gradient is not zero,
  # which is the only place a relative comparison against numDeriv means
  # anything
  fit0 <- statmod(formula, distributions7::gaussian1_distrib(), data,
                  hyper = hyper0, outer_criterion = NULL)
  spec <- fit0@spec
  design <- statmod_design(spec)
  blocks <- statmod_blocks(spec, design)
  idx <- outer_hyper_index(spec, blocks)
  basis <- integrated_basis(spec, design, method@kind)
  list(
    idx = idx,
    fn = function(eta) {
      hy <- eta_to_hyper(eta, idx, fit0@hyper)
      f <- statmod(formula, distributions7::gaussian1_distrib(), data,
                   hyper = as_hyper_list(hy))
      statmod_marginal(spec, design, f@coefficients, hy, method,
                       basis = basis)$value
    },
    gr = function(eta) {
      hy <- eta_to_hyper(eta, idx, fit0@hyper)
      f <- statmod(formula, distributions7::gaussian1_distrib(), data,
                   hyper = as_hyper_list(hy))
      statmod_marginal_grad(spec, design, f@coefficients, hy, method, idx,
                            basis)
    },
    eta0 = hyper_to_eta(fit0@hyper, idx)
  )
}

# statmod()'s hyper= takes a plain named list of named vectors, and a
# distribution parameter with no penalized term has no entry at all
as_hyper_list <- function(hy) {
  out <- lapply(hy, function(per_param)
    lapply(per_param, function(v) stats::setNames(as.numeric(v), names(v))))
  out[lengths(out) > 0L]
}

test_that("the gradient of a single smoothing parameter matches numDeriv", {
  skip_if_not_installed("numDeriv")
  h <- crit_of_eta(y ~ s(x, k = 10), dg, reml(hessian = "observed"))
  eta <- h$eta0 + 0.4
  expect_equal(h$gr(eta), numDeriv::grad(h$fn, eta), tolerance = 1e-6)
  # away from the start, where the mode has genuinely moved
  eta2 <- h$eta0 - 1.3
  expect_equal(h$gr(eta2), numDeriv::grad(h$fn, eta2), tolerance = 1e-6)
})

test_that("the gradient matches numDeriv with the scale modelled too", {
  # a second distribution parameter with its own design makes the third
  # derivative genuinely three-index: the mu_mu_sigma and mu_sigma_sigma
  # components enter u, and a version that summed only the diagonal would
  # pass the previous test and fail this one
  skip_if_not_installed("numDeriv")
  h <- crit_of_eta(y ~ s(x, k = 8) | sigma ~ z, dg,
                   reml(hessian = "observed"))
  eta <- h$eta0 + 0.2
  expect_equal(h$gr(eta), numDeriv::grad(h$fn, eta), tolerance = 1e-6)
})

test_that("the gradient matches numDeriv with several parameters at once", {
  skip_if_not_installed("numDeriv")
  set.seed(32)
  dt <- data.frame(x1 = runif(200, -1, 1), x2 = runif(200, -1, 1))
  dt$y <- dt$x1^2 + sin(3 * dt$x2) + stats::rnorm(200, sd = 0.3)
  h <- crit_of_eta(y ~ s(x1, k = 8) + s(x2, k = 8), dt,
                   reml(hessian = "observed"))
  eta <- h$eta0 + c(0.5, -0.6)
  g <- h$gr(eta)
  expect_length(g, 2L)
  # the reference refits the mode at every perturbation, so it carries the
  # inner tolerance; measured here, 2.2021349 against 2.2021371 and 2.1755400
  # against 2.1755434, which is 1.5e-6 relative. It passed at 1e-6 while the
  # fits started somewhere else, which was luck rather than accuracy.
  expect_equal(g, numDeriv::grad(h$fn, eta), tolerance = 1e-5)
})

test_that("the gradient matches numDeriv under ml, on the range space", {
  skip_if_not_installed("numDeriv")
  h <- crit_of_eta(y ~ s(x, k = 10), dg, ml(hessian = "observed"))
  eta <- h$eta0 + 0.3
  expect_equal(h$gr(eta), numDeriv::grad(h$fn, eta), tolerance = 1e-6)
})

test_that("the gradient vanishes at the reported optimum", {
  fit <- statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), dg,
                 outer_criterion = reml(hessian = "observed"))
  spec <- fit@spec
  design <- statmod_design(spec)
  idx <- outer_hyper_index(spec, statmod_blocks(spec, design))
  g <- statmod_marginal_grad(spec, design, fit@coefficients, fit@hyper,
                             reml(hessian = "observed"), idx)
  expect_lt(max(abs(g)), 1e-4)
})

test_that("the penalty is asked, not measured", {
  # the first version of this decided what a penalty was by testing whether
  # its Hessian happened to be linear in the hyperparameters, which excluded
  # every penalty built from a density. The question is now put to the
  # penalty, and a ridge answers it.
  q <- penalties7::quadratic_penalty(diag(4))
  expect_true(penalty_answers(q, 1L))
  expect_true(penalty_answers(q, 2L))
  expect_true(penalty_answers(penalties7::ridge_penalty(n_coef = 4L), 2L))
  # a kinked penalty has no such derivative and says so
  expect_false(penalty_answers(penalties7::lasso_penalty(n_coef = 4L), 1L))

  d <- penalties7::penalty_dhessian(q, rep(0.5, 4), list(lambda = 3))
  expect_equal(d$lambda, diag(4), tolerance = 1e-12)
  expect_equal(penalties7::penalty_hessian(q, rep(0.5, 4), list(lambda = 3)),
               3 * d$lambda, tolerance = 1e-12)
})

test_that("the exact route is taken only where it applies", {
  spec <- statmod_spec(y ~ s(x, k = 8), distributions7::gaussian1_distrib(),
                       dg)
  design <- statmod_design(spec)
  idx <- outer_hyper_index(spec, statmod_blocks(spec, design))
  expect_true(outer_gradient_ok(spec, design, idx, reml("observed")))
  # the expected information would need the derivative in beta of -E[l''],
  # which is not -E[l'''] and is not a generic of distributions7
  expect_false(outer_gradient_ok(spec, design, idx, reml("expected")))

  # a random effect's penalty is built from a density, and it is covered: its
  # derivatives come from penalties7's generics, which read the parent's
  # response surface rather than requiring the Hessian to be linear
  set.seed(33)
  dr <- data.frame(g = factor(rep(paste0("g", 1:12), each = 10)))
  dr$y <- stats::rnorm(120)
  sp2 <- statmod_spec(y ~ random(~ 1 | g),
                      distributions7::gaussian1_distrib(), dr)
  de2 <- statmod_design(sp2)
  ix2 <- outer_hyper_index(sp2, statmod_blocks(sp2, de2))
  expect_true(outer_gradient_ok(sp2, de2, ix2, reml("observed")))
  f <- statmod(y ~ random(~ 1 | g), distributions7::gaussian1_distrib(), dr,
               outer_criterion = reml("observed"))
  expect_true(is.finite(f@criterion))

  # a lasso is not, its penalty having a kink where a Laplace approximation
  # asks for a second derivative
  dl <- dr
  dl$n1 <- stats::rnorm(120)
  sp3 <- statmod_spec(y ~ lasso(~ n1),
                      distributions7::gaussian1_distrib(), dl)
  de3 <- statmod_design(sp3)
  ix3 <- outer_hyper_index(sp3, statmod_blocks(sp3, de3))
  expect_identical(nrow(ix3), 0L)
})

test_that("exact and derivative-free reach the same hyperparameter", {
  # the two routes share the criterion and nothing else, so agreeing on where
  # it is largest is a check of the gradient rather than of the optimizer
  fast <- statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), dg,
                  outer_criterion = reml(hessian = "observed"))
  slow <- statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), dg,
                  outer_criterion = reml(hessian = "observed"),
                  outer_optimizer = optimizers7::nelder_mead())
  lam_f <- fast@hyper$mu[["s(x, k = 10)"]][["lambda"]]
  lam_s <- slow@hyper$mu[["s(x, k = 10)"]][["lambda"]]
  expect_equal(lam_f, lam_s, tolerance = 1e-3)
  expect_equal(fast@criterion, slow@criterion, tolerance = 1e-6)
})

test_that("the gradient pays from the second hyperparameter on", {
  # ONE hyperparameter is where it does not pay, and the first version of this
  # test asserted the opposite from a guess. Measured, evaluations of the
  # criterion -- each a whole inner fit -- against the derivative-free search:
  #   1 parameter    40 against  32
  #   2 parameters   40 against 135
  #   3 parameters   41 against 269
  # a simplex needing one vertex per dimension and a quasi-Newton method not.
  skip_on_cran()
  set.seed(34)
  n2 <- 300
  d3 <- data.frame(a = runif(n2, -2, 2), b = runif(n2, -2, 2))
  d3$y <- sin(1.4 * d3$a) + d3$b^2 + stats::rnorm(n2, sd = 0.3)
  f <- y ~ s(a, k = 8) + s(b, k = 8)
  fast <- statmod(f, distributions7::gaussian1_distrib(), d3,
                  outer_criterion = reml(hessian = "observed"))
  slow <- statmod(f, distributions7::gaussian1_distrib(), d3,
                  outer_criterion = reml(hessian = "observed"),
                  outer_optimizer = optimizers7::nelder_mead())
  expect_equal(fast@criterion, slow@criterion, tolerance = 1e-5)
  expect_lt(nrow(fast@history$outer), nrow(slow@history$outer))
})
