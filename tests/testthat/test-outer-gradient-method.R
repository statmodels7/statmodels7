# The exact outer gradient and Hessian of a criterion on the expected
# information are the derivatives of the criterion that is evaluated, for every
# method of the expectation.
#
# The method object of expected_info() defines the matrix in the determinant
# and the derivatives that the gradient reads: the closed form with its
# closed-form derivative, the integrated expectation with the derivative of the
# expectation, and the outer product of the scores at the data with the exact
# derivative of that matrix. A gradient built on a derivative of another matrix
# than the one in the determinant is off by an amount that does not shrink with
# the step, so each check compares it with a central difference of the
# criterion value, the mode refitted and polished at every point.
#
# The step 1e-4 on the free scale of the hyperparameter leaves a truncation
# error of order 1e-9 in the difference, and the criterion (of order 1e2) is
# computed to about 1e-12, which divided by the step is a rounding error of
# 1e-8 at worst. The comparisons use 1e-7 for the gradient and 1e-6 for the
# Hessian, against measured agreement of 2e-10 and 1e-8.

crit_of_eta_m <- function(formula, data, distrib, method) {
  fit0 <- statmod(formula, distrib, data, outer_criterion = NULL)
  spec <- fit0@spec
  design <- statmod_design(spec)
  blocks <- statmod_blocks(spec, design)
  idx <- outer_hyper_index(spec, blocks)
  basis <- integrated_basis(spec, design, method@kind)
  emethod <- info_emethod(method@hessian)
  at <- function(eta) {
    hy <- eta_to_hyper(eta, idx, fit0@hyper)
    list(hy = hy, cf = fit_at_hyper(formula, distrib, data, hy,
                                    iwls(tol = 1e-9),
                                    polish = TRUE)$coefficients)
  }
  ctx_at <- function(a) outer_context(spec, design, a$cf, a$hy, emethod)
  list(
    eta0 = hyper_to_eta(fit0@hyper, idx),
    gate1 = outer_gradient_ok(spec, design, idx, method, 1L),
    gate2 = outer_gradient_ok(spec, design, idx, method, 2L),
    fn = function(eta) {
      a <- at(eta)
      statmod_marginal(spec, design, a$cf, a$hy, method, emethod, basis,
                       ctx_at(a))$value
    },
    gr = function(eta) {
      a <- at(eta)
      statmod_marginal_grad(spec, design, a$cf, a$hy, method, idx, basis,
                            ctx = ctx_at(a))
    },
    he = function(eta) {
      a <- at(eta)
      statmod_marginal_hess(spec, design, a$cf, a$hy, method, idx, basis,
                            ctx = ctx_at(a))
    })
}

fd_grad <- function(f, eta, h = 1e-4) {
  vapply(seq_along(eta), function(i) {
    e <- numeric(length(eta)); e[i] <- h
    (f(eta + e) - f(eta - e)) / (2 * h)
  }, numeric(1))
}

fd_hess <- function(g, eta, h = 1e-4) {
  vapply(seq_along(eta), function(i) {
    e <- numeric(length(eta)); e[i] <- h
    (g(eta + e) - g(eta - e)) / (2 * h)
  }, numeric(length(eta)))
}

pig_small <- function(n = 40L, seed = 11L) {
  set.seed(seed)
  dd <- data.frame(x = runif(n))
  dd$y <- distributions7::distrib_rng(distributions7::pig1_distrib(), n,
    list(mu = exp(0.8 + 0.5 * sin(2 * pi * dd$x)), sigma = 0.5))
  dd
}

check_consistency <- function(h, optimum) {
  expect_true(h$gate1)
  expect_true(h$gate2)
  for (eta in list(h$eta0 + 0.4, optimum)) {
    g <- h$gr(eta)
    expect_lt(max(abs(g - fd_grad(h$fn, eta))), 1e-7)
    H <- h$he(eta)
    expect_lt(max(abs(H - fd_hess(h$gr, eta))), 1e-6)
  }
}

test_that("the default method of a density-only family is the outer product of the scores, and the derivatives follow", {
  skip_if_not_installed("numDeriv")
  dd <- pig_small()
  fo <- y ~ s(x, bspline_smooth(k = 6))
  m <- reml(information = expected_info())
  h <- crit_of_eta_m(fo, dd, pig_bare_distrib(), m)
  fit <- statmod(fo, pig_bare_distrib(), dd, outer_criterion = m)
  idx <- outer_hyper_index(fit@spec, statmod_blocks(fit@spec, statmod_design(fit@spec)))
  check_consistency(h, hyper_to_eta(fit@hyper, idx))
  # naming the method is the same criterion as the default
  h2 <- crit_of_eta_m(fo, dd, pig_bare_distrib(),
                      reml(information = expected_info(
                        method = distributions7::expect_opg())))
  eta <- h$eta0 + 0.4
  expect_identical(h$fn(eta), h2$fn(eta))
  expect_identical(h$gr(eta), h2$gr(eta))
  # and ml() reads the same matrix on the range space
  hm <- crit_of_eta_m(fo, dd, pig_bare_distrib(),
                      ml(information = expected_info()))
  expect_lt(max(abs(hm$gr(eta) - fd_grad(hm$fn, eta))), 1e-7)
})

test_that("an explicit outer product of the scores is used also by a family with a closed form", {
  dd <- data.frame(x = seq(0, 1, length.out = 60))
  set.seed(5)
  dd$y <- stats::rgamma(60, shape = 4, rate = 4 / exp(0.4 * sin(3 * dd$x) + 1))
  fo <- y ~ s(x, bspline_smooth(k = 6))
  g <- distributions7::gamma1_distrib()
  hc <- crit_of_eta_m(fo, dd, g, reml(information = expected_info()))
  ho <- crit_of_eta_m(fo, dd, g, reml(information = expected_info(
    method = distributions7::expect_opg())))
  eta <- hc$eta0 + 0.4
  # the matrices differ, so the criteria do
  expect_false(isTRUE(all.equal(hc$fn(eta), ho$fn(eta), tolerance = 1e-8)))
  # and each gradient is the derivative of its own criterion
  expect_lt(max(abs(hc$gr(eta) - fd_grad(hc$fn, eta))), 1e-7)
  expect_lt(max(abs(ho$gr(eta) - fd_grad(ho$fn, eta))), 1e-7)
  expect_lt(max(abs(ho$he(eta) - fd_hess(ho$gr, eta))), 1e-6)
})

test_that("the integrated expectation of the score identity has consistent derivatives", {
  skip_on_cran()
  dd <- pig_small(40L)
  fo <- y ~ s(x, bspline_smooth(k = 6))
  m <- reml(information = expected_info(
    method = distributions7::expect_integrate(identity = "score")))
  h <- crit_of_eta_m(fo, dd, pig_bare_distrib(), m)
  expect_true(h$gate1)
  eta <- h$eta0 + 0.4
  expect_lt(max(abs(h$gr(eta) - fd_grad(h$fn, eta))), 1e-7)
  expect_lt(max(abs(h$he(eta) - fd_hess(h$gr, eta))), 1e-6)
})

test_that("a Monte Carlo expectation has a gradient that estimates the one of the expectation", {
  skip_on_cran()
  dd <- pig_small(40L)
  fo <- y ~ s(x, bspline_smooth(k = 6))
  mi <- reml(information = expected_info(
    method = distributions7::expect_integrate(identity = "score")))
  mm <- reml(information = expected_info(
    method = distributions7::expect_mc(nsim = 4000, identity = "score")))
  hi <- crit_of_eta_m(fo, dd, pig_bare_distrib(), mi)
  hm <- crit_of_eta_m(fo, dd, pig_bare_distrib(), mm)
  expect_true(hm$gate1)
  eta <- hi$eta0 + 0.4
  gi <- hi$gr(eta)
  set.seed(1)
  gm <- hm$gr(eta)
  # a sample is drawn at every call of the matrix and of each derivative, so
  # the gradient is that of the expectation up to the Monte Carlo error, and
  # not that of the sampled matrix
  expect_lt(max(abs(gm - gi)) / max(abs(gi)), 0.05)
})
