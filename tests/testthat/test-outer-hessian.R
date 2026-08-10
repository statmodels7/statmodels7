# The exact Hessian of the marginal criterion.

set.seed(51)
n <- 250
dh <- data.frame(x = runif(n, -2, 2), z = runif(n))
dh$y <- sin(1.4 * dh$x) + stats::rnorm(n, sd = 0.3)

# The criterion, its exact gradient and its exact Hessian as functions of the
# free hyperparameters, the coefficients refitted at each. Differentiating the
# GRADIENT numerically is one stencil on an analytic quantity, which is the
# reference the toolkit sanctions; the value under test shares none of its
# arithmetic.
outer_handles <- function(formula, data, method, distrib = NULL) {
  if (is.null(distrib)) distrib <- distributions7::gaussian1_distrib()
  # The reference differentiates a gradient read at the penalized MODE, so
  # whatever the inner fit leaves short of stationarity is noise the numerical
  # derivative amplifies, and that sets how tightly the two can be compared.
  # Asking the inner fit for 1e-12 makes it WORSE rather than better: the
  # tolerance is unreachable, the run always spends its budget, and where it
  # stops then varies with the hyperparameter -- which is a rougher function
  # of it than the one the default produces. Measured, five of these
  # comparisons went from agreeing to disagreeing.
  inner <- iwls()
  fit0 <- statmod(formula, distrib, data, inner_method = inner)
  spec <- fit0@spec
  design <- statmod_design(spec)
  idx <- outer_hyper_index(spec, statmod_blocks(spec, design))
  basis <- integrated_basis(spec, design, method@kind)
  at <- function(eta) {
    hy <- eta_to_hyper(eta, idx, fit0@hyper)
    hl <- lapply(hy, function(pp)
      lapply(pp, function(v) stats::setNames(as.numeric(v), names(v))))
    hl <- hl[lengths(hl) > 0L]
    list(hy = hy, fit = statmod(formula, distrib, data, hyper = hl,
                                inner_method = inner))
  }
  list(
    idx = idx,
    eta0 = hyper_to_eta(fit0@hyper, idx),
    gr = function(eta) {
      a <- at(eta)
      statmod_marginal_grad(spec, design, a$fit@coefficients, a$hy, method,
                            idx, basis)
    },
    he = function(eta) {
      a <- at(eta)
      statmod_marginal_hess(spec, design, a$fit@coefficients, a$hy, method,
                            idx, basis)
    })
}

test_that("the Hessian of one smoothing parameter matches numDeriv", {
  skip_if_not_installed("numDeriv")
  h <- outer_handles(y ~ s(x, k = 10), dh, reml(hessian = "observed"))
  for (shift in c(0.4, -1.1)) {
    eta <- h$eta0 + shift
    expect_equal(as.numeric(h$he(eta)),
                 as.numeric(numDeriv::jacobian(h$gr, eta)), tolerance = 1e-5)
  }
})

test_that("the Hessian matches numDeriv with the scale modelled", {
  skip_if_not_installed("numDeriv")
  h <- outer_handles(y ~ s(x, k = 8) | sigma ~ z, dh,
                     reml(hessian = "observed"))
  eta <- h$eta0 + 0.25
  expect_equal(as.numeric(h$he(eta)),
               as.numeric(numDeriv::jacobian(h$gr, eta)), tolerance = 1e-5)
})

test_that("the Hessian matches numDeriv with two smooths", {
  skip_if_not_installed("numDeriv")
  set.seed(52)
  n2 <- 250
  d2 <- data.frame(a = runif(n2, -2, 2), b = runif(n2, -2, 2))
  d2$y <- sin(1.4 * d2$a) + d2$b^2 + stats::rnorm(n2, sd = 0.3)
  h <- outer_handles(y ~ s(a, k = 8) + s(b, k = 8), d2,
                     reml(hessian = "observed"))
  eta <- h$eta0 + c(0.5, -0.4)
  got <- h$he(eta)
  expect_identical(dim(got), c(2L, 2L))
  # a second derivative of a scalar is symmetric, and the assembly is written
  # so that it comes out symmetric rather than being symmetrized afterwards
  expect_equal(got, t(got), tolerance = 1e-12)
  expect_equal(as.numeric(got), as.numeric(numDeriv::jacobian(h$gr, eta)),
               tolerance = 1e-5)
})

test_that("the Hessian of an anisotropic tensor matches numDeriv", {
  skip_if_not_installed("numDeriv")
  # two hyperparameters inside ONE penalty, where the mixed entry comes from
  # the penalty's own second derivatives rather than from two separate terms.
  # The intercept is dropped: a tensor block contains the constant.
  set.seed(53)
  n2 <- 250
  d2 <- data.frame(a = runif(n2, -1, 1), b = runif(n2, -1, 1))
  d2$y <- d2$a^2 + sin(3 * d2$b) + stats::rnorm(n2, sd = 0.3)
  h <- outer_handles(y ~ te(a, b, k = 4) - 1, d2, reml(hessian = "observed"))
  eta <- h$eta0 + c(0.3, -0.3)
  # the loosest tolerance here, and it is the reference's: a tensor design is
  # the worst conditioned of these, so the mode moves least cleanly with the
  # hyperparameter and the difference of the gradient carries most noise
  expect_equal(as.numeric(h$he(eta)),
               as.numeric(numDeriv::jacobian(h$gr, eta)), tolerance = 1e-3)
})

test_that("the Hessian matches numDeriv under ml", {
  skip_if_not_installed("numDeriv")
  h <- outer_handles(y ~ s(x, k = 10), dh, ml(hessian = "observed"))
  eta <- h$eta0 + 0.3
  expect_equal(as.numeric(h$he(eta)),
               as.numeric(numDeriv::jacobian(h$gr, eta)), tolerance = 1e-5)
})

test_that("a variance component is covered, penalty and all", {
  # a ridge is NOT linear in its hyperparameter -- its Hessian carries
  # 1/sigma^2 -- and it works because the derivatives are asked of the penalty
  # rather than recovered from its behaviour
  skip_if_not_installed("numDeriv")
  set.seed(54)
  m <- 20
  dr <- data.frame(g = factor(rep(paste0("g", seq_len(m)), each = 10)),
                   x = runif(200))
  u <- stats::rnorm(m, sd = 0.7)
  dr$y <- 1 + 2 * dr$x + u[as.integer(dr$g)] + stats::rnorm(200, sd = 0.4)
  h <- outer_handles(y ~ x + random(~ 1 | g), dr, reml(hessian = "observed"))
  eta <- h$eta0 + 0.2
  expect_equal(h$gr(eta), numDeriv::grad(function(e) {
    # the criterion itself, differenced, as an independent check that the
    # gradient this Hessian is built on is the right one for THIS penalty
    a <- statmod(y ~ x + random(~ 1 | g),
                 distributions7::gaussian1_distrib(), dr,
                 hyper = list(mu = list(`random(~1 | g)` =
                   c(sigma = exp(e)))))
    statmod_marginal(a@spec, statmod_design(a@spec), a@coefficients,
                     a@hyper, reml(hessian = "observed"))$value
  }, eta), tolerance = 1e-5)
  expect_equal(as.numeric(h$he(eta)),
               as.numeric(numDeriv::jacobian(h$gr, eta)), tolerance = 1e-4)
})

test_that("a Newton step on the exact pair lands where the search does", {
  # a check that owes nothing to numDeriv: the gradient and the Hessian are
  # used to take Newton steps, and they must arrive where an optimizer that
  # only ever sees the value and the gradient arrives
  ref <- statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), dh,
                 outer_method = reml(hessian = "observed"))
  spec <- ref@spec
  design <- statmod_design(spec)
  idx <- outer_hyper_index(spec, statmod_blocks(spec, design))
  method <- reml(hessian = "observed")

  at <- function(eta) {
    hy <- eta_to_hyper(eta, idx, ref@hyper)
    f <- statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), dh,
                 hyper = list(mu = list(`s(x, k = 10)` =
                   c(lambda = exp(eta)))))
    list(g = statmod_marginal_grad(spec, design, f@coefficients, hy, method,
                                   idx),
         H = statmod_marginal_hess(spec, design, f@coefficients, hy, method,
                                   idx))
  }

  eta <- log(ref@hyper$mu[["s(x, k = 10)"]][["lambda"]]) + 1.5
  for (i in 1:8) {
    p <- at(eta)
    # a maximum, so the Hessian is negative there and the step is uphill
    eta <- eta - as.numeric(p$g / p$H)
  }
  # the two agree to 4e-5 relative, and Newton's point is the more stationary
  # of them: its gradient is below 1e-6 while the search stopped on its own
  # rule, so the tolerance here is the SEARCH's accuracy and not this pair's
  expect_equal(exp(eta), ref@hyper$mu[["s(x, k = 10)"]][["lambda"]],
               tolerance = 1e-4)
  expect_lt(abs(at(eta)$g), 1e-6)
  expect_lt(at(eta)$H, 0)
})

test_that("the exact route now covers the penalties it used to refuse", {
  spec <- statmod_spec(y ~ x + random(~ 1 | g),
                       distributions7::gaussian1_distrib(),
                       within(data.frame(x = runif(60)), {
                         g <- factor(rep(1:6, 10)); y <- stats::rnorm(60)
                       }))
  design <- statmod_design(spec)
  idx <- outer_hyper_index(spec, statmod_blocks(spec, design))
  # a ridge answers penalty_dhessian(), so the gradient and the Hessian are
  # both available where the old test of linearity said they were not
  expect_true(outer_gradient_ok(spec, design, idx, reml("observed"), 1L))
  expect_true(outer_gradient_ok(spec, design, idx, reml("observed"), 2L))
  # the expected information still is not, for the reason it never was
  expect_false(outer_gradient_ok(spec, design, idx, reml("expected"), 1L))
})
