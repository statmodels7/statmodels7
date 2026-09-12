# A structural term no penalty covers is estimated beside the coefficients and
# moves with the hyperparameter all the same. The criterion's determinant is
# over the coefficients alone, on the static design, but the predictors it is
# read at carry the term's contribution, so its exact gradient reads the mode's
# JOINT movement and the information's movement along the TOTAL derivative of
# every predictor. For a filter that is written; for a term of the likelihood
# shape it would need the derivative of the posterior along a direction, which
# modelterms7 does not expose, and the gradient is refused there.

# a series whose level is driven by the model's OWN score, beside a smooth
sim_unpen_filter <- function(seed = 11, n = 300) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(n, -2, 2), t = seq_len(n))
  f <- 0
  y <- numeric(n)
  for (i in seq_len(n)) {
    mu <- sin(2 * d$x[i]) + f
    y[i] <- stats::rnorm(1, mu, 0.5)
    f <- 0.25 * (y[i] - mu) / 0.25 + 0.7 * f
  }
  d$y <- y
  d
}

unpen_form <- y ~ s(x, bspline_smooth(k = 8)) + gas(p = 1, q = 1, time = t)

# the criterion and its exact gradient with the mode refitted from a FIXED
# start at every point, read at the refitted design: a structural term's own
# parameters live in the design's state
unpen_harness <- function(form, data, method) {
  fit0 <- statmod(form, distributions7::gaussian1_distrib(), data,
                  outer_criterion = NULL)
  idx <- outer_hyper_index(fit0@spec,
                           statmod_blocks(fit0@spec,
                                          statmod_design(fit0@spec)))
  marginal <- method@kind %in% c("reml", "ml")
  refit <- function(eta) {
    hy <- eta_to_hyper(eta, idx, fit0@hyper)
    a <- fit_at_hyper(form, distributions7::gaussian1_distrib(), data, hy)
    a$hy <- hy
    a$basis <- integrated_basis(a$spec, a$design, method@kind)
    a
  }
  list(idx = idx, eta0 = hyper_to_eta(fit0@hyper, idx),
       at = function(eta) refit(eta),
       fn = function(eta) {
         a <- refit(eta)
         if (marginal) {
           statmod_marginal(a$spec, a$design, a$coefficients, a$hy, method,
                            basis = a$basis)$value
         } else {
           statmod_pe(a$spec, a$design, a$coefficients, a$hy, method)$value
         }
       },
       gr = function(eta) {
         a <- refit(eta)
         if (marginal) {
           statmod_marginal_grad(a$spec, a$design, a$coefficients, a$hy,
                                 method, idx, a$basis)
         } else {
           statmod_pe_derivs(a$spec, a$design, a$coefficients, a$hy, method,
                             idx)$grad
         }
       })
}

central_diff <- function(fn, eta, h) {
  vapply(seq_along(eta), function(j) {
    ep <- eta; em <- eta
    ep[j] <- ep[j] + h; em[j] <- em[j] - h
    (fn(ep) - fn(em)) / (2 * h)
  }, numeric(1))
}

# THE SIGNATURE, not a tolerance alone: a gradient missing a term is out by a
# fixed relative amount whatever the step, while a correct one converges
# O(h^2) onto the reference's floor. Measured on a smooth beside an
# unpenalized gas(1, 1) at 400 observations, 7.2e-03 flat before under
# reml(), and 5.6e-05, 5.3e-06, 1.2e-06 at h of 1e-2, 3e-3, 1e-3 after.
expect_converges <- function(h, eta, bound) {
  g <- h$gr(eta)
  rel <- vapply(c(1e-2, 1e-3), function(hh) {
    fd <- central_diff(h$fn, eta, hh)
    max(abs(g - fd) / abs(fd))
  }, 0)
  expect_lt(rel[[2L]], rel[[1L]] / 10)
  expect_lt(rel[[2L]], bound)
}

test_that("the gradient beside an unpenalized filter reads the joint movement", {
  skip_on_cran()
  d <- sim_unpen_filter()
  h <- unpen_harness(unpen_form, d, reml())
  a <- h$at(h$eta0)
  expect_false(structural_penalized(a$spec, a$design))
  expect_true(outer_gradient_ok(a$spec, a$design, h$idx, reml(), 1L))
  expect_converges(h, h$eta0 + 0.3, 1e-4)
})

test_that("a prediction-error criterion beside a filter reads it too", {
  # the same two pieces, which tau reads through the information alone. Read
  # on the coefficients the gradient was out by 9.48 under aic() and 2.22
  # under bic(), relative and flat in the step
  skip_on_cran()
  d <- sim_unpen_filter()
  for (crit in list(aic(), bic())) {
    h <- unpen_harness(unpen_form, d, crit)
    expect_converges(h, h$eta0 + 0.3, 1e-3)
  }
})

test_that("an unpenalized filter's order 2 is the stencil of the exact gradient", {
  # the joint assembly is the second derivative of the criterion whose
  # determinant spans the term's parameters, which is not this one: measured,
  # -3.97422 against a second difference of -4.00272, 0.71 per cent out
  skip_on_cran()
  d <- sim_unpen_filter()
  fit <- statmod(unpen_form, distributions7::gaussian1_distrib(), d,
                 outer_criterion = reml())
  des <- statmod_design(fit@spec)
  idx <- outer_hyper_index(fit@spec, statmod_blocks(fit@spec, des))
  expect_false(outer_gradient_ok(fit@spec, des, idx, reml(), 2L))
  expect_false(outer_gradient_ok(fit@spec, des, idx, aic(), 2L))
  H <- statmod_marginal_hess(fit@spec, des, fit@coefficients, fit@hyper,
                             fit@methods$outer, idx)
  S <- statmod_hess_stencil(fit@spec, des, fit@coefficients, fit@hyper,
                            fit@methods$outer, idx, NULL, iwls())
  expect_identical(H, S)

  h <- unpen_harness(unpen_form, d, reml())
  e <- hyper_to_eta(fit@hyper, idx)
  hh <- 1e-2
  d2 <- (h$fn(e + hh) - 2 * h$fn(e) + h$fn(e - hh)) / hh^2
  expect_equal(as.numeric(H), d2, tolerance = 1e-3)
  # and the prediction-error route's own assembly says so rather than
  # answering over the coefficients
  expect_error(statmod_pe_derivs(fit@spec, des, fit@coefficients, fit@hyper,
                                 aic(), idx, 2L), "no exact Hessian")
})

test_that("a term of the likelihood shape is refused and says what it costs", {
  # regime()'s exact gradient needs how the smoothed posterior moves along the
  # direction the mode moves in, which the term does not supply. Read on the
  # coefficients alone it was out by 2.05e-04, flat in the step, so it is
  # refused: the search is derivative-free and the certificate unknown
  skip_on_cran()
  set.seed(5)
  n <- 300L
  Z <- matrix(stats::rnorm(n * 4), n, 4)
  st <- integer(n)
  st[1L] <- 1L
  for (i in 2:n) st[i] <- if (stats::runif(1) < 0.05) 3L - st[i - 1L] else st[i - 1L]
  dr <- data.frame(y = as.numeric(Z %*% c(0.9, -0.6, 0, 0)) +
                     ifelse(st == 1L, -1.5, 1.5) + stats::rnorm(n, sd = 0.8))
  dr$Z <- Z
  fit <- statmod(y ~ ridge(Z) + regime(k = 2),
                 distributions7::gaussian1_distrib(), dr,
                 outer_criterion = reml())
  des <- statmod_design(fit@spec)
  idx <- outer_hyper_index(fit@spec, statmod_blocks(fit@spec, des))
  expect_false(outer_gradient_ok(fit@spec, des, idx, reml(), 1L))
  expect_false(outer_gradient_ok(fit@spec, des, idx, reml("expected"), 1L))
  expect_identical(class(fit@methods$search)[[1L]], "optimizers7::NelderMead")
  ce <- statmod_certificate(fit)
  expect_identical(ce$state, "unknown")
  expect_match(ce$reason, "no exact outer gradient")
})
