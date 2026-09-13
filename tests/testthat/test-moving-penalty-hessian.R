# The outer Hessian beside a penalty whose Hessian moves with the coefficients.
#
# A heavy-tailed prior on a random effect has S = -D' diag(l_yy(D beta)) D,
# which depends on the effects, so differentiating the criterion's gradient
# once more reads how dS/dbeta[v] moves: along a second direction
# (penalties7::penalty_d2hessian_beta()) and in each hyperparameter
# (penalties7::penalty_dhessian_beta_theta()). Without those pieces the
# assembled Hessian was 122 per cent out and flat in the step.
#
# The reference is a central difference of the EXACT GRADIENT, validated in
# test-outer-gradient.R, with the mode refitted from the same start at every
# probe. Its own floor is the mode's reproducibility, about 2e-4 relative on
# this panel, so what is asserted is the rate above that floor and a bound at
# the step where the floor has not yet been reached.

mp_panel <- function(seed = 1L, m = 30L, ni = 12L) {
  set.seed(seed)
  b <- stats::rnorm(m, 0, 0.3)
  b[1:3] <- c(2.5, -2.8, 3.1)
  g <- factor(rep(seq_len(m), each = ni))
  x <- stats::rnorm(m * ni)
  data.frame(g = g, x = x,
             y = 1 + 0.5 * x + b[as.integer(g)] + stats::rnorm(m * ni))
}

mp_formula <- y ~ x + random(~ 1 | g, distrib = distributions7::fixed(
  distributions7::student_t1_distrib(), mu = 0))

mp_parts <- function(oc) {
  d <- mp_panel()
  fit <- statmod(mp_formula, distributions7::gaussian1_distrib(), d,
                 outer_criterion = oc)
  spec <- fit@spec
  design <- statmod_design(spec)
  blocks <- statmod_blocks(spec, design)
  method <- fit@methods$outer
  pe <- outer_minimize(method)
  list(fit = fit, spec = spec, design = design, blocks = blocks,
       idx = outer_hyper_index(spec, blocks), method = method, pe = pe,
       basis = if (pe) NULL else integrated_basis(spec, design, method@kind))
}

mp_derivs <- function(p, coef, hy, order) {
  if (p$pe) {
    statmod_pe_derivs(p$spec, p$design, coef, hy, p$method, p$idx, order)
  } else {
    g <- statmod_marginal_grad(p$spec, p$design, coef, hy, p$method, p$idx,
                               p$basis)
    if (order < 2L) list(grad = g) else
      list(grad = g, hess = statmod_marginal_hess(p$spec, p$design, coef, hy,
                                                  p$method, p$idx, p$basis))
  }
}

mp_check <- function(p) {
  inner <- iwls()
  cfg <- inner_settings(inner, p$spec@distrib)
  beta0 <- unlist(p$fit@coefficients[p$spec@distrib@params], use.names = FALSE)
  at <- function(eta) {
    hy <- eta_to_hyper(eta, p$idx, p$fit@hyper)
    r <- statmod_alternate(p$spec, p$design, p$blocks, hy, inner, beta0,
                           cfg$expected, cfg$approx, cfg$maxit, 1e-10,
                           verbosity(0))
    list(coef = r$obj$split(r$par), hy = hy)
  }
  eta0 <- hyper_to_eta(p$fit@hyper, p$idx)
  r0 <- at(eta0)
  H <- mp_derivs(p, r0$coef, r0$hy, 2L)$hess
  nh <- nrow(p$idx)
  gaps <- vapply(c(8e-2, 4e-2), function(h) {
    Hn <- matrix(0, nh, nh)
    for (m in seq_len(nh)) {
      ep <- eta0; ep[[m]] <- ep[[m]] + h
      em <- eta0; em[[m]] <- em[[m]] - h
      Hn[, m] <- (mp_derivs(p, at(ep)$coef, at(ep)$hy, 1L)$grad -
                    mp_derivs(p, at(em)$coef, at(em)$hy, 1L)$grad) / (2 * h)
    }
    max(abs(H - (Hn + t(Hn)) / 2)) / max(abs(H))
  }, numeric(1))
  list(H = H, gaps = gaps)
}

test_that("order 2 is answered beside a heavy-tailed prior", {
  p <- mp_parts(reml())
  expect_true(outer_gradient_ok(p$spec, p$design, p$idx, p$method, 1L))
  expect_true(outer_gradient_ok(p$spec, p$design, p$idx, p$method, 2L))
  # the premise: the prior's Hessian really moves with the effects here
  u <- Filter(function(z) !is.null(z$penalty),
              statmod_penalized(p$spec, p$design))[[1L]]
  th <- as.list(p$fit@hyper[[u$param]][[u$key]])
  expect_false(isTRUE(penalties7::beta_quadratic(u$penalty, th)))
})

test_that("the REML Hessian converges on a difference of the exact gradient", {
  skip_on_cran()
  p <- mp_parts(reml())
  r <- mp_check(p)
  expect_identical(r$H, t(r$H))
  # the RATE: 3.7e-03 and 5.6e-04 measured, O(h^2), where the assembly
  # without the penalty's second movement is 1.22 at every step
  expect_gt(r$gaps[[1L]] / r$gaps[[2L]], 4)
  expect_lt(r$gaps[[2L]], 2e-3)
  # and what a reader reads now exists: a verdict on an analytic curvature
  # and a variance for the prior's hyperparameters
  cert <- statmod_certificate(p$fit)
  expect_identical(cert$curvature, "analytic")
  expect_false(is.null(statmod_hyper_vcov(p$spec, p$design, p$fit@coefficients,
                                          p$fit@hyper, p$method)))
})

test_that("the prediction-error Hessian converges on its exact gradient", {
  skip_on_cran()
  p <- mp_parts(aic())
  expect_true(outer_gradient_ok(p$spec, p$design, p$idx, p$method, 2L))
  r <- mp_check(p)
  # 1.6e-03 and 3.9e-04 measured
  expect_gt(r$gaps[[1L]] / r$gaps[[2L]], 3)
  expect_lt(r$gaps[[2L]], 2e-3)
})

test_that("beside a structural term a moving penalty keeps the stencil", {
  # statmod_structural_hess() reads no movement of the penalty's Hessian, so
  # order 2 is refused there rather than answered with a piece missing
  set.seed(9)
  m <- 6L; ni <- 20L
  d <- data.frame(g = factor(rep(seq_len(m), each = ni)))
  d$y <- stats::rnorm(m * ni) + stats::rnorm(m, 0, 0.5)[as.integer(d$g)]
  spec <- statmod_spec(y ~ random(~ 1 | g, distrib = distributions7::fixed(
    distributions7::student_t1_distrib(), mu = 0)) +
      gas(p = 1, q = 1, by = g, alpha1 ~ 1 + random(~ 1 | g)),
    distributions7::gaussian1_distrib(), d)
  design <- statmod_design(spec)
  idx <- outer_hyper_index(spec, statmod_blocks(spec, design))
  expect_false(outer_gradient_ok(spec, design, idx, reml(), 2L))
})
