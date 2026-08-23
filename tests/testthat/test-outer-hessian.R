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
  # Asking the inner fit for more is NOT the remedy: the stall guard on the
  # objective (a decrease under 1e-12 of its own magnitude is rounding, not
  # progress) fires before a tighter score rule can, so a tighter tolerance
  # only makes where the run stops vary with the hyperparameter -- measured
  # twice, at 1e-12 under the absolute rule and again at 1e-8 under a
  # dimensionless one.
  inner <- iwls()
  fit0 <- statmod(formula, distrib, data, inner_optimizer = inner)
  spec <- fit0@spec
  design <- statmod_design(spec)
  idx <- outer_hyper_index(spec, statmod_blocks(spec, design))
  basis <- integrated_basis(spec, design, method@kind)
  at <- function(eta) {
    hy <- eta_to_hyper(eta, idx, fit0@hyper)
    list(hy = hy, fit = fit_at_hyper(formula, distrib, data, hy, inner))
  }
  list(
    idx = idx,
    eta0 = hyper_to_eta(fit0@hyper, idx),
    gr = function(eta) {
      a <- at(eta)
      statmod_marginal_grad(spec, design, a$fit$coefficients, a$hy, method,
                            idx, basis)
    },
    he = function(eta) {
      a <- at(eta)
      statmod_marginal_hess(spec, design, a$fit$coefficients, a$hy, method,
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
  set.seed(53)
  n2 <- 250
  d2 <- data.frame(a = runif(n2, -1, 1), b = runif(n2, -1, 1))
  d2$y <- d2$a^2 + sin(3 * d2$b) + stats::rnorm(n2, sd = 0.3)
  h <- outer_handles(y ~ te(a, b, k = 4), d2, reml(hessian = "observed"))
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
  # the effects' penalty is a ridge, and it works because the derivatives
  # are asked of the penalty rather than recovered from its behaviour
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
    a <- statmod(y ~ x + random(~ 1 | g, hyper = c(sigma = exp(e))),
                 distributions7::gaussian1_distrib(), dr)
    statmod_marginal(a@spec, statmod_design(a@spec), a@coefficients,
                     a@hyper, reml(hessian = "observed"))$value
    # the reference REFITS the mode at every probe, so its own accuracy is
    # what the tolerance has to allow: measured, the two agree to 1.1e-5
    # relative, and the defect this exists for -- a derivative read on the
    # wrong chart -- is off by a factor, not by a rounding
  }, eta), tolerance = 1e-4)
  expect_equal(as.numeric(h$he(eta)),
               as.numeric(numDeriv::jacobian(h$gr, eta)), tolerance = 1e-4)
})

test_that("a Newton step on the exact pair lands where the search does", {
  # a check that owes nothing to numDeriv: the gradient and the Hessian are
  # used to take Newton steps, and they must arrive where an optimizer that
  # only ever sees the value and the gradient arrives
  ref <- statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), dh,
                 outer_criterion = reml(hessian = "observed"))
  spec <- ref@spec
  design <- statmod_design(spec)
  idx <- outer_hyper_index(spec, statmod_blocks(spec, design))
  method <- reml(hessian = "observed")

  at <- function(eta) {
    hy <- eta_to_hyper(eta, idx, ref@hyper)
    f <- statmod(y ~ s(x, k = 10, lambda = exp(eta)),
                 distributions7::gaussian1_distrib(), dh)
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
  # the expected information carries the GRADIENT now, through
  # distrib_dexpected_hessian(), and not the Hessian: order 2 would ask for
  # the next order of dE[l'']/deta
  expect_true(outer_gradient_ok(spec, design, idx, reml("expected"), 1L))
  expect_false(outer_gradient_ok(spec, design, idx, reml("expected"), 2L))
})

test_that("a variance component's Hessian is exact now, not differenced", {
  # penalties7 used to difference the second-order pieces of a separable
  # penalty; distributions7 supplies them closed for a gaussian parent, so a
  # ridge and a random effect go through no difference at all. The tolerance
  # is what says so: it is the reference's, not ours.
  skip_if_not_installed("numDeriv")
  set.seed(55)
  m <- 20
  dr <- data.frame(g = factor(rep(paste0("g", seq_len(m)), each = 10)),
                   x = runif(200))
  u <- stats::rnorm(m, sd = 0.7)
  dr$y <- 1 + 2 * dr$x + u[as.integer(dr$g)] + stats::rnorm(200, sd = 0.4)
  h <- outer_handles(y ~ x + random(~ 1 | g), dr, reml(hessian = "observed"))
  for (shift in c(0.2, -0.4)) {
    eta <- h$eta0 + shift
    expect_equal(as.numeric(h$he(eta)),
                 as.numeric(numDeriv::jacobian(h$gr, eta)), tolerance = 1e-3)
  }
  # and the fit runs through newton(), which needs the Hessian at every step
  f <- statmod(y ~ x + random(~ 1 | g),
               distributions7::gaussian1_distrib(), dr,
               outer_criterion = reml(hessian = "observed"))
  expect_true(is.finite(f@criterion))
  expect_gt(f@hyper$mu[["random(~1 | g)"]][[1L]], 0.3)
})


# --- a block that moves with its coefficients -------------------------------

# The pieces the refresh corrections are assembled from, at one point, so that
# each can be asked its own question rather than only the Hessian's.
refresh_bits <- function(formula, data, shift = 0.3) {
  inner <- iwls()
  fit0 <- statmod(formula, distributions7::gaussian1_distrib(), data,
                  inner_optimizer = inner,
                  outer_criterion = reml(hessian = "observed"))
  spec <- fit0@spec
  design0 <- statmod_design(spec)
  idx <- outer_hyper_index(spec, statmod_blocks(spec, design0))
  hy <- eta_to_hyper(hyper_to_eta(fit0@hyper, idx) + shift, idx, fit0@hyper)
  cf <- fit_at_hyper(formula, distributions7::gaussian1_distrib(), data, hy,
                     inner)$coefficients
  design <- statmod_design_at(spec, cf, design0)
  params <- spec@distrib@params
  npar <- vapply(design, function(z) z$npar, integer(1))
  offs <- cumsum(npar) - npar
  total <- sum(npar)
  pen <- ctx_penalized(NULL, spec, design, cf, hy, FALSE)
  M <- ctx_trace_matrix(NULL, pen, integrated_basis(spec, design, "reml"),
                        FALSE)
  list(spec = spec, design = design, coef = cf, hyper = hy, M = M,
       params = params, npar = npar, offs = offs, total = total,
       d3 = ctx_deriv(NULL, spec, design, cf, hy, 3L),
       Hl = refresh_hessian(spec, design, cf, FALSE, "bartlett"),
       units = refresh_units(spec, design, cf, params, npar, offs))
}

test_that("the Hessian covers a block that moves with its coefficients", {
  skip_if_not_installed("numDeriv")
  # nl()'s block is the Jacobian, so dX/dbeta reaches the assembly in three
  # places -- the matrix dK/dt, the trace of dK_m/dt_l against M, and the
  # twice-contracted fourth derivative -- and the mode moves by the penalized
  # likelihood's own curvature rather than by K. Beside those, the twice-
  # contracted fourth derivative reads the block's SECOND derivative and the
  # direction's own predictor, and the mode's SECOND movement reads both:
  # measured against a central difference of the exact gradient with the mode
  # refitted, the whole set takes this cell from 2.19e-05 to 1.63e-07, which
  # is the reference's own floor, and a nearly straight one (r*x_max = 1.5,
  # where the criterion is flat and the Hessian small) from 5.48e-04 to
  # 1.37e-07 with its hyperparameter's standard error from 2.74e-04 to
  # 6.85e-08.
  set.seed(9)
  np <- 40; m <- 8
  dn <- data.frame(x = rep(seq(0.2, 4, length.out = np), m),
                   grp = factor(rep(sprintf("g%d", seq_len(m)), each = np)))
  a_g <- 3 + stats::rnorm(m, sd = 0.4)
  dn$y <- a_g[dn$grp] * exp(-0.6 * dn$x) + stats::rnorm(nrow(dn), sd = 0.15)
  h <- outer_handles(y ~ nl(~ a * exp(-r * x), a ~ 0 + ridge(~ grp)), dn,
                     reml(hessian = "observed"))
  eta <- h$eta0 + 0.35
  expect_equal(as.numeric(h$he(eta)),
               as.numeric(numDeriv::jacobian(h$gr, eta)), tolerance = 1e-3)
})

test_that("a break-point term's Hessian is covered too", {
  skip_if_not_installed("numDeriv")
  # seg()'s block is the Jacobian of a CONTINUOUS construction, so it carries
  # term_block_deriv() and the same four corrections apply; jump() and jseg()
  # answer zeros, their position being read off a product of coefficients.
  set.seed(12)
  ns <- 300
  ds <- data.frame(x = runif(ns, 0, 10),
                   id = factor(rep(seq_len(6), length.out = ns)))
  ds$y <- 1 + 0.3 * ds$x + 1.5 * pmax(ds$x - 5, 0) + stats::rnorm(ns, sd = 0.4)
  h <- outer_handles(y ~ seg(x, gamma1 ~ 0 + ridge(~ id)), ds,
                     reml(hessian = "observed"))
  eta <- h$eta0 + 0.3
  expect_equal(as.numeric(h$he(eta)),
               as.numeric(numDeriv::jacobian(h$gr, eta)), tolerance = 1e-4)
})

test_that("a fixed design gets exactly zero from every refresh correction", {
  # the negative control, and it is what says the corrections cannot move a
  # model with no block that moves: they are not small there, they are the
  # zero matrix and the number zero.
  b <- refresh_bits(y ~ s(x, k = 8), dh)
  expect_length(b$units, 0L)
  R <- contract3_refresh(b$spec, b$design, b$params, b$npar, b$offs, b$total,
                         list(), b$Hl, b$units)
  expect_true(all(as.matrix(R) == 0))
  expect_identical(trace_refresh4(b$spec, b$M, b$params, b$npar, b$Hl, list(),
                                  list(), b$units, NULL, b$d3, NULL, list(),
                                  list()), 0)
  expect_true(all(u_refresh(b$spec, b$design, b$coef, b$M, b$params, b$npar,
                            b$offs, b$total, units = b$units) == 0))
  expect_true(all(refresh_mode_third(b$spec, b$params, b$npar, b$units, b$Hl,
                                     NULL, NULL, list(), list(),
                                     b$total) == 0))
})

test_that("the trace of the refresh correction is read off the adjoint", {
  # tr(M (R + R')) is v'u with u the contraction the GRADIENT already forms,
  # because term_block_contract() is the adjoint of term_block_deriv(). The
  # Hessian uses that shortcut for one of its three places; here it is checked
  # against the matrix assembled and traced, which shares none of its
  # arithmetic. A term whose block does not move would satisfy this with zeros
  # on both sides, so the nl term is asked and its own correction is required
  # to be non-trivial.
  set.seed(9)
  np <- 30; m <- 6
  dn <- data.frame(x = rep(seq(0.2, 4, length.out = np), m),
                   grp = factor(rep(sprintf("g%d", seq_len(m)), each = np)))
  a_g <- 3 + stats::rnorm(m, sd = 0.4)
  dn$y <- a_g[dn$grp] * exp(-0.6 * dn$x) + stats::rnorm(nrow(dn), sd = 0.15)
  b <- refresh_bits(y ~ nl(~ a * exp(-r * x), a ~ 0 + ridge(~ grp)), dn)
  expect_length(b$units, 1L)
  set.seed(4)
  v <- stats::rnorm(b$total)
  uref <- u_refresh(b$spec, b$design, b$coef, b$M, b$params, b$npar, b$offs,
                    b$total, units = b$units, Hl = b$Hl)
  dir <- refresh_direction(b$spec, b$design, b$M, b$params, b$npar, b$offs,
                           b$d3,
                           block_predictors(b$design, b$params, b$npar,
                                            b$offs, v),
                           v, b$units)
  R <- contract3_refresh(b$spec, b$design, b$params, b$npar, b$offs, b$total,
                         dir, b$Hl, b$units)
  direct <- sum(as.matrix(b$M) * as.matrix(R + t(R)))
  expect_gt(abs(direct), 1e-8)
  expect_equal(sum(uref * v), direct, tolerance = 1e-10)
})


# The twice-contracted fourth derivative on its own, against a mixed second
# difference of H. It is a far sharper instrument than the Hessian end to end,
# where one addend's error is diluted among the others, and it is what the
# corrections reading the block's second derivative were written against.
u_trace_bits <- function(formula, data, shift = 0.3) {
  b <- refresh_bits(formula, data, shift)
  spec <- b$spec
  design0 <- statmod_design(spec)
  join <- function(cf) unlist(cf[spec@distrib@params], use.names = FALSE)
  splt <- function(bv) {
    out <- list()
    for (a in seq_along(b$params)) {
      out[[b$params[a]]] <- if (b$npar[a] == 0L) numeric(0) else
        bv[b$offs[a] + seq_len(b$npar[a])]
    }
    out
  }
  # H alone: the order-2 route assumes the penalty is quadratic in beta, so
  # d2S/dbeta2 is zero and the penalty does not enter U
  Hof <- function(bv) {
    cf <- splt(bv)
    as_dense(statmod_information_at(spec, cf,
                                    statmod_design_at(spec, cf, design0),
                                    expected = FALSE))
  }
  c(b, list(bflat = join(b$coef), Hof = Hof))
}

u_trace_gap <- function(b, v, u, h = 3e-4) {
  M <- as_dense(b$M)
  ref <- (b$Hof(b$bflat + h * v + h * u) - b$Hof(b$bflat + h * v - h * u) -
          b$Hof(b$bflat - h * v + h * u) + b$Hof(b$bflat - h * v - h * u)) /
    (4 * h^2)
  d4 <- ctx_deriv(NULL, b$spec, b$design, b$coef, b$hyper, 4L)
  G <- block_leverage(b$design, M, b$params, b$npar, b$offs, 1L)
  tv <- block_predictors(b$design, b$params, b$npar, b$offs, v)
  tu <- block_predictors(b$design, b$params, b$npar, b$offs, u)
  got <- trace_design_form(b$spec, G, d4, b$params, b$npar, tv, tu)
  if (length(b$units)) {
    dv <- refresh_direction(b$spec, b$design, M, b$params, b$npar, b$offs,
                            b$d3, tv, v, b$units)
    du <- refresh_direction(b$spec, b$design, M, b$params, b$npar, b$offs,
                            b$d3, tu, u, b$units)
    ac <- refresh_curv_amat(b$spec, b$design, M, b$params, b$npar, b$offs,
                            b$units, b$Hl)
    f2 <- lapply(b$units, function(un)
      refresh_dblock2(un, v[un$ra], u[un$ra], b$spec@n_obs))
    got <- got + trace_refresh4(b$spec, M, b$params, b$npar, b$Hl, dv, du,
                                b$units, G, b$d3, v, f2, ac)
  }
  trace <- sum(M * ref)
  abs(got - trace) / max(1e-12, abs(trace))
}

test_that("the twice-contracted fourth derivative matches a second difference", {
  # Three cells, and the middle one is what identifies each correction. On a
  # bilinear f the block's SECOND derivative is exactly zero, so whatever gap
  # remains there belongs to the term reading only the first -- the direction's
  # own predictor, which moves with beta because the block is a Jacobian.
  # Measured before these were written: 5.04e-08 on the fixed design, 2.32e-02
  # on the bilinear one and 7.60e-02 on the curved one.
  set.seed(12)
  n1 <- 200
  d1 <- data.frame(x = runif(n1, 0, 10))
  d1$y <- 1 + 0.3 * d1$x + stats::rnorm(n1, sd = 0.4)

  set.seed(21)
  nb <- 200
  db <- data.frame(x = runif(nb, 0.5, 3), z = runif(nb, -1, 1),
                   w = runif(nb, 0.5, 2),
                   grp = factor(rep(sprintf("g%d", 1:5), length.out = nb)))
  db$y <- 1.4 * db$x + 0.8 * db$z + 1.4 * 0.8 * db$w +
    stats::rnorm(nb, sd = 0.2)

  set.seed(9)
  nn <- 200
  dn <- data.frame(x = seq(0.2, 4, length.out = nn),
                   grp = factor(rep(sprintf("g%d", 1:5), length.out = nn)))
  dn$y <- 3 * exp(-0.6 * dn$x) + stats::rnorm(nn, sd = 0.15)

  cases <- list(
    list(y ~ s(x, k = 6), d1, 0L),
    list(y ~ nl(~ a * x + b * z + a * b * w, a ~ 0 + ridge(~ grp),
                start = list(b = 0.8)), db, 1L),
    list(y ~ 0 + nl(~ a * exp(-r * x), a ~ 0 + ridge(~ grp)), dn, 1L))
  for (cs in cases) {
    b <- u_trace_bits(cs[[1L]], cs[[2L]])
    expect_length(b$units, cs[[3L]])
    set.seed(2)
    v <- stats::rnorm(b$total) * 0.1
    u <- stats::rnorm(b$total) * 0.1
    expect_lt(u_trace_gap(b, v, u), 1e-5)
  }
})
