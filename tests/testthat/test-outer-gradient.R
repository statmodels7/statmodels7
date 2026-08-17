# The exact gradient of the marginal criterion.

set.seed(31)
n <- 250
dg <- data.frame(x = runif(n, -2, 2), z = runif(n))
dg$y <- sin(1.4 * dg$x) + stats::rnorm(n, sd = 0.3)

# The criterion as a function of the free hyperparameters, refitting the
# coefficients at each: this is what the outer search minimizes, and
# differentiating it numerically shares no arithmetic with the assembly the
# gradient uses.
crit_of_eta <- function(formula, data, method) {
  # the base point is the PROBE value, not the estimated optimum: this
  # differentiates the criterion at a point where its gradient is not zero,
  # which is the only place a relative comparison against numDeriv means
  # anything
  fit0 <- statmod(formula, distributions7::gaussian1_distrib(), data,
                  outer_criterion = NULL)
  spec <- fit0@spec
  design <- statmod_design(spec)
  blocks <- statmod_blocks(spec, design)
  idx <- outer_hyper_index(spec, blocks)
  basis <- integrated_basis(spec, design, method@kind)
  list(
    idx = idx,
    fn = function(eta) {
      hy <- eta_to_hyper(eta, idx, fit0@hyper)
      cf <- fit_at_hyper(formula, distributions7::gaussian1_distrib(),
                         data, hy)$coefficients
      statmod_marginal(spec, design, cf, hy, method, basis = basis)$value
    },
    gr = function(eta) {
      hy <- eta_to_hyper(eta, idx, fit0@hyper)
      cf <- fit_at_hyper(formula, distributions7::gaussian1_distrib(),
                         data, hy)$coefficients
      statmod_marginal_grad(spec, design, cf, hy, method, idx, basis)
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

test_that("the gradient is exact where the block moves with the coefficients", {
  skip_if_not_installed("numDeriv")
  # nl()'s block is the Jacobian, so it moves, and two pieces the layer used to
  # miss are supplied by the TERM rather than differenced here:
  # term_block_contract() gives dX/dbeta, which enters dK/dbeta (u_refresh) and
  # the mode's own curvature (mode_curvature). Measured against a finite
  # difference of the criterion with the mode refitted, this went 4.6e-03 ->
  # 1.4e-03 (the stale block) -> 1.5e-04 (dK/dbeta) -> 4.0e-09 (the mode).
  set.seed(9)
  np <- 40; m <- 8
  dn <- data.frame(x = rep(seq(0.2, 4, length.out = np), m),
                   grp = factor(rep(sprintf("g%d", seq_len(m)), each = np)))
  a_g <- 3 + stats::rnorm(m, sd = 0.4)
  dn$y <- a_g[dn$grp] * exp(-0.6 * dn$x) + stats::rnorm(nrow(dn), sd = 0.15)
  f <- y ~ nl(~ a * exp(-r * x), a ~ 0 + ridge(~ grp))
  for (hh in c("observed", "expected")) {
    h <- crit_of_eta(f, dn, reml(hessian = hh))
    eta <- h$eta0 + 0.35
    # the reference refits the mode at the inner default, so it carries that
    # tolerance; the gradient itself measures 4.0e-09 against a mode located
    # to 1e-10
    expect_equal(h$gr(eta), numDeriv::grad(h$fn, eta), tolerance = 1e-4)
  }
})

test_that("a penalty inside a refreshable term reaches the same optimum", {
  # ⚠️ KNOWN LIMITATION, pinned by its CONSEQUENCE rather than by its size.
  # nl(), seg(), jump() and jseg() register term_refresh(), so their design
  # block moves with the coefficients; u_vector() assembles tr(M dK/dbeta) as
  # though X were constant, and with X = X(beta) there is a second
  # contribution through dX/dbeta that nothing computes. Measured on
  # nl(a ~ 0 + ridge(~grp)), the exact gradient disagrees with a finite
  # difference of its own criterion by 5.3e-03 at n = 320 and 3.9e-04 at
  # n = 960 -- FLAT in the inner tolerance across three decades while the
  # mode's score fell an order, and the SAME ridge on a fixed block agrees to
  # 3.6e-08, which is what localizes it.
  #
  # It is left in place because the consequence is small and refusing it would
  # cost more than it saves: the search reaches the same hyperparameter in a
  # third of the evaluations. This test guards THAT, so a gap that grew would
  # fail here even though the gradient's own error is not asserted.
  set.seed(9)
  np <- 40; m <- 8
  dn <- data.frame(x = rep(seq(0.2, 4, length.out = np), m),
                   grp = factor(rep(sprintf("g%d", seq_len(m)), each = np)))
  a_g <- 3 + stats::rnorm(m, sd = 0.4)
  dn$y <- a_g[dn$grp] * exp(-0.6 * dn$x) + stats::rnorm(nrow(dn), sd = 0.15)
  f <- y ~ nl(~ a * exp(-r * x), a ~ 0 + ridge(~ grp))
  fe <- statmod(f, distributions7::gaussian1_distrib(), dn,
                outer_criterion = reml("observed"))
  fd <- statmod(f, distributions7::gaussian1_distrib(), dn,
                outer_criterion = reml("observed"),
                outer_optimizer = optimizers7::nelder_mead())
  expect_equal(unlist(lapply(fe@hyper, unlist)),
               unlist(lapply(fd@hyper, unlist)), tolerance = 5e-3)
  expect_equal(as.numeric(logLik(fe)), as.numeric(logLik(fd)),
               tolerance = 1e-6)
  # and the exact route is what it is for: fewer evaluations
  expect_lt(nrow(fe@history$outer), nrow(fd@history$outer))
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
  # The expected information needs the derivative in beta of -E[l''], which is
  # NOT -E[l''']: differentiating an expectation moves the measure as well as
  # the integrand. distributions7 carries it as distrib_dexpected_hessian()
  # and the route is available wherever the family answers -- a gaussian does.
  expect_true(outer_gradient_ok(spec, design, idx, reml("expected")))
  # but only at order 1: the criterion's own second derivative would want the
  # next order of the same object
  expect_false(outer_gradient_ok(spec, design, idx, reml("expected"), 2L))

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

# A panel whose level is driven by the model's OWN score. Simulating it from
# an exogenous series instead leaves the fit at a degenerate point -- no
# dynamics, the persistence at the unit root -- where the criterion is flat
# and nothing is worth differentiating.
sim_gas_panel <- function(seed, m_grp, ti, a = 0.25, b = 0.55, beta = 0.5) {
  set.seed(seed)
  id <- factor(rep(seq_len(m_grp), each = ti))
  n <- length(id)
  x <- stats::rnorm(n)
  om <- stats::rnorm(m_grp, 0.2, 0.35)
  y <- numeric(n)
  for (g in seq_len(m_grp)) {
    rows <- which(as.integer(id) == g)
    f <- om[g] / (1 - b)
    s <- 0
    for (k in seq_along(rows)) {
      f <- om[g] + a * s + b * f
      eta <- beta * x[rows[k]] + f
      y[rows[k]] <- eta + stats::rnorm(1)
      s <- y[rows[k]] - eta
    }
  }
  data.frame(id = id, t = rep(seq_len(ti), m_grp), x = x, y = y)
}

# the criterion and its gradient with the mode refitted, read at the REFITTED
# design: a structural term's own parameters live in the design's state, so a
# design built once and reused would read them at a non-stationary point
struct_harness <- function(form, data, method) {
  fit0 <- statmod(form, distributions7::gaussian1_distrib(), data,
                  outer_criterion = NULL)
  idx <- outer_hyper_index(fit0@spec,
                           statmod_blocks(fit0@spec,
                                          statmod_design(fit0@spec)))
  refit <- function(eta) {
    hy <- eta_to_hyper(eta, idx, fit0@hyper)
    a <- fit_at_hyper(form, distributions7::gaussian1_distrib(), data, hy)
    list(f = a, d = a$design, hy = hy,
         basis = integrated_basis(a$spec, a$design, method@kind))
  }
  list(idx = idx, fit0 = fit0, eta0 = hyper_to_eta(fit0@hyper, idx),
       fn = function(eta) {
         r <- refit(eta)
         statmod_marginal(r$f$spec, r$d, r$f$coefficients, r$hy, method,
                          basis = r$basis)$value
       },
       gr = function(eta) {
         r <- refit(eta)
         statmod_marginal_grad(r$f$spec, r$d, r$f$coefficients, r$hy, method,
                               idx, r$basis)
       })
}

# A CENTRAL DIFFERENCE and not numDeriv::grad, which is Richardson and pays
# eight refits per hyperparameter where two will do. It is a legitimate
# reference here because the criterion carries no difference of its own:
# measured, it agrees with the exact gradient to 1.2e-07 at h = 1e-3 and
# converges O(h^2), where a structural error would be O(1) relative.
central <- function(fn, eta, h = 1e-3) {
  vapply(seq_along(eta), function(j) {
    ep <- eta; em <- eta
    ep[j] <- ep[j] + h; em[j] <- em[j] - h
    a <- fn(ep); b <- fn(em)
    # the criterion does not exist everywhere: a hyperparameter far enough
    # out takes the inner fit somewhere it cannot reach a mode, and
    # statmod_marginal() returns NULL there. A reference that cannot be
    # formed is reported as such rather than raising a length error three
    # frames down, which is what the first version of this did.
    if (!length(a) || !length(b)) return(NA_real_)
    (a - b) / (2 * h)
  }, numeric(1))
}

test_that("the gradient reaches a penalty over a filter's own parameters", {
  # Where the penalty covers a structural term's parameters the determinant
  # spans them too, so the chain term reads the recursion's THIRD derivative
  # -- contracted in the one direction the mode moves in, which is what
  # modelterms7::term_third() propagates.
  skip_on_cran()
  dp <- sim_gas_panel(7, 3L, 35L)
  form <- y ~ x +
    gas(p = 1, q = 1, omega ~ random(~1 | id), by = id, time = t)
  h <- struct_harness(form, dp, reml(hessian = "observed"))
  # the fit must be a sane one, or the reference is differentiating a flat
  # direction rather than the criterion
  expect_true(h$fit0@converged)
  for (d in c(0, -0.5)) {
    eta <- h$eta0 + d
    expect_equal(h$gr(eta), central(h$fn, eta), tolerance = 1e-4,
                 info = sprintf("eta0 %+0.1f", d))
  }
})

test_that("the structural gradient is right under ml and beside a smooth", {
  skip_on_cran()
  dp <- sim_gas_panel(7, 3L, 35L)
  form <- y ~ x +
    gas(p = 1, q = 1, omega ~ random(~1 | id), by = id, time = t)
  hm <- struct_harness(form, dp, ml(hessian = "observed"))
  expect_equal(hm$gr(hm$eta0 + 0.2), central(hm$fn, hm$eta0 + 0.2),
               tolerance = 1e-3)

  # two hyperparameters of DIFFERENT kinds: once a structural penalty is
  # present the determinant is the joint one, so the ordinary smooth's
  # gradient goes through the same assembly and must still be right
  dp$z <- stats::runif(nrow(dp), -2, 2)
  dp$y <- dp$y + sin(1.6 * dp$z)
  f2 <- y ~ x + s(z, k = 6) +
    gas(p = 1, q = 1, omega ~ random(~1 | id), by = id, time = t)
  h2 <- struct_harness(f2, dp, reml(hessian = "observed"))
  expect_identical(nrow(h2$idx), 2L)
  eta <- h2$eta0 + c(0.15, -0.15)
  ref <- central(h2$fn, eta)
  # on a small panel the criterion is not available at every probe, and a
  # reference that could not be formed is not evidence either way
  skip_if(anyNA(ref), "the criterion is unavailable at a probe point")
  expect_equal(h2$gr(eta), ref, tolerance = 1e-3)
})

test_that("an outer step the search takes back moves nothing", {
  # The coefficients were always protected -- `state$beta` is written only
  # once a point is known to be usable -- and a structural term's own
  # parameters were NOT: they live in the design's structural state, an
  # environment the inner fit writes into as it goes, and the structural
  # sub-fit stores its optimizer's last point whether it converged or not.
  # An unavailable point therefore moved the filter and the next evaluation
  # started from wherever the failure had left it, which ratchets.
  #
  # The exact gradient is what made it reachable: lbfgs takes its first step
  # as the full gradient, which grows with the number of penalized
  # coordinates, so it lands far out where a simplex never goes. Measured
  # before the fix, a panel of twenty groups reported thirty consecutive
  # unavailable points and no fit; after it, ten evaluations.
  # Two things had to be right before this holds, and both were found here.
  # The restore stopped the ratchet: without it the search lost every
  # backtracked point -- thirty consecutive unavailable ones and no fit at
  # twenty groups. It was not sufficient, because bfgs and lbfgs took their
  # FIRST direction as the raw gradient, which on this criterion grows with
  # the number of penalized coordinates, so the first trial point was an
  # arbitrary distance out; optimizers7 0.3.0 scales it to a step of order
  # one in the parameters. Measured after both, on panels of 3, 8, 12 and 20
  # groups: 8, 6, 7 and 5 evaluations against 33, 55, 23 and 41, the
  # criterion agreeing to the printed digit every time.
  skip_on_cran()
  dp <- sim_gas_panel(7, 4L, 45L)
  form <- y ~ x +
    gas(p = 1, q = 1, omega ~ random(~1 | id), by = id, time = t)
  fast <- statmod(form, distributions7::gaussian1_distrib(), dp,
                  outer_criterion = reml(hessian = "observed"))
  slow <- statmod(form, distributions7::gaussian1_distrib(), dp,
                  outer_criterion = reml(hessian = "observed"),
                  outer_optimizer = optimizers7::nelder_mead())
  expect_true(fast@converged)
  expect_equal(fast@criterion, slow@criterion, tolerance = 1e-5)
  expect_equal(fast@hyper$mu[[1L]][["sigma"]],
               slow@hyper$mu[[1L]][["sigma"]], tolerance = 1e-3)
  expect_lt(nrow(fast@history$outer), nrow(slow@history$outer))
})

test_that("the exact route is asked of the term, not of its class", {
  dp <- sim_gas_panel(9, 3L, 30L)
  sp <- statmod_spec(y ~ x +
                       gas(p = 1, q = 1, omega ~ random(~1 | id),
                           by = id, time = t),
                     distributions7::gaussian1_distrib(), dp)
  de <- statmod_design(sp)
  ix <- outer_hyper_index(sp, statmod_blocks(sp, de))
  expect_true(structural_penalized(sp, de))
  expect_true(outer_gradient_ok(sp, de, ix, reml("observed")))
  # the criterion's own SECOND derivative would ask for a fourth order
  # through the recursion, which is not written: newton is not offered there
  expect_false(outer_gradient_ok(sp, de, ix, reml("observed"), order = 2L))

  # and the question is put to the term. A gas term answers term_third; an
  # additive one inherits the zero that is right for it; a regime term bends
  # the predictor and has not written one, so it must NOT be mistaken for
  # either
  tm <- sp@terms[["mu"]][[grep("^gas", names(sp@terms[["mu"]]))]]
  expect_true(answers_term_third(tm))
  expect_true(answers_term_third(modelterms7::term_build(
    modelterms7::linpar(~x), dp)))
  expect_false(answers_term_third(modelterms7::term_build(
    modelterms7::regime(k = 2), dp)))
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
