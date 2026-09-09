# The outer Hessian of a model carrying a structural term, which is one
# central difference of the exact gradient rather than an assembly. What the
# checks below rest on is that the SAME routine, run on a model where the
# analytic route exists, reproduces it: that comparison is against an exact
# quantity and needs no reference of its own.

stencil_smooth <- function(n = 300, seed = 4L) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(n, -2, 2))
  d$y <- sin(2 * d$x) + stats::rnorm(n, 0, 0.4)
  d
}

# a score-driven filter over a panel, its level developed over a ridge, which
# is the shape whose hyperparameter had no standard error at all
stencil_panel <- function(m = 20L, ni = 30L, seed = 5L) {
  set.seed(seed)
  b <- stats::rnorm(m, 0, 0.4)
  d <- data.frame(id = factor(rep(seq_len(m), each = ni)))
  d$y <- as.numeric(unlist(lapply(seq_len(m), function(i) {
    e <- numeric(ni)
    f <- 0.5 + b[i]
    for (t in seq_len(ni)) {
      e[t] <- stats::rnorm(1, f, 1)
      f <- 0.5 + b[i] + 0.3 * (e[t] - f) + 0.5 * (f - 0.5 - b[i])
    }
    e
  })))
  d
}

stencil_parts <- function(fit) {
  spec <- fit@spec
  design <- statmod_design(spec)
  list(spec = spec, design = design,
       idx = outer_hyper_index(spec, statmod_blocks(spec, design)),
       method = fit@methods$outer,
       basis = integrated_basis(spec, design, fit@methods$outer@kind),
       coef = fit@coefficients, hyper = fit@hyper)
}

test_that("the stencil reproduces the analytic Hessian where one exists", {
  # THE CHECK THAT RESTS ON NOTHING ELSE: on a model with no structural term
  # the assembly is exact, so the difference of the exact gradient has an
  # exact quantity to be compared against rather than a reference of its own.
  fit <- statmod(y ~ s(x, k = 10), gaussian1_distrib(), stencil_smooth(),
                 outer_criterion = reml())
  p <- stencil_parts(fit)
  A <- as.matrix(statmod_marginal_hess(p$spec, p$design, p$coef, p$hyper,
                                       p$method, p$idx, p$basis))
  B <- statmod_hess_stencil(p$spec, p$design, p$coef, p$hyper, p$method,
                            p$idx, p$basis)
  expect_false(is.null(B))
  expect_equal(as.numeric(B), as.numeric(A), tolerance = 1e-3)
  # and the two-step guard passes there, which is what says the reading is
  # resolved rather than merely close on one step
  expect_false(is.null(statmod_hess_stencil(p$spec, p$design, p$coef, p$hyper,
                                            p$method, p$idx, p$basis,
                                            h = 3 * hess_stencil_step())))
})

test_that("a model with no structural term is untouched", {
  # the route is taken on the design's own answer, so a model carrying no
  # such term must reach the assembly and reach it identically
  fit <- statmod(y ~ s(x, k = 10), gaussian1_distrib(), stencil_smooth(),
                 outer_criterion = reml())
  p <- stencil_parts(fit)
  expect_identical(length(attr(p$design, "structural")), 0L)
  expect_true(outer_gradient_ok(p$spec, p$design, p$idx, p$method, 2L))
  a <- statmod_marginal_hess(p$spec, p$design, p$coef, p$hyper, p$method,
                             p$idx, p$basis)
  b <- statmod_marginal_hess(p$spec, p$design, p$coef, p$hyper, p$method,
                             p$idx, p$basis, inner = iwls())
  expect_identical(a, b)
  expect_s3_class(fit@methods$search, class(optimizers7::newton())[[1L]])
})

test_that("a structural term takes the stencil, and the number is real", {
  fit <- statmod(y ~ 0 + gas(p = 1, q = 1, by = id, omega ~ 0 + ridge(~ id)),
                 gaussian1_distrib(), stencil_panel(),
                 outer_criterion = reml())
  p <- stencil_parts(fit)
  expect_gt(length(attr(p$design, "structural")), 0L)
  # order 2 is refused for the SEARCH, so the default stays the gradient-only
  # optimizer and no stencil is paid per iteration
  expect_true(outer_gradient_ok(p$spec, p$design, p$idx, p$method, 1L))
  expect_false(outer_gradient_ok(p$spec, p$design, p$idx, p$method, 2L))
  expect_false(inherits(fit@methods$search, class(optimizers7::newton())[[1L]]))

  H <- statmod_marginal_hess(p$spec, p$design, p$coef, p$hyper, p$method,
                             p$idx, p$basis)
  expect_false(is.null(H))
  expect_identical(dim(as.matrix(H)), c(1L, 1L))
  # THE NEGATIVE CONTROL, and it is what the assembled route would fail: over
  # the stacked coefficients the curvature of this criterion comes out at the
  # rounding, 1.09e-06 where it is of order one, so a reading near zero is
  # the defect rather than a small number
  expect_gt(abs(as.numeric(H)), 0.5)
  expect_lt(as.numeric(H), 0)
})

test_that("the stencil agrees with a second difference of the criterion", {
  skip_on_cran()
  fit <- statmod(y ~ 0 + gas(p = 1, q = 1, by = id, omega ~ 0 + ridge(~ id)),
                 gaussian1_distrib(), stencil_panel(),
                 outer_criterion = reml())
  p <- stencil_parts(fit)
  H <- as.numeric(statmod_marginal_hess(p$spec, p$design, p$coef, p$hyper,
                                        p$method, p$idx, p$basis))
  # the reference shares no arithmetic with the gradient: it reads the
  # criterion alone, twice differenced, with the mode refitted at each point
  cfg <- inner_settings(iwls())
  blocks <- statmod_blocks(p$spec, p$design)
  eta0 <- hyper_to_eta(p$hyper, p$idx)
  beta0 <- unlist(p$coef[p$spec@distrib@params], use.names = FALSE)
  sst <- statmod_structural_state(p$design)
  z0 <- if (is.null(sst)) NULL else sst$zeta
  cval <- function(eta) {
    if (!is.null(z0)) sst$zeta <- z0
    hy <- eta_to_hyper(eta, p$idx, p$hyper)
    r <- statmod_alternate(p$spec, p$design, blocks, hy, iwls(), beta0,
                           cfg$expected, cfg$approx, cfg$maxit, cfg$tol,
                           verbosity(0), hold_refresh = TRUE)
    statmod_marginal(p$spec, p$design, r$obj$split(r$par), hy, p$method,
                     cfg$approx, p$basis, NULL)$value
  }
  h <- 1e-2
  ref <- (cval(eta0 + h) - 2 * cval(eta0) + cval(eta0 - h)) / h^2
  expect_equal(H, ref, tolerance = 1e-3)
})

test_that("the stencil refuses where the curvature is not resolved", {
  skip_on_cran()
  # a mixed covariance class whose correlation the search leaves at the
  # boundary of the spherical chart. The gradient itself is 1e-3 there, so a
  # difference of it divided by 2h is noise, and the two steps disagree by
  # 6e-01 against 1e-06 where the mode is well located.
  set.seed(77L)
  m <- 8L; ni <- 30L
  z <- matrix(stats::rnorm(2 * m), m, 2) %*% chol(diag(c(0.36, 0.25)))
  g <- factor(rep(seq_len(m), each = ni))
  y <- unlist(lapply(seq_len(m), function(i) {
    f <- 0
    out <- numeric(ni)
    for (t in seq_len(ni)) {
      mu <- z[i, 1L] + f
      out[t] <- stats::rnorm(1, mu, 1)
      f <- 0.25 * exp(z[i, 2L]) * (out[t] - mu) + 0.6 * f
    }
    out
  }))
  d <- data.frame(g = g, y = y)
  fit <- statmod(y ~ random(~ 1 | u | g) +
                   gas(p = 1, q = 1, by = g, alpha1 ~ 1 + random(~ 1 | u | g)),
                 gaussian1_distrib(), d, outer_criterion = reml())
  p <- stencil_parts(fit)
  # the premise of the case, asserted so that a fit landing elsewhere skips
  # rather than passing for the wrong reason
  eta <- hyper_to_eta(p$hyper, p$idx)
  skip_if(max(abs(eta)) < statmod_certificate(fit)$edge,
          "this fit did not reach the chart's boundary")
  expect_null(statmod_marginal_hess(p$spec, p$design, p$coef, p$hyper,
                                    p$method, p$idx, p$basis))
  expect_null(statmod_hyper_vcov(p$spec, p$design, p$coef, p$hyper, p$method))
})

test_that("a filter's hyperparameter carries a standard error", {
  skip_on_cran()
  fit <- statmod(y ~ 0 + gas(p = 1, q = 1, by = id, omega ~ 0 + ridge(~ id)),
                 gaussian1_distrib(), stencil_panel(),
                 outer_criterion = reml())
  spec <- fit@spec
  V <- statmod_hyper_vcov(spec, statmod_design(spec), fit@coefficients,
                          fit@hyper, fit@methods$outer,
                          inner = fit@methods$smooth)
  expect_false(is.null(V))
  expect_identical(dim(V), c(1L, 1L))
  expect_gt(V[[1L, 1L]], 0)
  # and the summary prints it rather than an empty cell. A structural term's
  # own penalty is reported in the block of the parameter it develops, so the
  # row is under that component and not in the term's own table.
  cmp <- summary(fit)@tables$mu[[1L]]$components
  expect_identical(names(cmp), "omega")
  r <- cmp$omega$table
  r <- r[r$role == "estimated", ]
  expect_identical(nrow(r), 1L)
  expect_identical(r$name, "precision")
  expect_true(is.finite(r$se))
  expect_lt(r$lower, r$estimate)
  expect_gt(r$upper, r$estimate)
})
