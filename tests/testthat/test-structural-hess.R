# The exact outer Hessian of a model carrying a structural term:
# statmod_structural_hess(), the joint twin of statmod_marginal_hess().
#
# What is asserted, in the order the plan asks for it:
#   - the matrix against a central difference of the EXACT gradient with the
#     mode refitted, at three steps, reading the CONVERGENCE RATE and not one
#     tolerance -- clean O(h^2) is the signature that separates a correct
#     derivative from one missing a term, which is flat in h;
#   - the same matrix against statmod_hess_stencil(), which shares no
#     arithmetic with it, being four refits of the gradient;
#   - an IDENTITY between the assembled derivative and the traced one, which
#     no tolerance chooses;
#   - and a model with NO structural term untouched.

sh_panel <- function(seed = 51L, m = 6L, ni = 20L) {
  set.seed(seed)
  g <- factor(rep(seq_len(m), each = ni))
  n <- m * ni
  dev <- stats::rnorm(m, 0, 0.6)
  y <- numeric(n)
  for (l in seq_len(m)) {
    rows <- which(as.integer(g) == l)
    f <- 0.2
    a <- exp(log(0.3) + dev[l])
    for (t in seq_along(rows)) {
      y[rows[t]] <- stats::rnorm(1, f, 1)
      f <- 0.2 + a * (y[rows[t]] - f) + 0.5 * f
    }
  }
  data.frame(y = y, g = g)
}

# the same panel with the LEVEL developed as well, so the criterion carries
# two hyperparameters and the Hessian an off-diagonal. A one-hyperparameter
# check cannot see a term placed on the wrong direction, the two directions
# coinciding there -- the trap the symmetry of an index pair at p = q = 1
# already records one layer down.
sh_panel2 <- function(seed = 51L, m = 8L, ni = 25L) {
  set.seed(seed)
  g <- factor(rep(seq_len(m), each = ni))
  n <- m * ni
  dv <- stats::rnorm(m, 0, 0.6)
  do <- stats::rnorm(m, 0, 0.4)
  y <- numeric(n)
  for (l in seq_len(m)) {
    rows <- which(as.integer(g) == l)
    f <- 0.2
    a <- exp(log(0.3) + dv[l])
    om <- 0.2 + do[l]
    for (t in seq_along(rows)) {
      y[rows[t]] <- stats::rnorm(1, f, 1)
      f <- om + a * (y[rows[t]] - f) + 0.5 * f
    }
  }
  data.frame(y = y, g = g)
}

sh_parts <- function(fit) {
  spec <- fit@spec
  design <- statmod_design(spec)
  blocks <- statmod_blocks(spec, design)
  method <- fit@methods$outer
  list(spec = spec, design = design, blocks = blocks,
       idx = outer_hyper_index(spec, blocks), method = method,
       basis = integrated_basis(spec, design, method@kind),
       coef = fit@coefficients, hyper = fit@hyper)
}

# the exact gradient at a hyperparameter, with the mode refitted from the SAME
# start and the structural state restored before every probe -- the ratchet
# guard statmod_hess_stencil() records, and the reason two probes differ in
# the hyperparameter alone
sh_grad <- function(fit, p) {
  inner <- iwls()
  cfg <- inner_settings(inner)
  beta0 <- unlist(p$coef[p$spec@distrib@params], use.names = FALSE)
  sst <- statmod_structural_state(p$design)
  z0 <- if (is.null(sst)) NULL else sst$zeta
  function(eta) {
    if (!is.null(z0)) sst$zeta <- z0
    hy <- eta_to_hyper(eta, p$idx, p$hyper)
    r <- statmod_alternate(p$spec, p$design, p$blocks, hy, inner, beta0,
                           cfg$expected, cfg$approx, cfg$maxit, cfg$tol,
                           verbosity(0), hold_refresh = TRUE)
    statmod_marginal_grad(p$spec, p$design, r$obj$split(r$par), hy, p$method,
                          p$idx, p$basis)
  }
}

test_that("the structural Hessian converges O(h^2) onto the exact gradient", {
  skip_on_cran()
  d <- sh_panel()
  fit <- statmod(y ~ gas(p = 1, q = 1, by = g, alpha1 ~ 1 + random(~ 1 | g)),
                 gaussian1_distrib(), d, outer_criterion = reml())
  p <- sh_parts(fit)
  expect_gt(length(attr(p$design, "structural")), 0L)
  A <- statmod_structural_hess(p$spec, p$design, p$coef, p$hyper, p$method,
                               p$idx, p$basis)
  expect_false(is.null(A))
  expect_identical(A, t(A))

  gf <- sh_grad(fit, p)
  eta0 <- hyper_to_eta(p$hyper, p$idx)
  nh <- length(eta0)
  fd <- function(h) {
    Hn <- matrix(0, nh, nh)
    for (m in seq_len(nh)) {
      ep <- eta0; ep[[m]] <- ep[[m]] + h
      em <- eta0; em[[m]] <- em[[m]] - h
      Hn[, m] <- (gf(ep) - gf(em)) / (2 * h)
    }
    (Hn + t(Hn)) / 2
  }
  sc <- max(abs(A))
  gaps <- vapply(c(1e-2, 3e-3, 1e-3),
                 function(h) max(abs(A - fd(h))) / sc, numeric(1))
  # the RATE, not the size: a missing term is flat in h and a badly located
  # mode grows as 1/h, so what says the derivative is right is that the gap
  # falls like h^2. The two ratios are 11.1 and 9.0 to within the reference's
  # own noise.
  expect_gt(gaps[[1L]] / gaps[[2L]], 5)
  expect_gt(gaps[[2L]] / gaps[[3L]], 4)
  expect_lt(gaps[[3L]], 1e-4)

  # and against the stencil, which is four refits of the same gradient and
  # shares no arithmetic with the assembly
  S <- statmod_hess_stencil(p$spec, p$design, p$coef, p$hyper, p$method,
                            p$idx, p$basis, inner = iwls())
  expect_false(is.null(S))
  expect_equal(as.numeric(A), as.numeric(S), tolerance = 1e-4)
})

test_that("the off-diagonal of the structural Hessian is exact too", {
  skip_on_cran()
  # TWO hyperparameters, so the pair loop's own arithmetic is exercised: at
  # one hyperparameter the two directions coincide and a term placed on the
  # wrong one is invisible, which is why this case exists beside the last.
  d <- sh_panel2()
  fit <- statmod(y ~ gas(p = 1, q = 1, by = g, omega ~ 1 + random(~ 1 | g),
                         alpha1 ~ 1 + random(~ 1 | g)),
                 gaussian1_distrib(), d, outer_criterion = reml())
  p <- sh_parts(fit)
  expect_identical(nrow(p$idx), 2L)
  gf <- sh_grad(fit, p)
  eta0 <- hyper_to_eta(p$hyper, p$idx)
  # THE PREMISE, asserted so that a fit landing elsewhere fails here and not
  # on the arithmetic below. A criterion's second derivative is a statement
  # about the point it is read at, and away from a mode nothing verifies it:
  # measured on a panel of six groups the search stops with the outer
  # gradient at -734, the assembled curvature reads 1.8e+06, the stencil
  # refuses and a difference of the gradient is 80 per cent asymmetric at
  # every step.
  expect_lt(max(abs(gf(eta0))), 1e-3)
  A <- statmod_structural_hess(p$spec, p$design, p$coef, p$hyper, p$method,
                               p$idx, p$basis)
  expect_false(is.null(A))
  # the off-diagonal is not zero, so the assertion below is about something
  expect_gt(abs(A[1L, 2L]) / max(abs(A)), 1e-3)

  fd <- function(h) {
    Hn <- matrix(0, 2L, 2L)
    for (m in 1:2) {
      ep <- eta0; ep[[m]] <- ep[[m]] + h
      em <- eta0; em[[m]] <- em[[m]] - h
      Hn[, m] <- (gf(ep) - gf(em)) / (2 * h)
    }
    (Hn + t(Hn)) / 2
  }
  sc <- max(abs(A))
  g1 <- max(abs(A - fd(1e-2))) / sc
  g2 <- max(abs(A - fd(3e-3))) / sc
  expect_gt(g1 / g2, 5)
  expect_lt(g2, 1e-3)
  # and the stencil, which refuses on the small panel above, agrees here
  S <- statmod_hess_stencil(p$spec, p$design, p$coef, p$hyper, p$method,
                            p$idx, p$basis, inner = iwls())
  expect_false(is.null(S))
  expect_equal(as.numeric(A), as.numeric(S), tolerance = 1e-4)
})

test_that("the assembled joint derivative and its trace are one quantity", {
  skip_on_cran()
  # structural_dk_matrix() assembles dK/du[v] and structural_chain_extra()
  # traces it without assembling. The two share the recursion's own pieces
  # and nothing else, so the identity between them is what says the assembly
  # places every term where the trace expects it -- and it is an identity
  # rather than a tolerance because both are exact.
  d <- sh_panel(seed = 12L, m = 5L, ni = 18L)
  fit <- statmod(y ~ gas(p = 1, q = 1, by = g, alpha1 ~ 1 + random(~ 1 | g)),
                 gaussian1_distrib(), d, outer_criterion = reml())
  p <- sh_parts(fit)
  jd <- joint_design_rows(p$spec, p$design, p$coef)
  K <- statmod_marginal_full(p$spec, p$design, p$coef, p$hyper, NULL)
  M <- chol2inv(chol(K))
  st <- structural_grad_parts(p$spec, p$design, p$coef, jd, M)
  set.seed(3)
  for (i in 1:3) {
    v <- stats::rnorm(length(jd$keep), 0, 0.3)
    Tv <- structural_dk_matrix(p$spec, p$design, jd, st, v)
    tr <- sum(st$u * v) +
      structural_chain_extra(p$spec, p$design, jd, M, st, v)
    expect_equal(sum(M * Tv), tr, tolerance = 1e-10)
    expect_identical(Tv, t(Tv))
  }
})

test_that("the second-order trace is symmetric in its two directions", {
  skip_on_cran()
  # tr(M d2K/du2[v, w]) is a second derivative, so it cannot depend on the
  # order the two directions are given in. Nine terms contribute and the
  # swap maps each onto another, so a term placed on the wrong direction
  # breaks it.
  d <- sh_panel(seed = 12L, m = 5L, ni = 18L)
  fit <- statmod(y ~ gas(p = 1, q = 1, by = g, alpha1 ~ 1 + random(~ 1 | g)),
                 gaussian1_distrib(), d, outer_criterion = reml())
  p <- sh_parts(fit)
  jd <- joint_design_rows(p$spec, p$design, p$coef)
  K <- statmod_marginal_full(p$spec, p$design, p$coef, p$hyper, NULL)
  M <- chol2inv(chol(K))
  st <- structural_grad_parts(p$spec, p$design, p$coef, jd, M)
  D5 <- distributions7::distrib_deriv5(p$spec@distrib, p$spec@response,
                                       jd$ev$theta, scale = "link")
  blk4 <- .structural_blocks(jd$params, jd$ap, jd$V, st$H, st$D3, st$D4,
                             jd$n, D5)
  set.seed(5)
  v <- stats::rnorm(length(jd$keep), 0, 0.3)
  w <- stats::rnorm(length(jd$keep), 0, 0.3)
  a <- structural_chain_extra2(p$spec, p$design, jd, M, st, blk4, v, w)
  b <- structural_chain_extra2(p$spec, p$design, jd, M, st, blk4, w, v)
  expect_equal(a, b, tolerance = 1e-9)
  expect_true(is.finite(a))
})

test_that("a model with no structural term is untouched by the joint route", {
  # the route is taken on the design's own answer, so a model carrying no
  # such term must not reach it at all
  set.seed(8)
  n <- 200L
  d <- data.frame(x = stats::runif(n, -2, 2))
  d$y <- sin(d$x) + stats::rnorm(n, 0, 0.5)
  fit <- statmod(y ~ s(x, k = 8), gaussian1_distrib(), d,
                 outer_criterion = reml())
  p <- sh_parts(fit)
  expect_identical(length(attr(p$design, "structural")), 0L)
  expect_null(structural_term_of(p$spec, p$design))
  # statmod_marginal_hess() reaches its own assembly, and calling the joint
  # one directly returns NULL rather than a number
  expect_null(statmod_structural_hess(p$spec, p$design, p$coef, p$hyper,
                                      p$method, p$idx, p$basis))
  expect_false(is.null(statmod_marginal_hess(p$spec, p$design, p$coef,
                                             p$hyper, p$method, p$idx,
                                             p$basis)))
})

test_that("a filter's hyperparameter keeps its standard error", {
  skip_on_cran()
  # the consumer this route exists for: the interval summary() prints for a
  # hyperparameter estimated by a marginal criterion, which used to come
  # from the stencil and now comes from the assembly
  d <- sh_panel()
  fit <- statmod(y ~ gas(p = 1, q = 1, by = g, alpha1 ~ 1 + random(~ 1 | g)),
                 gaussian1_distrib(), d, outer_criterion = reml())
  p <- sh_parts(fit)
  V <- statmod_hyper_vcov(p$spec, p$design, p$coef, p$hyper, p$method,
                          inner = fit@methods$smooth)
  expect_false(is.null(V))
  expect_gt(V[[1L, 1L]], 0)
  txt <- paste(capture.output(print(summary(fit))), collapse = "\n")
  expect_match(txt, "effect sd")
})
