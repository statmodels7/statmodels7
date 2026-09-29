# The coefficients a marginal criterion estimates: reml(marginal =) and
# ml(marginal =).

marginal_data <- function(n = 60, seed = 1) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n), z = stats::runif(n))
  d$y <- 1 + 2 * d$x - d$z + stats::rnorm(n, sd = 1.5)
  d
}

test_that("reml() returns the sigma and the criterion of lm(REML)", {
  d <- marginal_data()
  g <- distributions7::gaussian1_distrib()
  l <- stats::lm(y ~ x + z, d)
  f <- statmod(y ~ x + z, g, d)
  expect_equal(exp(f@coefficients$sigma), summary(l)$sigma, tolerance = 1e-7)
  expect_equal(f@criterion, as.numeric(stats::logLik(l, REML = TRUE)),
               tolerance = 1e-8)
  expect_equal(unname(f@coefficients$mu), unname(stats::coef(l)),
               tolerance = 1e-8)
})

test_that("ml() returns sqrt(rss/n), and marginal = 'none' the joint mode", {
  d <- marginal_data()
  g <- distributions7::gaussian1_distrib()
  l <- stats::lm(y ~ x + z, d)
  ml_sigma <- sqrt(sum(stats::residuals(l)^2) / nrow(d))
  f_ml <- statmod(y ~ x + z, g, d, outer_criterion = ml())
  f_none <- statmod(y ~ x + z, g, d, outer_criterion = reml(marginal = "none"))
  expect_equal(exp(f_ml@coefficients$sigma), ml_sigma, tolerance = 1e-7)
  expect_equal(exp(f_none@coefficients$sigma), ml_sigma, tolerance = 1e-7)
  # "none" with no penalty runs no criterion at all
  expect_true(is.na(f_none@criterion))
})

test_that("with every parameter named, reml() and ml() are one criterion", {
  set.seed(11)
  m <- 20L
  gr <- factor(rep(seq_len(m), each = 8L))
  x <- stats::runif(length(gr))
  eta <- 0.3 + 0.8 * x + stats::rnorm(m, 0, 0.5)[gr]
  dat <- data.frame(g = gr, x = x,
                    y = stats::rgamma(length(gr), shape = 3,
                                      rate = 3 / exp(eta)))
  d <- distributions7::gamma1_distrib()
  fr <- statmod(y ~ x + random(~ 1 | g), d, dat,
                outer_criterion = reml(marginal = "all"))
  fm <- statmod(y ~ x + random(~ 1 | g), d, dat,
                outer_criterion = ml(marginal = "all"))
  expect_equal(fr@criterion, fm@criterion, tolerance = 1e-10)
  expect_equal(unlist(fr@coefficients), unlist(fm@coefficients),
               tolerance = 1e-8)
})

test_that("the variance carries the criterion's curvature in sigma", {
  d <- marginal_data()
  g <- distributions7::gaussian1_distrib()
  l <- stats::lm(y ~ x + z, d)
  f <- statmod(y ~ x + z, g, d)
  V <- vcov(f, readable = FALSE)
  expect_equal(unname(V[1:3, 1:3]), unname(stats::vcov(l)), tolerance = 1e-8)
  # -d2 l_M / d(log sigma)^2 is 2(n - p) at the REML optimum
  expect_equal(V[4, 4], 1 / (2 * (nrow(d) - 3)), tolerance = 1e-8)
  expect_identical(statmod_certificate(f)$state, "converged")
})

test_that("the exact gradient and Hessian over (eta, gamma) match differences", {
  set.seed(2)
  n <- 200
  d <- data.frame(x = stats::runif(n), z = stats::runif(n),
                  g = factor(sample(15, n, TRUE)))
  d$y <- sin(3 * d$x) + stats::rnorm(15, 0, 0.4)[d$g] +
    stats::rnorm(n, sd = exp(-0.5 + 0.8 * d$z))
  G <- distributions7::gaussian1_distrib()
  method <- reml()
  fit <- statmod(y ~ x + random(~ 1 | g) | sigma ~ z, G, d,
                 outer_criterion = method)
  spec <- fit@spec
  design <- statmod_design(spec)
  blocks <- statmod_blocks(spec, design)
  idx <- outer_hyper_index(spec, blocks)
  gam <- marginal_coords(spec, design, method)
  nh <- nrow(idx)
  ng <- length(gam$where)
  expect_identical(ng, 2L)
  basis <- integrated_basis(spec, design, method@kind, gamma = TRUE)
  beta0 <- unlist(fit@coefficients[spec@distrib@params], use.names = FALSE)
  set.seed(5)
  v0 <- c(hyper_to_eta(fit@hyper, idx), beta0[gam$where]) +
    stats::rnorm(nh + ng, 0, 0.15)
  at <- function(v) {
    hy <- eta_to_hyper(v[seq_len(nh)], idx, fit@hyper)
    sp <- statmod_hold(spec, gam, v[nh + seq_len(ng)])
    b <- beta0
    b[gam$where] <- v[nh + seq_len(ng)]
    res <- statmod_alternate(sp, design, blocks, hy,
                             iwls(tol = 1e-12, maxit = 500), b, FALSE, "opg",
                             500, 1e-12, verbosity(0))
    # the mode polished by Newton steps on the free coordinates: a reference
    # refitted to a stall is biased, not noisy
    par <- res$par
    free <- setdiff(seq_along(par), gam$where)
    for (it in 1:4) {
      gr <- res$obj$gr(par)
      he <- as.matrix(res$obj$he(par))
      par[free] <- par[free] - solve(he[free, free], gr[free])
    }
    cf <- res$obj$split(par)
    list(sp = sp, hy = hy, cf = cf,
         val = statmod_marginal(sp, design, cf, hy, method, "opg", basis)$value,
         g = statmod_marginal_grad(sp, design, cf, hy, method, idx, basis,
                                   gam = gam))
  }
  a <- at(v0)
  H <- statmod_marginal_hess(a$sp, design, a$cf, a$hy, method, idx, basis,
                             gam = gam)
  gap <- function(h) {
    fd_g <- vapply(seq_along(v0), function(j) {
      e <- replace(numeric(length(v0)), j, h)
      (at(v0 + e)$val - at(v0 - e)$val) / (2 * h)
    }, numeric(1))
    fd_H <- sapply(seq_along(v0), function(j) {
      e <- replace(numeric(length(v0)), j, h)
      (at(v0 + e)$g - at(v0 - e)$g) / (2 * h)
    })
    fd_H <- (fd_H + t(fd_H)) / 2
    c(g = max(abs(a$g - fd_g)) / max(abs(a$g)),
      H = max(abs(H - fd_H)) / max(abs(H)))
  }
  g1 <- gap(1e-2)
  g2 <- gap(1e-3)
  expect_lt(g2[["g"]], 1e-5)
  expect_lt(g2[["H"]], 1e-5)
  # O(h^2): a missing term would be flat in the step
  expect_gt(g1[["g"]] / g2[["g"]], 30)
  expect_gt(g1[["H"]] / g2[["H"]], 30)
})

test_that("marginal is checked against the family", {
  d <- marginal_data()
  g <- distributions7::gaussian1_distrib()
  expect_error(statmod(y ~ x, g, d, outer_criterion = reml(marginal = "nu")),
               "Its parameters are: mu, sigma")
  expect_error(reml(marginal = c("none", "sigma")), "alone")
  expect_error(reml(marginal = NA_character_), "'marginal' must be")
  d$w <- stats::runif(nrow(d))
  expect_error(statmod(y ~ x | sigma ~ te(z, w), g, d,
                       outer_criterion = reml(marginal = "sigma")),
               "null space")
})

test_that("a parameter the default names keeps the joint mode where not covered", {
  d <- marginal_data()
  d$w <- stats::runif(nrow(d))
  g <- distributions7::gaussian1_distrib()
  spec <- statmod_spec(y ~ x | sigma ~ te(z, w), g, d)
  design <- statmod_design(spec)
  mc <- marginal_coords(spec, design, reml())
  expect_length(mc$where, 0L)
  expect_match(mc$skipped[["sigma"]], "null space")
})

test_that("a block that is a working linearization keeps the joint mode", {
  # A sharp jseg()'s block is not a Jacobian, so a determinant over its
  # columns reads no curvature: measured, two fits reaching the same
  # break-point read -132.9 and -109.8 according to where they began.
  set.seed(2)
  d <- data.frame(x = sort(stats::runif(200, 0, 10)))
  d$y <- 1 + 0.4 * d$x + 1.5 * (d$x > 6) + stats::rnorm(200, sd = 0.3)
  g <- distributions7::gaussian1_distrib()
  spec <- statmod_spec(y ~ jseg(x), g, d)
  design <- statmod_design(spec)
  mc <- marginal_coords(spec, design, reml())
  expect_length(mc$where, 0L)
  expect_match(mc$skipped[["sigma"]], "working linearization")
  expect_error(statmod(y ~ jseg(x, n_boot = 0), g, d,
                       outer_criterion = reml(marginal = "sigma")),
               "working linearization")
})

test_that("an estimated coefficient at its limit is named a boundary", {
  # nu estimated on the criterion runs to its gaussian limit on gaussian
  # data, where the criterion's curvature in it vanishes: the certificate
  # names it by equation and coefficient, and gives it no hyperparameter key
  set.seed(11)
  d <- data.frame(x = stats::runif(300))
  d$y <- 1 + 2 * d$x + stats::rnorm(300, sd = 0.5)
  f <- statmod(y ~ x, distributions7::student_t1_distrib(), d)
  ct <- statmod_certificate(f)
  expect_identical(ct$state, "boundary")
  expect_identical(ct$boundary, "nu/(Intercept)")
  expect_length(ct$boundary_key, 0L)
})

# --- Verbyla (1993) and a kinked penalty beside the estimated coefficients --

verbyla_data <- function(n = 120, seed = 4, pz = 2) {
  set.seed(seed)
  d <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n))
  for (j in seq_len(pz)) d[[paste0("z", j)]] <- stats::rnorm(n)
  ls <- 0.1 + 0.4 * d$z1 - 0.3 * d$z2
  d$y <- 1 + d$x1 - d$x2 + exp(ls) * stats::rnorm(n)
  d
}

# the restricted log-likelihood of y ~ N(X beta, exp(2 Z g)), written out:
# beta by generalized least squares at g, and the determinant over X alone
verbyla_lr <- function(g, X, Z, y) {
  eta <- drop(Z %*% g)
  w <- exp(-2 * eta)
  A <- crossprod(X * w, X)
  beta <- solve(A, crossprod(X * w, y))
  r <- drop(y - X %*% beta)
  -sum(eta) - 0.5 * sum(w * r^2) -
    0.5 * as.numeric(determinant(A)$modulus) -
    (length(y) - ncol(X)) / 2 * log(2 * pi)
}

test_that("reml() on a modelled dispersion is the REML of Verbyla (1993)", {
  d <- verbyla_data()
  X <- cbind(1, d$x1, d$x2)
  Z <- cbind(1, d$z1, d$z2)
  f <- statmod(y ~ x1 + x2 | sigma ~ z1 + z2,
               distributions7::gaussian1_distrib(), d)
  g <- unname(f@coefficients$sigma)
  # the criterion the fit reports is the restricted log-likelihood at its g
  expect_equal(f@criterion, verbyla_lr(g, X, Z, d$y), tolerance = 1e-10)
  # and g maximizes it: the gradient of the written-out function vanishes
  # there, and its curvature is negative
  skip_if_not_installed("numDeriv")
  lr <- function(v) verbyla_lr(v, X, Z, d$y)
  expect_lt(max(abs(numDeriv::grad(lr, g))), 1e-4)
  expect_true(all(eigen(numDeriv::hessian(lr, g))$values < 0))
})

test_that("a kinked penalty no longer takes the dispersion off the criterion", {
  d <- verbyla_data(n = 200, seed = 7, pz = 6)
  g <- distributions7::gaussian1_distrib()
  fml <- y ~ x1 + x2 | sigma ~ lasso(~ z1 + z2 + z3 + z4 + z5 + z6, lambda = 20)
  spec <- statmod_spec(fml, g, d)
  design <- statmod_design(spec)
  mc <- marginal_coords(spec, design, reml())
  # the dispersion's intercept, and not its six penalized slopes
  expect_identical(mc$param, "sigma")
  expect_identical(mc$name, "(Intercept)")
  expect_length(mc$skipped, 0L)
  f <- statmod(fml, g, d)
  expect_true(f@converged)
  # a name the caller writes is accepted rather than refused
  f2 <- statmod(fml, g, d, outer_criterion = reml(marginal = "sigma"))
  expect_equal(f2@coefficients$sigma, f@coefficients$sigma, tolerance = 1e-8)
  # the intercept maximizes the criterion, the rest of the model refitted at
  # every value it is held at: the kinked slopes go back to the joint mode
  crit_at <- function(a) {
    sp <- statmod_hold(spec, mc, a)
    ds <- statmod_design(sp)
    cf <- fit_at_hyper_spec(sp, ds, f@hyper, g)
    statmod_marginal(sp, ds, cf, f@hyper, reml())$value
  }
  a0 <- f@coefficients$sigma[[1]]
  c0 <- crit_at(a0)
  expect_equal(c0, f@criterion, tolerance = 1e-8)
  for (h in c(-1e-2, 1e-2)) expect_lt(crit_at(a0 + h), c0)
})

test_that("the determinant leaves out every coordinate a kink covers", {
  d <- verbyla_data(n = 200, seed = 7, pz = 6)
  g <- distributions7::gaussian1_distrib()
  fml <- y ~ x1 + x2 | sigma ~ lasso(~ z1 + z2 + z3 + z4 + z5 + z6, lambda = 20)
  f <- statmod(fml, g, d, outer_criterion = reml(marginal = "none"))
  spec <- f@spec
  design <- statmod_design(spec)
  kc <- kinked_coords(spec, design)
  expect_length(kc, 6L)
  expect_true(all(kc %in% laplace_pinned(spec, design)))
  m <- statmod_marginal(spec, design, f@coefficients, f@hyper,
                        reml(marginal = "none"))
  nb <- length(unlist(f@coefficients))
  expect_identical(m$q, nb - length(kc))
  H <- statmod_information_at(spec, f@coefficients, design, FALSE, "opg")
  keep <- setdiff(seq_len(nb), kc)
  expect_equal(m$logdet,
               as.numeric(determinant(as.matrix(H)[keep, keep])$modulus),
               tolerance = 1e-10)
})

test_that("a coordinate a kink holds at zero adds nothing to the mode error", {
  # its entry of the smooth part's gradient is the likelihood's score, which
  # the kink's interval contains: read as a residual it put a fit at its mode
  # log-likelihood units above it
  d <- verbyla_data(n = 400, seed = 3, pz = 6)
  g <- distributions7::gaussian1_distrib()
  f <- statmod(y ~ x1 + x2 | sigma ~ lasso(~ z1 + z2 + z3 + z4 + z5 + z6,
                                           lambda = 60),
               g, d, outer_criterion = reml(marginal = "none"))
  spec <- f@spec
  design <- statmod_design(spec)
  zk <- zero_kinked(spec, design, f@coefficients)
  expect_gt(length(zk), 0L)
  obj <- statmod_objective(spec, f@hyper, design, FALSE, "opg")
  sc <- obj$gr(obj$stack(f@coefficients))
  expect_gt(max(abs(sc[zk])), 1)
  q <- inner_mode_error(NULL, spec, design, f@coefficients, f@hyper, sc)
  expect_lt(q, 1e-4)
  expect_identical(free_of_kinks(sc, spec, design, f@coefficients)[zk],
                   rep(0, length(zk)))
})

test_that("the exact outer gradient holds a kinked zero where it is", {
  skip_on_cran()
  # the mode's movement with a smoothing parameter leaves a coordinate at zero
  # at zero; solved as if it moved, the gradient was out by 3.3e-05 relative
  # and flat in the step
  set.seed(31)
  n <- 250
  d <- data.frame(x = stats::runif(n, -2, 2))
  for (j in 1:6) d[[paste0("z", j)]] <- stats::rnorm(n)
  d$y <- sin(1.4 * d$x) + 0.6 * d$z1 - 0.4 * d$z2 + stats::rnorm(n, sd = 0.3)
  g <- distributions7::gaussian1_distrib()
  fml <- y ~ s(x, bspline_smooth(k = 10)) +
    lasso(~ z1 + z2 + z3 + z4 + z5 + z6, lambda = 20)
  method <- reml(marginal = "none")
  fit0 <- statmod(fml, g, d, outer_criterion = NULL)
  spec <- fit0@spec
  design <- statmod_design(spec)
  idx <- outer_hyper_index(spec, statmod_blocks(spec, design))
  inner <- iwls(tol = 1e-11, maxit = 500)
  at <- function(eta) {
    hy <- eta_to_hyper(eta, idx, fit0@hyper)
    list(hy = hy, cf = fit_at_hyper(fml, g, d, hy, inner = inner)$coefficients)
  }
  eta <- hyper_to_eta(fit0@hyper, idx) + 0.4
  a <- at(eta)
  expect_gt(length(zero_kinked(spec, design, a$cf)), 0L)
  ex <- statmod_marginal_grad(spec, design, a$cf, a$hy, method, idx)
  crit <- function(e) {
    b <- at(e)
    statmod_marginal(spec, design, b$cf, b$hy, method)$value
  }
  h <- 1e-3
  fd <- (crit(eta + h) - crit(eta - h)) / (2 * h)
  expect_equal(ex, fd, tolerance = 1e-6)
})

test_that("a path point is scored where its nested search stands at the optimum", {
  skip_on_cran()
  # beside a kinked block the criterion carries the coordinate descent's own
  # rounding, and a nested search can run out of backtracks at the optimum:
  # on this model the path used to leave 10 of its 25 points unscored, and
  # one once the mode error read the free coordinates alone
  set.seed(4)
  n <- 500
  d <- data.frame(x = stats::runif(n, -2, 2))
  for (j in 1:8) d[[paste0("z", j)]] <- stats::rnorm(n)
  d$y <- sin(1.4 * d$x) + 0.5 * d$z1 - 0.4 * d$z2 + stats::rnorm(n, sd = 0.4)
  f <- suppressWarnings(statmod(
    y ~ s(x, bspline_smooth(k = 10)) +
      lasso(~ z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8) | sigma ~ 1,
    distributions7::gaussian1_distrib(), d,
    outer_criterion = reml(marginal = "none")))
  h <- as.data.frame(f@history$outer)
  h <- h[h$name == "lambda" & grepl("lasso", h$term), ]
  expect_identical(nrow(h), 25L)
  expect_identical(sum(is.na(h$criterion)), 0L)
  expect_true(f@converged)
})
