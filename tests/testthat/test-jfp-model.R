# jfp() in statmod(): the penalized fit, its consumers, and the solve with a
# matrix plus a product of low rank.

jfp_sep_data <- function() {
  set.seed(11)
  n <- 60
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  # x1 separates the outcome quasi-completely
  d$y <- as.numeric(d$x1 > 0)
  d$y[1:6] <- 1 - d$y[1:6]
  d$y[d$x1 > 1] <- 1
  d
}

test_that("statmod() with jfp() over every coefficient is Firth's estimate", {
  skip_if_not_installed("numDeriv")
  d <- jfp_sep_data()
  f <- statmod(y ~ jfp(~ x1 + x2), distrib = distributions7::bernoulli_distrib(),
               data = d)
  b <- unlist(coef(f), use.names = FALSE)
  X <- cbind(1, d$x1, d$x2)
  closed <- function(b) {
    mu <- plogis(drop(X %*% b))
    sum(dbinom(d$y, 1, mu, log = TRUE)) +
      0.5 * as.numeric(determinant(crossprod(X * (mu * (1 - mu)), X))$modulus)
  }
  expect_lt(max(abs(numDeriv::grad(closed, b))), 1e-6)
  expect_true(all(is.finite(b)))
  expect_identical(statmod_certificate(f)$state, "converged")
  # the same estimate with the block written as a matrix
  d$X <- cbind(x1 = d$x1, x2 = d$x2)
  f2 <- statmod(y ~ jfp(X), distrib = distributions7::bernoulli_distrib(),
                data = d)
  expect_equal(unlist(coef(f2), use.names = FALSE), b, tolerance = 1e-8)
})

test_that("a factor in jfp() is coded with its contrasts", {
  set.seed(3)
  d <- data.frame(g = factor(rep(c("A", "B"), each = 80)), x = rnorm(160))
  lam <- exp(1.2 + 0.3 * d$x)
  d$y <- ifelse(runif(160) < 0.3, 0L, rpois(160, lam))
  # group B has no zeros: the zero-inflation of group B goes to zero
  d$y[d$g == "B" & d$y == 0] <- 1L
  d$gB <- as.numeric(d$g == "B")
  Z <- distributions7::zero_inflated(distributions7::poisson_distrib())
  f1 <- statmod(y ~ x | zi ~ jfp(~ g), distrib = Z, data = d)
  f2 <- statmod(y ~ x | zi ~ jfp(~ gB), distrib = Z, data = d)
  expect_equal(unlist(coef(f1), use.names = FALSE),
               unlist(coef(f2), use.names = FALSE), tolerance = 1e-8)
  expect_length(f1@aliased, 0L)
  expect_true(all(is.finite(unlist(coef(f1)))))
  expect_identical(statmod_certificate(f1)$state, "converged")
  # without the penalty the same coefficient runs away and is a boundary
  f0 <- statmod(y ~ x | zi ~ g, distrib = Z, data = d)
  expect_identical(statmod_certificate(f0)$state, "boundary")
})

test_that("jfp() beside a random intercept is fitted with reml()", {
  set.seed(4)
  d <- data.frame(id = factor(rep(1:20, each = 8)), x = rnorm(160),
                  z = rbinom(160, 1, 0.2))
  u <- rnorm(20, 0, 0.8)
  d$y <- rbinom(160, 1, plogis(-0.3 + 0.7 * d$x + u[d$id]))
  d$y[d$z == 1] <- 1L
  f <- statmod(y ~ x + jfp(~ z) + random(~ 1 | id),
               distrib = distributions7::bernoulli_distrib(), data = d)
  expect_true(is.finite(coef(f)$mu[["jfp.z"]]))
  expect_lt(coef(f)$mu[["jfp.z"]], 10)
  expect_identical(statmod_certificate(f)$state, "converged")
})

test_that("the consumers of a fit read the Jeffreys penalty", {
  d <- jfp_sep_data()
  f <- statmod(y ~ x2 + jfp(~ x1), distrib = distributions7::bernoulli_distrib(),
               data = d)
  for (vt in list(bayesian_vcov(), frequentist_vcov(), unconditional_vcov(),
                  frequentist_vcov(sandwich = hc0_sandwich()))) {
    V <- vcov(f, vcov_type = vt)
    expect_true(all(is.finite(diag(V))))
  }
  # the bayesian variance is the inverse of the penalized curvature
  sp <- f@spec
  des <- statmod_design(sp)
  K <- as_dense(statmod_information_at(sp, f@coefficients, des, TRUE)) +
    as_dense(statmod_penalty_at(sp, f@coefficients, f@hyper, des, "hessian"))
  expect_equal(unname(vcov(f, readable = FALSE)), unname(solve(K)),
               tolerance = 1e-6)
  s <- summary(f)
  expect_true(S7::S7_inherits(statmod_test(f, param = "mu",
                                           coefname = "jfp.x1",
                                           test = lr_test()), StatmodTest))
  pr <- statmod_profile(f, which = "mu:jfp.x1")
  expect_true(all(is.finite(pr@points$profile_penalty)))
  sim <- simulate(f, nsim = 1, seed = 1)
  expect_equal(nrow(sim), nrow(d))
  e <- f@edf
  expect_true(all(is.finite(e$edf)))
  expect_true(is.finite(AIC(f)))
})

test_that("the effective degrees of freedom can exceed p under separation", {
  d <- jfp_sep_data()
  f <- statmod(y ~ jfp(~ x1 + x2), distrib = distributions7::bernoulli_distrib(),
               data = d)
  expect_true(is.finite(sum(f@edf$edf)))
})

test_that("a block whose information is singular has no usable value", {
  d <- jfp_sep_data()
  d$z0 <- 0
  sp <- statmod_spec(y ~ x1 + z0, distributions7::bernoulli_distrib(), d)
  des <- statmod_design(sp)
  # a block on a column of zeros: the information is exactly zero, so the
  # factorization fails on every platform (a repeated column left a tiny
  # positive pivot on Linux's BLAS)
  st <- jfp_state(sp, des, list(mu = c(0.1, 0.2, 0)),
                  list(jfp_block("mu", 3L)))
  expect_identical(jfp_value(st), -Inf)
  expect_error(jfp_gradient(st), "not positive definite")
})

test_that("lowrank_solve() equals the dense solve", {
  set.seed(2)
  p <- 30
  M <- crossprod(matrix(rnorm(60 * p), 60, p))
  V <- matrix(rnorm(p * 3), p, 3)
  u <- rnorm(p)
  ref <- solve(M + tcrossprod(V) + 0.5 * diag(p), u)
  expect_equal(lowrank_solve(M, V, u, damp = 0.5), ref, tolerance = 1e-10)
  Ms <- Matrix::Matrix(M, sparse = TRUE)
  expect_equal(lowrank_solve(Ms, V, u, damp = 0.5), ref, tolerance = 1e-10)
  # an indefinite A has no factor and the caller assembles the sum
  expect_null(lowrank_solve(M - 1e4 * diag(p), V, u))
  pc <- list(A = M, V = V)
  expect_equal(iwls_solve(pc, u, "chol")$delta, solve(M + tcrossprod(V), u),
               tolerance = 1e-10)
})

test_that("the third and fourth derivatives match differences of the lower", {
  skip_if_not_installed("numDeriv")
  set.seed(7)
  n <- 80
  d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:8, 10)))
  d$y <- rbinom(n, 1, plogis(0.3 + d$x - 0.5 * d$z))
  sp <- statmod_spec(y ~ x + z + random(~1 | g),
                     distributions7::bernoulli_distrib(), d)
  des <- statmod_design(sp)
  m <- sum(vapply(des, function(z) z$npar, integer(1)))
  bl <- list(jfp_block("mu", 2:3, const = TRUE))
  set.seed(3)
  u0 <- rnorm(m, sd = 0.2)
  v <- rnorm(m)
  w <- rnorm(m)
  split <- function(u) list(mu = u)
  st_at <- function(u) jfp_state(sp, des, split(u), bl)
  n3 <- matrix(numDeriv::jacobian(function(t)
    as.vector(jfp_hessian(st_at(u0 + t * v))$hessian), 0), m)
  n4 <- matrix(numDeriv::jacobian(function(t)
    as.vector(jfp_deriv3_dir(st_at(u0 + t * w), v)), 0), m)
  st <- st_at(u0)
  D3 <- jfp_deriv3_dir(st, v)
  D4 <- jfp_deriv4_dir(st, v, w)
  expect_lt(max(abs(D3 - n3)) / max(abs(n3)), 1e-7)
  expect_lt(max(abs(D4 - n4)) / max(abs(n4)), 1e-7)
  M <- crossprod(matrix(rnorm(m * m), m)) / m
  expect_equal(jfp_deriv3(st, M, v), sum(M * D3), tolerance = 1e-10)
  expect_equal(jfp_deriv4(st, M, v, w), sum(M * D4), tolerance = 1e-10)
  # the negative control: the low-rank terms carry a share of D3
  t3 <- jfp_terms3(st, v)
  t3$lr <- list()
  expect_gt(max(abs(jfp_terms_matrix(st, t3) - n3)) / max(abs(n3)), 1e-2)
})

test_that("the outer gradient and Hessian of reml() and ml() are exact", {
  skip_if_not_installed("numDeriv")
  set.seed(4)
  d <- data.frame(id = factor(rep(1:20, each = 8)), x = rnorm(160),
                  z = rbinom(160, 1, 0.2))
  d$y <- rbinom(160, 1, plogis(-0.3 + 0.7 * d$x +
                                 rnorm(20, 0, 0.8)[d$id]))
  d$y[d$z == 1] <- 1L
  fm <- y ~ x + jfp(~ z) + random(~1 | id)
  B <- distributions7::bernoulli_distrib()
  for (method in list(reml(information = observed_info()),
                      ml(information = observed_info()))) {
    fit0 <- statmod(fm, distrib = B, data = d, outer_criterion = NULL)
    spec <- fit0@spec
    design <- statmod_design(spec)
    idx <- outer_hyper_index(spec, statmod_blocks(spec, design))
    expect_true(outer_gradient_ok(spec, design, idx, method, 2L))
    basis <- integrated_basis(spec, design, method@kind)
    at <- function(eta) {
      hy <- eta_to_hyper(eta, idx, fit0@hyper)
      list(hy = hy, cf = fit_at_hyper(fm, B, d, hy,
                                      polish = TRUE)$coefficients)
    }
    fn <- function(eta) {
      a <- at(eta)
      statmod_marginal(spec, design, a$cf, a$hy, method, basis = basis)$value
    }
    gr <- function(eta) {
      a <- at(eta)
      statmod_marginal_grad(spec, design, a$cf, a$hy, method, idx, basis)
    }
    eta <- hyper_to_eta(fit0@hyper, idx) + 0.3
    a <- at(eta)
    g <- gr(eta)
    gn <- (fn(eta + 1e-4) - fn(eta - 1e-4)) / 2e-4
    expect_lt(abs(g - gn) / abs(gn), 1e-7)
    H <- statmod_marginal_hess(spec, design, a$cf, a$hy, method, idx, basis)
    Hn <- (gr(eta + 1e-4) - gr(eta - 1e-4)) / 2e-4
    expect_lt(abs(H - Hn) / abs(Hn), 1e-7)
  }
  # a fit under ml() counts the jfp coefficient among the estimated ones
  f <- statmod(fm, distrib = B, data = d, outer_criterion = ml())
  expect_identical(statmod_certificate(f)$state, "converged")
  lab <- coef_labels(f@spec, statmod_design(f@spec))
  expect_equal(attr(logLik(f, loglik_type = "marginal"), "df"),
               sum(!lab$penalized) + 1 + 1)
})

test_that("a block in a subformula of nl() moves, and its derivatives are exact", {
  skip_if_not_installed("numDeriv")
  set.seed(2)
  n <- 150
  d <- data.frame(x = runif(n, 0.2, 3), g = rnorm(n))
  d$y <- rpois(n, exp(log(3 + 0.5 * d$g) - 0.6 * d$x))
  fm <- y ~ 0 + nl(~ exp(a - r * x), a ~ 1 + jfp(~ g))
  sp <- statmod_spec(fm, distributions7::poisson_distrib(), d)
  des <- statmod_design(sp)
  units <- statmod_penalized(sp, des)
  expect_identical(units[[1]]$kind, "refresh")
  expect_true(units[[1]]$const)
  bl <- jfp_blocks(sp, des)
  u0 <- c(1.1, 0.15, 0.55)
  ref <- function(u) {
    cf <- list(mu = u)
    th <- statmod_eta(sp, des, cf)$theta
    X <- as.matrix(statmod_design_at(sp, cf, des)$mu$X)[, bl[[1]]$basis] %*%
      bl[[1]]$C
    E <- distributions7::distrib_expected_hessian(sp@distrib, sp@response, th,
                                                  scale = "link")
    0.5 * as.numeric(determinant(crossprod(X * (-E$mu_mu), X))$modulus)
  }
  st_at <- function(u) jfp_state(sp, des, list(mu = u), bl)
  g <- jfp_gradient(st_at(u0))
  gn <- numDeriv::grad(ref, u0)
  expect_lt(max(abs(g - gn)) / max(abs(gn)), 1e-7)
  hj <- numDeriv::jacobian(function(u) jfp_gradient(st_at(u)), u0)
  h <- jfp_hessian(st_at(u0))$hessian
  expect_lt(max(abs(h - (hj + t(hj)) / 2)) / max(abs(hj)), 1e-7)
  # the negative control: the same block read as though its rows were fixed
  st0 <- st_at(u0)
  st0$blocks[[1]]$block$kind <- "static"
  expect_gt(max(abs(jfp_gradient(st0) - gn)) / max(abs(gn)), 1e-2)
})

test_that("jfp() inside nl() keeps a separated slope finite", {
  set.seed(6)
  n <- 200
  d <- data.frame(x = runif(n, -2, 2), h = rbinom(n, 1, 0.3))
  d$y <- rbinom(n, 1, plogis(2 * (d$x - 0.3)))
  d$y[d$h == 1] <- as.numeric(d$x[d$h == 1] > 0.5)
  B <- distributions7::bernoulli_distrib()
  f <- statmod(y ~ 0 + nl(~ b * (x - e), b ~ 1 + jfp(~ h),
                          start = list(b = 1, e = 0)), distrib = B, data = d)
  cf <- unlist(coef(f))
  expect_true(all(is.finite(cf)))
  expect_lt(abs(cf[["mu.nl.b.jfp.h"]]), 50)
  expect_identical(statmod_certificate(f)$state, "converged")
  obj <- statmod_objective(f@spec, f@hyper, statmod_design(f@spec))
  expect_lt(max(abs(obj$gr(obj$stack(f@coefficients)))), 1e-3)
})

test_that("the subformulas that move discontinuously are rejected", {
  set.seed(6)
  n <- 120
  d <- data.frame(x = runif(n, -2, 2), h = rbinom(n, 1, 0.3))
  d$y <- 1 + 0.5 * d$x + 1.5 * pmax(d$x - 0.2 - 0.6 * d$h, 0) +
    rnorm(n, sd = 0.3)
  G <- distributions7::gaussian1_distrib()
  expect_error(statmod(y ~ seg(x, psi ~ 1 + jfp(~ h)), distrib = G, data = d),
               "sharp break-point")
  f <- suppressWarnings(statmod(y ~ seg(x, psi ~ 1 + jfp(~ h),
                                         smoothed = numericals7::smooth_probit()),
                                distrib = G, data = d))
  expect_true(all(is.finite(unlist(coef(f)))))
})

test_that("the outer criterion reaches the gaussian closed form beside jfp()", {
  # with the full information the Jeffreys value is
  # const - 2 a p_tilde log sigma, so the REML optimum is
  # sigma^2 = RSS / (n - p + 2 a p_tilde)
  set.seed(1)
  n <- 200
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  d$y <- 1 + 0.5 * d$x + 0.3 * d$z + rnorm(n)
  f <- statmod(y ~ x + jfp(~ z) | sigma ~ 1,
               distrib = distributions7::gaussian1_distrib(), data = d)
  rss <- sum(stats::residuals(stats::lm(y ~ x + z, d))^2)
  an <- sqrt(rss / (n - 3 + 2 * 0.5 * 2))
  expect_equal(exp(f@coefficients$sigma[[1L]]), an, tolerance = 1e-7)
  expect_identical(statmod_certificate(f)$state, "converged")
  # the gradient at the optimum, and without the value's derivative (the
  # defect): 2 a p_tilde
  spec <- f@spec
  design <- statmod_design(spec)
  idx <- outer_hyper_index(spec, statmod_blocks(spec, design))
  gam <- marginal_coords(spec, design, reml())
  basis <- integrated_basis(spec, design, "reml", gamma = TRUE)
  # on the specification that holds sigma, as the outer search reads it
  sp <- statmod_hold(spec, gam, f@coefficients$sigma)
  g <- statmod_marginal_grad(sp, design, f@coefficients, f@hyper, reml(),
                             idx, basis, gam = gam)
  jg <- jfp_penalty(spec, design, f@coefficients, "gradient")
  expect_lt(abs(g), 1e-5)
  expect_equal(g + jg[gam$where], 2, tolerance = 1e-6)
})

test_that("a jfp block keeps its parameter on the outer criterion", {
  set.seed(1)
  n <- 100
  d <- data.frame(x = rnorm(n), w = rnorm(n))
  d$y <- 1 + 0.5 * d$x + rnorm(n, sd = exp(0.3 * d$w))
  G <- distributions7::gaussian1_distrib()
  route <- function(fm) {
    sp <- statmod_spec(fm, G, d)
    des <- statmod_design(sp)
    marginal_coords(sp, des, reml())
  }
  a <- route(y ~ x | sigma ~ jfp(~ w))
  b <- route(y ~ x | sigma ~ w)
  # the same parameter on the criterion, nothing skipped; the block's own
  # coordinate stays at the mode
  expect_identical(unique(a$param), unique(b$param))
  expect_length(a$skipped, 0L)
  expect_identical(a$name, "(Intercept)")
  expect_identical(b$name, c("(Intercept)", "w"))
  f <- statmod(y ~ x | sigma ~ jfp(~ w), distrib = G, data = d)
  expect_identical(statmod_certificate(f)$state, "converged")
})

test_that("a model without jfp() does not enumerate the penalized units", {
  set.seed(4)
  d <- data.frame(x = stats::runif(40, 0.2, 3), g = stats::rnorm(40),
                  z = stats::rbinom(40, 1, 0.4))
  d$y <- stats::rpois(40, exp(1 - 0.4 * d$x + 0.3 * d$g))
  P <- distributions7::poisson_distrib()
  s0 <- statmod_spec(y ~ s(x) + g, P, d)
  expect_false(s0@jeffreys)
  expect_true(statmod_spec(y ~ s(x) + jfp(~ g), P, d)@jeffreys)
  # in the subformula of a term, through term_penalties()
  s2 <- statmod_spec(y ~ 0 + nl(~ exp(a - r * x), a ~ 1 + jfp(~ g)), P, d)
  expect_true(s2@jeffreys)
  # where the property was not set the terms are read
  expect_true(spec_jeffreys(S7::set_props(s2, jeffreys = NA)))
  expect_false(spec_jeffreys(S7::set_props(s0, jeffreys = NA)))
  # with the property FALSE the Jeffreys state is NULL without enumerating
  local_mocked_bindings(statmod_penalized = function(...) stop("enumerated"))
  ds <- statmod_design(s0)
  expect_null(jfp_state_at(s0, ds, list(mu = rep(0, ds$mu$npar))))
  expect_identical(jfp_blocks(s0, ds), list())
  expect_error(jfp_blocks(s2, statmod_design(s2)), "enumerated")
})
