# The HAC variance. The references are written in the test: the textbook
# Newey-West estimator of a linear model, and per-observation scores of a
# score-driven or a regime model taken by central differences of the
# per-observation log-likelihood.

sim_gas_hac <- function(n, omega, kappa1, b1, seed = 21) {
  set.seed(seed)
  f <- numeric(n)
  y <- numeric(n)
  s <- 0
  f0 <- omega / (1 - b1)
  for (t in seq_len(n)) {
    f[t] <- omega + kappa1 * (if (t > 1) s else 0) +
      b1 * (if (t > 1) f[t - 1] else f0)
    y[t] <- f[t] + stats::rnorm(1)
    s <- y[t] - f[t]
  }
  data.frame(t = seq_len(n), y = y, x = stats::rnorm(n))
}

# the textbook meat: Gamma_0 + sum_l k(l/(L+1)) (Gamma_l + Gamma_l')
nw_meat <- function(u, L, kernel = "bartlett") {
  k <- switch(kernel,
    bartlett = function(x) 1 - x,
    parzen = function(x) ifelse(x <= 0.5, 1 - 6 * x^2 + 6 * x^3,
                                2 * (1 - x)^3),
    truncated = function(x) 1)
  n <- nrow(u)
  S <- crossprod(u)
  for (l in seq_len(L)) {
    G <- crossprod(u[(l + 1):n, , drop = FALSE], u[1:(n - l), , drop = FALSE])
    S <- S + k(l / (L + 1)) * (G + t(G))
  }
  S
}

test_that("HAC reduces to the textbook Newey-West estimator of a linear model", {
  set.seed(1)
  n <- 150
  d <- data.frame(t = sample(n), x = stats::rnorm(n))
  d <- d[order(d$t), ]
  e <- as.numeric(stats::filter(stats::rnorm(n), 0.6, method = "recursive"))
  d$y <- 1 + 0.5 * d$x + e
  d <- d[sample(n), ]   # the rows are shuffled: the order comes from t
  fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), d)
  o <- order(d$t)
  X <- cbind(1, d$x)[o, ]
  res <- stats::residuals(stats::lm(y ~ x, d))[o]
  XtXi <- solve(crossprod(X))
  L0 <- floor(4 * (n / 100)^(2 / 9))
  for (kern in c("bartlett", "parzen", "truncated")) {
    for (L in c(L0, 7L)) {
      Vnw <- XtXi %*% nw_meat(X * res, L, kern) %*% XtXi
      V <- vcov(fit, vcov_type = frequentist_vcov(
        sandwich = hac_sandwich(kernel = kern, bandwidth = L, time = ~ t)),
        readable = FALSE)
      expect_equal(unname(V[1:2, 1:2]), unname(Vnw), tolerance = 1e-8)
    }
  }
  # the default bandwidth is floor(4 (T/100)^(2/9))
  expect_equal(
    vcov(fit, vcov_type = frequentist_vcov(sandwich = hac_sandwich(
      time = ~ t)), readable = FALSE),
    vcov(fit, vcov_type = frequentist_vcov(sandwich = hac_sandwich(
      bandwidth = L0, time = ~ t)), readable = FALSE))
})

test_that("HAC at bandwidth zero is HC0, and close to it without serial correlation", {
  set.seed(2)
  n <- 500
  d <- data.frame(t = seq_len(n), x = stats::rnorm(n))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(n, sd = 1 + abs(d$x))
  fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), d)
  v0 <- vcov(fit, vcov_type = frequentist_vcov(sandwich = hc0_sandwich()),
             readable = FALSE)
  h0 <- vcov(fit, vcov_type = frequentist_vcov(
    sandwich = hac_sandwich(bandwidth = 0, time = ~ t)), readable = FALSE)
  expect_equal(h0, v0, tolerance = 1e-12)
  h <- vcov(fit, vcov_type = frequentist_vcov(
    sandwich = hac_sandwich(time = ~ t)), readable = FALSE)
  expect_equal(sqrt(diag(h)), sqrt(diag(v0)), tolerance = 0.15)
})

test_that("series are separated by the cluster variable in a model without a structural term", {
  set.seed(3)
  G <- 12
  m <- 15
  d <- data.frame(id = rep(seq_len(G), each = m), t = rep(seq_len(m), G),
                  x = stats::rnorm(G * m))
  e <- unlist(lapply(seq_len(G), function(g)
    as.numeric(stats::filter(stats::rnorm(m), 0.5, method = "recursive"))))
  d$y <- 1 + d$x + e
  fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), d)
  X <- cbind(1, d$x)
  res <- stats::residuals(stats::lm(y ~ x, d))
  u <- X * res
  L <- 3L
  S <- matrix(0, 2, 2)
  for (g in seq_len(G)) S <- S + nw_meat(u[d$id == g, ], L)
  XtXi <- solve(crossprod(X))
  Vref <- XtXi %*% S %*% XtXi
  V <- vcov(fit, vcov_type = frequentist_vcov(sandwich = hac_sandwich(
    bandwidth = L, time = ~ t, cluster = ~ id)), readable = FALSE)
  expect_equal(unname(V[1:2, 1:2]), unname(Vref), tolerance = 1e-8)
  # not separating the series changes the answer
  V1 <- vcov(fit, vcov_type = frequentist_vcov(sandwich = hac_sandwich(
    bandwidth = L, time = ~ t)), readable = FALSE)
  expect_gt(max(abs(V1 - V)), 1e-6)
})

test_that("the scores of a gas() fit are the derivatives of the observation's log-likelihood", {
  dd <- sim_gas_hac(120, 0.3, 0.4, 0.7)
  fit <- statmod(y ~ x + gas(p = 1, q = 1, time = t) - 1,
                 distributions7::gaussian1_distrib(), dd)
  spec <- fit@spec
  design <- statmod_design(spec)
  sst <- statmod_structural_state(design)
  tn <- names(sst$zeta)
  cf0 <- fit@coefficients
  z0 <- sst$zeta[[tn]]
  # a point off the optimum, where the scores are not small
  u1 <- c(unlist(cf0, use.names = FALSE), unlist(z0)) +
    c(0.1, 0.1, 0.15, -0.2, 0.3)
  place <- function(u) {
    cf <- cf0
    cf$mu[] <- u[1]
    cf$sigma[] <- u[2]
    sst$zeta[[tn]][] <- u[3:5]
    cf
  }
  ll_t <- function(u) {
    cf <- place(u)
    ev <- statmod_eta(spec, design, cf)
    as.numeric(distributions7::distrib_pdf(spec@distrib, spec@response,
                                           ev$theta, log = TRUE))
  }
  Un <- vapply(1:5, function(j) {
    h <- 1e-6
    up <- um <- u1
    up[j] <- up[j] + h
    um[j] <- um[j] - h
    (ll_t(up) - ll_t(um)) / (2 * h)
  }, numeric(nrow(dd)))
  U <- joint_scores(spec, place(u1), design)
  expect_lt(max(abs(U - Un)) / max(abs(Un)), 1e-7)

  # at the fit the scores sum to zero; the HAC variance is the sandwich of
  # the meat of those scores, with the textbook weights
  ll_fit <- function(j, h = 1e-6) {
    u0 <- c(unlist(cf0, use.names = FALSE), unlist(z0))
    up <- um <- u0
    up[j] <- up[j] + h
    um[j] <- um[j] - h
    (ll_t(up) - ll_t(um)) / (2 * h)
  }
  Uf <- vapply(1:5, ll_fit, numeric(nrow(dd)))
  expect_lt(max(abs(colSums(Uf))), 1e-5)
  Vb <- vcov(fit, vcov_type = bayesian_vcov(), readable = FALSE)
  # joint order: the coefficients, then the filter's parameters; the
  # variance is reported by equation, mu first
  ord <- c("mu:x", "sigma:(Intercept)", "mu:gas.omega", "mu:gas.kappa1",
           "mu:gas.pacf1")
  Hb <- Vb[ord, ord]
  L <- floor(4 * (120 / 100)^(2 / 9))
  Vref <- Hb %*% nw_meat(Uf, L) %*% Hb
  Vh <- vcov(fit, vcov_type = frequentist_vcov(sandwich = hac_sandwich()),
             readable = FALSE)
  expect_equal(unname(Vh[ord, ord]), unname(Vref), tolerance = 1e-5)
  # the summary reports the term's own parameters with the same variance
  expect_no_error(summary(fit, vcov_type = frequentist_vcov(
    sandwich = hac_sandwich())))
})

test_that("lags run within a group in a gas() panel", {
  set.seed(4)
  G <- 6
  m <- 40
  d <- data.frame(id = rep(seq_len(G), each = m), t = rep(seq_len(m), G),
                  x = stats::rnorm(G * m))
  d$y <- 0.5 * d$x + unlist(lapply(seq_len(G), function(g)
    as.numeric(stats::filter(stats::rnorm(m), 0.4, method = "recursive"))))
  fit <- statmod(y ~ x + gas(p = 1, q = 1, by = id, time = t) - 1,
                 distributions7::gaussian1_distrib(), d)
  spec <- fit@spec
  design <- statmod_design(spec)
  U <- joint_scores(spec, fit@coefficients, design)
  series <- hac_series(fit, design, hac_sandwich())
  expect_length(series, G)
  expect_identical(sort(unlist(series)), seq_len(G * m))
  hm <- hac_meat(U, series, "bartlett", 3L)$meat
  ref <- Reduce(`+`, lapply(series, function(ix) nw_meat(U[ix, ], 3L)))
  expect_equal(hm, ref, tolerance = 1e-12)
  # the whole sample as one series would add the products across groups
  one <- hac_meat(U, list(seq_len(G * m)), "bartlett", 3L)$meat
  expect_gt(max(abs(one - hm)), 1e-8)
})

test_that("the scores of a regime() fit are the derivatives of the increments of its likelihood", {
  set.seed(3)
  n <- 120
  P <- matrix(c(0.93, 0.07, 0.06, 0.94), 2, 2, byrow = TRUE)
  s <- integer(n)
  s[1L] <- 1L
  for (t in 2:n) s[t] <- sample.int(2L, 1L, prob = P[s[t - 1L], ])
  dd <- data.frame(t = seq_len(n), x = stats::rnorm(n))
  dd$y <- stats::rnorm(n, c(0, 3)[s]) + 0.5 * dd$x
  fit <- statmod(y ~ x + regime(2, time = t) - 1,
                 distributions7::gaussian1_distrib(), dd)
  spec <- fit@spec
  design <- statmod_design(spec)
  sst <- statmod_structural_state(design)
  tn <- names(sst$zeta)
  cf0 <- fit@coefficients
  z0 <- unlist(sst$zeta[[tn]])
  place <- function(u) {
    cf <- cf0
    cf$mu[] <- u[1]
    cf$sigma[] <- u[2]
    sst$zeta[[tn]][] <- u[3:6]
    sst$key <- NULL
    cf
  }
  inc <- function(u) {
    cf <- place(u)
    ev <- statmod_eta(spec, design, cf)
    r <- ev$regimes[[1L]]
    modelterms7::term_loglik(r$tm, r$eta_static, spec@response,
                             r$cb$logdens, r$cb$score, r$psi)$loglik
  }
  u1 <- c(0.4, 0.1, z0 + c(0.2, 0.1, -0.3, 0.2))
  Un <- vapply(seq_along(u1), function(j) {
    h <- 1e-6
    up <- um <- u1
    up[j] <- up[j] + h
    um[j] <- um[j] - h
    (inc(up) - inc(um)) / (2 * h)
  }, numeric(n))
  U <- joint_scores(spec, place(u1), design)
  expect_lt(max(abs(U - Un)) / max(abs(Un)), 1e-7)
  # and the HAC variance of the fit runs, with the group order of the term
  Vh <- vcov(fit, vcov_type = frequentist_vcov(sandwich = hac_sandwich()),
             readable = FALSE)
  expect_true(all(is.finite(diag(Vh))))
  expect_identical(dim(Vh), c(6L, 6L))
})

test_that("the HAC variance rejects what it cannot cover", {
  set.seed(5)
  n <- 100
  d <- data.frame(t = seq_len(n), x = sort(stats::runif(n, 0, 10)))
  d$y <- 1 + 0.5 * d$x + 1.5 * pmax(d$x - 6, 0) + stats::rnorm(n, sd = 0.4)
  fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), d)
  expect_error(vcov(fit, vcov_type = frequentist_vcov(
    sandwich = hac_sandwich())), "time = ~ t")
  fs <- statmod(y ~ seg(x), distributions7::gaussian1_distrib(), d)
  # a term whose block moves is covered, see test-sandwich-refresh.R
  Vs <- vcov(fs, vcov_type = frequentist_vcov(
    sandwich = hac_sandwich(time = ~ t)), readable = FALSE)
  expect_true(all(is.finite(diag(Vs))))
  dg <- sim_gas_hac(60, 0.3, 0.4, 0.7)
  fg <- statmod(y ~ x + gas(p = 1, q = 1, time = t) - 1,
                distributions7::gaussian1_distrib(), dg)
  expect_error(vcov(fg, vcov_type = frequentist_vcov(
    sandwich = hac_sandwich(time = ~ t))), "structural term")
  expect_error(vcov(fg, vcov_type = frequentist_vcov(
    sandwich = hc0_sandwich())), "structural term")
  expect_error(hac_sandwich(kernel = "gauss"), "should be one of")
  expect_error(hac_sandwich(bandwidth = -1), "non-negative")
  expect_error(hac_sandwich(bandwidth = 1.5), "non-negative")
})
