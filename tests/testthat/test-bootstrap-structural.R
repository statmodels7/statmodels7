# The bootstrap of a model with a break-point term (replicates refitted with
# the restart cycle) and of a model with a gas() or regime() term (replicates
# that are new paths of the series). The references are written in the test:
# a brute-force bootstrap with the same draws, and the recursion of the model.

sim_gas_bs <- function(n, omega, kappa1, b1, sd = 1, seed = 21) {
  set.seed(seed)
  f <- numeric(n)
  y <- numeric(n)
  s <- 0
  f0 <- omega / (1 - b1)
  for (t in seq_len(n)) {
    f[t] <- omega + kappa1 * (if (t > 1) s else 0) +
      b1 * (if (t > 1) f[t - 1] else f0)
    y[t] <- f[t] + stats::rnorm(1, sd = sd)
    s <- (y[t] - f[t]) / sd^2
  }
  data.frame(t = seq_len(n), y = y, x = stats::rnorm(n))
}

seg_data_bs <- function(n = 150, seed = 3) {
  set.seed(seed)
  x <- sort(stats::runif(n, 0, 10))
  y <- 1 + 0.5 * x + 1.5 * pmax(x - 6, 0) + stats::rnorm(n, sd = 0.4)
  data.frame(x = x, y = y)
}

test_that("the case bootstrap of a seg() fit equals a brute-force bootstrap with the same draws", {
  dd <- seg_data_bs()
  n <- nrow(dd)
  fit <- statmod(y ~ seg(x), distributions7::gaussian1_distrib(), dd)
  B <- 8L
  set.seed(42)
  t_eng <- system.time(
    fit <- add_bootstrap(fit, n = B, resample = case_resample()))[["elapsed"]]
  eng <- coef(statmod_bootstrap(fit, resample = case_resample()))
  expect_identical(dim(eng), c(B, 5L))

  # the same stream: the counts of a replicate, then the seed of its refit
  # (the engine draws every replicate before it refits any)
  set.seed(42)
  draws <- lapply(seq_len(B), function(b) {
    list(w = tabulate(sample.int(n, n, TRUE), n),
         seed = sample.int(.Machine$integer.max, 1L))
  })
  brute <- t(vapply(draws, function(d) {
    set.seed(d$seed)
    w <- d$w
    f <- statmod(y ~ seg(x), distributions7::gaussian1_distrib(), dd,
                 weights = w)
    unlist(coef(f, readable = FALSE), use.names = FALSE)
  }, numeric(5)))
  # the position and the other coefficients agree replicate by replicate
  expect_lt(max(abs(eng - brute)), 0.05)
  psi <- grep("psi1", colnames(eng))
  expect_equal(stats::sd(eng[, psi]), stats::sd(brute[, psi]), tolerance = 0.1)
  # the engine's variance of the position is the covariance of its replicates
  V <- vcov(fit, vcov_type = bootstrap_vcov(), readable = FALSE)
  expect_equal(unname(V[psi, psi]), stats::var(eng[, psi]), tolerance = 1e-12)
  expect_lt(t_eng / B, 5)
})

test_that("the restart cycle rescues a replicate whose refit from the fit's position stalls", {
  n <- 200
  set.seed(7)
  x <- sort(stats::runif(n, 0, 10))
  set.seed(1)
  d1 <- data.frame(x = x, y = 1 + 0.2 * x + 3 * pmax(x - 3, 0) +
                     stats::rnorm(n, sd = 0.3))
  f <- statmod(y ~ seg(x), distributions7::gaussian1_distrib(), d1)
  spec0 <- f@spec
  cfg <- inner_settings(f@methods$smooth)
  beta0 <- unlist(f@coefficients, use.names = FALSE)
  # the replicate has a small change of slope where the fit has its
  # break-point and a large one at 8
  set.seed(104)
  y2 <- 1 + 0.2 * x + 0.4 * pmax(x - 3, 0) - 4 * pmax(x - 8, 0) +
    stats::rnorm(n, sd = 0.3)
  sb <- S7::set_props(spec0, response = y2)
  plain <- bootstrap_refit(f, sb, "coefficients", f@methods$smooth, cfg, beta0)
  set.seed(4)
  rescued <- boot_refit_breakpoint(f, sb, "coefficients", f@methods$smooth,
                                   cfg, beta0)
  fresh <- statmod(y ~ seg(x), distributions7::gaussian1_distrib(),
                   data.frame(x = x, y = y2))
  psi_fresh <- coef(fresh, readable = FALSE)$mu[[4L]]
  expect_gt(abs(plain$coefficients$mu[[4L]] - psi_fresh), 1)
  expect_equal(rescued$coefficients$mu[[4L]], psi_fresh, tolerance = 0.01)
})

test_that("a gas() replicate is a new path of the filter, reproducible and extendable", {
  dd <- sim_gas_bs(200, 0.3, 0.4, 0.7)
  fit <- statmod(y ~ x + gas(p = 1, q = 1, time = t) - 1,
                 distributions7::gaussian1_distrib(), dd)
  spec0 <- fit@spec
  design0 <- statmod_design(spec0)

  # one draw: a new path whose innovations, read through the fitted filter,
  # are independent normal draws of the fitted scale
  set.seed(9)
  draw <- bootstrap_plan(fit, parametric_resample(), spec0, design0,
                         random_modes(spec0, "zero"), TRUE)
  y1 <- draw()$y
  y2 <- draw()$y
  expect_length(y1, 200L)
  expect_lt(abs(stats::cor(y1, dd$y)), 0.3)
  expect_lt(abs(stats::cor(y1, y2)), 0.3)
  ev <- statmod_eta(S7::set_props(spec0, response = y1), design0,
                    fit@coefficients)
  z <- y1 - ev$eta$mu
  sigma <- exp(fit@coefficients$sigma[[1L]])
  expect_equal(stats::sd(z), sigma, tolerance = 0.15)
  expect_lt(abs(stats::cor(z[-1L], z[-200L])), 0.2)

  # the same seed gives the same replicates, and 3 are the first 3 of 6
  run <- function(n, seed = NULL) {
    set.seed(5)
    bootstrap_engine(fit, n, parametric_resample(), "coefficients")
  }
  r6 <- run(6)
  r6b <- run(6)
  r3 <- run(3)
  expect_identical(r6$coefficients, r6b$coefficients)
  expect_identical(r3$coefficients, r6$coefficients[1:3, ])
  expect_identical(colnames(r6$coefficients),
                   c("mu:x", "sigma:(Intercept)", "mu:gas.omega",
                     "mu:gas.kappa1", "mu:gas.pacf1"))

  # an entry extended from 3 to 6 equals the 6 of one run
  set.seed(5)
  fit <- add_bootstrap(fit, n = 3, resample = parametric_resample())
  fit <- add_bootstrap(fit, n = 6, resample = parametric_resample())
  e6 <- coef(statmod_bootstrap(fit, resample = parametric_resample()))
  expect_equal(unname(e6), unname(r6$coefficients), tolerance = 1e-12)

  # the variance is the covariance of the replicates, by equation
  V <- vcov(fit, vcov_type = bootstrap_vcov(), readable = FALSE)
  expect_identical(dim(V), c(5L, 5L))
  expect_equal(V["mu:x", "mu:gas.omega"], stats::cov(e6)["mu:x", "mu:gas.omega"],
               tolerance = 1e-12)
  expect_identical(rownames(V)[1:4], c("mu:x", "mu:gas.omega", "mu:gas.kappa1",
                                       "mu:gas.pacf1"))
})

test_that("case and cluster resamples are rejected beside a structural term", {
  dd <- sim_gas_bs(80, 0.3, 0.4, 0.7)
  dd$id <- 1
  fit <- statmod(y ~ x + gas(p = 1, q = 1, time = t) - 1,
                 distributions7::gaussian1_distrib(), dd)
  expect_error(add_bootstrap(fit, n = 3, resample = case_resample()),
               "parametric_resample\\(\\)")
  expect_error(add_bootstrap(fit, n = 3, resample = cluster_resample(~ id)),
               "parametric_resample\\(\\)")
})

test_that("a regime() replicate is a new path and the bootstrap runs", {
  set.seed(3)
  n <- 150
  P <- matrix(c(0.93, 0.07, 0.06, 0.94), 2, 2, byrow = TRUE)
  s <- integer(n)
  s[1L] <- 1L
  for (t in 2:n) s[t] <- sample.int(2L, 1L, prob = P[s[t - 1L], ])
  dd <- data.frame(t = seq_len(n), y = stats::rnorm(n, c(0, 3)[s]))
  fit <- statmod(y ~ regime(2, time = t) - 1,
                 distributions7::gaussian1_distrib(), dd)
  spec0 <- fit@spec
  design0 <- statmod_design(spec0)
  set.seed(2)
  draw <- bootstrap_plan(fit, parametric_resample(), spec0, design0,
                         random_modes(spec0, "zero"), TRUE)
  y1 <- draw()$y
  expect_lt(abs(stats::cor(y1, dd$y)), 0.5)
  expect_gt(stats::sd(y1), 1.1)
  run <- function(n) {
    set.seed(8)
    bootstrap_engine(fit, n, parametric_resample(), "coefficients")
  }
  expect_identical(run(4)$coefficients, run(6)$coefficients[1:4, ])
  e <- run(6)
  expect_identical(ncol(e$coefficients), 5L)
  expect_true(all(is.finite(e$coefficients)))
  fit <- { set.seed(8); add_bootstrap(fit, n = 6) }
  V <- vcov(fit, vcov_type = bootstrap_vcov())
  expect_true(all(c("mu:regime.p1.1", "mu:regime.p2.2") %in% rownames(V)))
  expect_true(all(is.finite(diag(V))))
})

test_that("the replicates of a sharp jump() carry the position of the break-point", {
  set.seed(3)
  n <- 150
  x <- sort(stats::runif(n, 0, 10))
  y <- 1 + 2 * (x > 6) + 0.3 * x + stats::rnorm(n, sd = 0.4)
  dd <- data.frame(x = x, y = y)
  fit <- statmod(y ~ jump(x), distributions7::gaussian1_distrib(), dd)
  tm <- fit@spec@terms$mu[[2L]]
  held <- fit@coefficients$mu[[3L]]
  pos <- modelterms7::term_readable(tm, fit@coefficients$mu[2:3])$value[[2L]]
  expect_lt(abs(pos - 6), 0.2)
  set.seed(1)
  fit <- add_bootstrap(fit, n = 6, resample = parametric_resample())
  B <- coef(statmod_bootstrap(fit, resample = parametric_resample()))
  # the held slot of the fit is a working value; the replicates store the
  # position instead
  expect_true(all(abs(B[, "mu:jump.g1"] - pos) < 0.5))
  expect_gt(abs(held - pos), 1)
})
