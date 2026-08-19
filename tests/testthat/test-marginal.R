# The marginal break-point terms through statmod(): the likelihood-shape
# branch regime() rides, with the mixture levels read off term_levels() --
# a vector of pattern shifts for the step kind, a per-observation matrix
# of node shifts for the continuous ones -- and the start off the exact
# profile of the predictor target. Everything here is plain maximum
# likelihood: the prior is part of the likelihood, so no outer criterion
# is involved.

marg_panel <- function(mI = 12L, nI = 20L, seed = 61) {
  set.seed(seed)
  id <- factor(rep(seq_len(mI), each = nI))
  x <- as.numeric(replicate(mI, sort(runif(nI, 0, 10))))
  psi <- 5 + rnorm(mI, 0, 0.5)
  y <- 1 + 2 * (x >= psi[as.integer(id)]) + rnorm(mI * nI, 0, 0.4)
  list(dd = data.frame(id = id, x = x, y = y), psi = psi)
}

test_that("a marginal jump fits end to end and reports what it estimated", {
  pn <- marg_panel()
  fit <- statmod(y ~ jump(x, psi ~ random(~1 | id), marginal = TRUE),
                 distributions7::gaussian1_distrib(), pn$dd)
  expect_true(fit@converged)
  st <- fit@structural[[1L]]$parameter
  expect_lt(abs(st[["m1"]] - 5), 0.6)
  expect_lt(abs(st[["tau1"]] - 0.5), 0.3)
  expect_lt(abs(st[["delta1"]] - 2), 0.3)
  # nothing is held: m does not enter the predictor as a level, so the
  # equation's intercept confounds nothing of the term's
  expect_identical(fit@structural[[1L]]$held, character(0))

  # the structural table carries (m, tau, delta) with standard errors from
  # the joint block, tau's interval built on the log scale and so positive
  s <- summary(fit)
  tab <- s@structural
  expect_identical(sort(tab$name), sort(c("m1", "tau1", "delta1")))
  expect_true(all(is.finite(tab$se)))
  expect_gt(tab$lower[tab$name == "tau1"], 0)

  V <- vcov(fit)
  expect_true(all(is.finite(diag(V))))

  # the latent posterior tracks the truth
  lat <- statmod_latent(fit)
  expect_identical(nrow(lat), 12L)
  expect_gt(cor(lat$mean, pn$psi), 0.75)

  # a model with no likelihood-shaped term has no latent to report
  f0 <- statmod(y ~ x, distributions7::gaussian1_distrib(), pn$dd)
  expect_error(statmod_latent(f0), "no structural term")
})

test_that("the layer's exact gradient matches numDeriv away from the optimum", {
  skip_if_not_installed("numDeriv")
  pn <- marg_panel(6L, 12L, seed = 3)
  spec <- statmod_spec(y ~ jump(x, psi ~ random(~1 | id), marginal = TRUE),
                       distributions7::gaussian1_distrib(), pn$dd)
  design <- statmod_design(spec)
  sst <- attr(design, "structure")
  key <- names(sst$zeta)[1L]
  z0 <- sst$zeta[[key]] + c(0.2, 0.1, -0.15)
  cf <- list(mu = 0.8, sigma = -0.7)
  obj <- function(u) {
    sst$zeta[[key]][] <- u[seq_len(3L)]
    sst$key <- NULL
    sst$value <- NULL
    statmod_loglik_at(spec, list(mu = u[4L], sigma = u[5L]), design)
  }
  u0 <- c(z0, cf$mu, cf$sigma)
  gn <- numDeriv::grad(obj, u0)
  sst$zeta[[key]][] <- z0
  sst$key <- NULL
  sst$value <- NULL
  gs <- statmod_structural_score(spec, cf, design)[[key]]
  gc <- statmod_score_at(spec, cf, design)
  gx <- c(gs, gc$mu, gc$sigma)
  expect_lt(max(abs(gx - gn)), 1e-6 * max(1, max(abs(gn))))
})

test_that("the fresh start reads the exact profile of the target", {
  pn <- marg_panel(8L, 16L, seed = 5)
  spec <- statmod_spec(y ~ jump(x, psi ~ random(~1 | id), marginal = TRUE),
                       distributions7::gaussian1_distrib(), pn$dd)
  design <- statmod_design(spec)
  z <- attr(design, "structure")$zeta[[1L]]
  # the profile puts m near the population position and the change of level
  # near the truth, where a conventional start has nothing to read them off
  expect_lt(abs(z[["m1"]] - 5), 1)
  expect_gt(z[["delta1"]], 1)
})

test_that("two latent break-points fit end to end", {
  set.seed(61)
  mI <- 8L
  nI <- 22L
  id <- factor(rep(seq_len(mI), each = nI))
  x <- as.numeric(replicate(mI, sort(runif(nI, 0, 10))))
  p1 <- 3 + rnorm(mI, 0, 0.4)
  p2 <- 7 + rnorm(mI, 0, 0.4)
  y <- 1 + 2 * (x >= p1[as.integer(id)]) - 1.5 * (x >= p2[as.integer(id)]) +
    rnorm(mI * nI, 0, 0.4)
  dd <- data.frame(id = id, x = x, y = y)
  fit <- statmod(y ~ jump(x, psi ~ random(~1 | id), npsi = 2,
                          marginal = TRUE),
                 distributions7::gaussian1_distrib(), dd)
  expect_true(fit@converged)
  st <- fit@structural[[1L]]$parameter
  expect_lt(abs(st[["m1"]] - 3), 0.7)
  expect_lt(abs(st[["m2"]] - 7), 0.7)
  expect_lt(abs(st[["delta1"]] - 2), 0.4)
  expect_lt(abs(st[["delta2"]] + 1.5), 0.4)
  lat <- statmod_latent(fit)
  expect_identical(nrow(lat), 16L)
  expect_gt(cor(lat$mean[lat$psi == 1L], p1), 0.5)
  expect_gt(cor(lat$mean[lat$psi == 2L], p2), 0.5)
})

test_that("the seg marginal fits end to end beside an intercept", {
  set.seed(61)
  mI <- 8L
  nI <- 16L
  id <- factor(rep(seq_len(mI), each = nI))
  x <- as.numeric(replicate(mI, sort(runif(nI, 0, 10))))
  ps <- 5 + rnorm(mI, 0, 0.5)
  y <- 1 + 0.5 * x - 1.2 * pmax(x - ps[as.integer(id)], 0) +
    rnorm(mI * nI, 0, 0.4)
  dd <- data.frame(id = id, x = x, y = y)
  fit <- statmod(y ~ seg(x, psi ~ random(~1 | id), marginal = TRUE),
                 distributions7::gaussian1_distrib(), dd)
  expect_true(fit@converged)
  st <- fit@structural[[1L]]$parameter
  expect_lt(abs(st[["beta"]] - 0.5), 0.2)
  expect_lt(abs(st[["m1"]] - 5), 0.7)
  expect_lt(abs(st[["gamma1"]] + 1.2), 0.3)
  lat <- statmod_latent(fit)
  expect_gt(cor(lat$mean, ps), 0.8)
})

test_that("a t prior rides the cdf surface end to end", {
  set.seed(17)
  mI <- 10L
  nI <- 14L
  id <- factor(rep(seq_len(mI), each = nI))
  x <- as.numeric(replicate(mI, sort(runif(nI, 0, 10))))
  ps <- pmin(pmax(5 + c(rnorm(8, 0, 0.3), 2.5, -2.8), 0.5), 9.5)
  y <- 1 + 2 * (x >= ps[as.integer(id)]) + rnorm(mI * nI, 0, 0.4)
  dd <- data.frame(id = id, x = x, y = y)
  pr <- distributions7::fixed(distributions7::student_t1_distrib(), mu = 0)
  fit <- statmod(y ~ jump(x, psi ~ random(~1 | id, distrib = pr),
                          marginal = TRUE),
                 distributions7::gaussian1_distrib(), dd)
  st <- fit@structural[[1L]]$parameter
  expect_identical(sort(names(st)),
                   sort(c("m1", "sigma", "nu", "delta1")))
  expect_lt(abs(st[["m1"]] - 5), 0.7)
  expect_lt(abs(st[["delta1"]] - 2), 0.4)
  lat <- statmod_latent(fit)
  # a heavy-tailed prior's edge moments can fail to exist, and what cannot
  # be computed is NA; most groups' posteriors sit inside the data
  expect_gt(mean(is.finite(lat$mean)), 0.7)
})

test_that("a kinked penalty fits beside the marginal term", {
  set.seed(3)
  mI <- 8L
  nI <- 14L
  id <- factor(rep(seq_len(mI), each = nI))
  x <- as.numeric(replicate(mI, sort(runif(nI, 0, 10))))
  ps <- 5 + rnorm(mI, 0, 0.5)
  z1 <- rnorm(mI * nI)
  z2 <- rnorm(mI * nI)
  z3 <- rnorm(mI * nI)
  y <- 1 + 2 * (x >= ps[as.integer(id)]) + 0.8 * z1 +
    rnorm(mI * nI, 0, 0.4)
  dd <- data.frame(id = id, x = x, y = y, z1 = z1, z2 = z2, z3 = z3)
  fit <- statmod(y ~ lasso(~ z1 + z2 + z3, n_lambda = 8) +
                   jump(x, psi ~ random(~1 | id), marginal = TRUE),
                 distributions7::gaussian1_distrib(), dd)
  st <- fit@structural[[1L]]$parameter
  expect_lt(abs(st[["m1"]] - 5), 0.7)
  expect_lt(abs(st[["delta1"]] - 2), 0.4)
  cf <- sort(abs(fit@coefficients$mu), decreasing = TRUE)
  # the intercept and the one signal column both survive with their size;
  # a selection that killed the signal would leave one large coefficient
  expect_gt(cf[2L], 0.4)
})
