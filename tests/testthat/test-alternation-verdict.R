# The alternation's verdict is read at the point it reached: the Newton
# decrement of the smooth part and the KKT decrement of each coefficient a
# kinked penalty holds at zero, both in log-likelihood units. Before, beside
# a filter it read the joint block's own flag, which the stall guard sets to
# FALSE exactly when that block is already at its mode.

verdict_data <- function() {
  n <- 400L
  set.seed(11)
  z1 <- stats::rnorm(n)
  z2 <- stats::rnorm(n)
  z3 <- stats::rnorm(n)
  f <- numeric(n)
  y <- numeric(n)
  for (t in seq_len(n)) {
    y[t] <- 1 + 0.8 * z1[t] + f[t] + stats::rnorm(1, 0, 0.5)
    if (t < n) f[t + 1L] <- 0.3 * (y[t] - 1 - 0.8 * z1[t] - f[t]) + 0.7 * f[t]
  }
  data.frame(y = y, z1 = z1, z2 = z2, z3 = z3)
}

# the design at the fitted structural state, which a fresh design does not
# carry
fitted_pieces <- function(fit) {
  sp <- fit@spec
  des <- statmod_design(sp)
  if (length(fit@structural)) {
    sst <- statmod_structural_state(des)
    for (k in names(fit@structural)) {
      sst$zeta[[k]] <- fit@structural[[k]]$unconstrained
      sst$held[[k]] <- fit@structural[[k]]$held
    }
  }
  obj <- statmod_objective(sp, fit@hyper, des, FALSE, "bartlett")
  list(spec = sp, des = des, obj = obj, b = obj$stack(fit@coefficients))
}

read_at <- function(p, fit, b = p$b) {
  alternation_readings(p$spec, p$des, p$obj, b, fit@hyper, FALSE, "bartlett")
}

test_that("a lasso beside a filter reports convergence at its KKT point", {
  skip_on_cran()
  fit <- statmod(y ~ lasso(~ z1 + z2 + z3, lambda = 15) + gas(1, 1),
                 distributions7::gaussian1_distrib(), verdict_data())
  expect_true(fit@converged)
  p <- fitted_pieces(fit)
  r <- read_at(p, fit)
  expect_lt(r[["mode"]], mode_error_limit())
  expect_lte(r[["kkt"]], mode_error_limit())
})

test_that("the readings refuse a point that is not the answer", {
  skip_on_cran()
  fit <- statmod(y ~ lasso(~ z1 + z2 + z3, lambda = 15) + gas(1, 1),
                 distributions7::gaussian1_distrib(), verdict_data())
  p <- fitted_pieces(fit)
  nm <- fit@spec@distrib@params
  j <- match("lasso.z1", unlist(lapply(p$des[nm], `[[`, "coef_names")))
  expect_false(is.na(j))
  expect_gt(abs(p$b[j]), 0.3)
  # the smooth part displaced: the intercept moved by 0.2
  b1 <- p$b
  b1[1L] <- b1[1L] + 0.2
  expect_gt(read_at(p, fit, b1)[["mode"]], 100 * mode_error_limit())
  # a coefficient with signal forced onto the kink: the data pull it off
  b2 <- p$b
  b2[j] <- 0
  expect_gt(read_at(p, fit, b2)[["kkt"]], 100 * mode_error_limit())
  expect_false(alternation_settled(p$spec, p$des, p$obj, b2, fit@hyper,
                                   FALSE, "bartlett"))
})

test_that("a mixture over latent states is read on the joint vector", {
  skip_on_cran()
  d <- verdict_data()
  set.seed(5)
  s <- numeric(nrow(d))
  s[1L] <- 1
  for (t in 2:nrow(d)) {
    s[t] <- if (stats::runif(1) < 0.9) s[t - 1L] else 3 - s[t - 1L]
  }
  d$y <- stats::rpois(nrow(d), exp(0.5 + 0.5 * d$z1 + ifelse(s == 2, 1, 0)))
  fit <- statmod(y ~ lasso(~ z1 + z2 + z3, lambda = 5) + regime(2),
                 distributions7::poisson_distrib(), d)
  expect_true(fit@converged)
  p <- fitted_pieces(fit)
  r <- read_at(p, fit)
  expect_lt(r[["mode"]], mode_error_limit())
  # the regime's own parameters are in the reading. The gap is displaced and
  # the coefficients re-solved at it, so that only the gradient in the
  # term's own parameters is left away from zero: a reading over the
  # coefficients alone would call that point settled.
  k <- names(fit@structural)[1L]
  sst <- statmod_structural_state(p$des)
  sst$zeta[[k]][["gap2"]] <- sst$zeta[[k]][["gap2"]] + 0.3
  jp <- statmod_joint_pieces(p$spec, p$des, p$obj, fit@hyper,
                             kinds = c("filter", "loglik"))
  b <- p$b
  act <- which(b != 0)
  for (it in 1:20) {
    u <- c(b, jp$zeta())
    g <- jp$gr(u)[act]
    K <- as.matrix(jp$he(u))[act, act]
    b[act] <- b[act] - solve(K, g)
  }
  u <- c(b, jp$zeta())
  expect_lt(max(abs(jp$gr(u)[act])), 1e-6)
  expect_gt(max(abs(jp$gr(u)[jp$ix])), 1)
  expect_gt(read_at(p, fit, b)[["mode"]], 100 * mode_error_limit())
})

test_that("either reading above the limit is enough to refuse", {
  lim <- mode_error_limit()
  settled <- function(mode, kkt) {
    testthat::local_mocked_bindings(
      alternation_readings = function(...) c(mode = mode, kkt = kkt))
    alternation_settled(NULL, NULL, NULL, NULL, NULL, FALSE, "bartlett")
  }
  expect_true(settled(0, NA))
  expect_true(settled(lim / 2, lim / 2))
  expect_false(settled(0, 10 * lim))
  expect_false(settled(10 * lim, 0))
  expect_false(settled(NA, 0))
})
