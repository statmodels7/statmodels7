# A kinked block beside a structural term is solved on the working response
# of the model. The reference is the derivative of the log-likelihood taken by
# numDeriv, which shares no arithmetic with the working response or with
# statmod_score_obs().

kkt_gap <- function(fit, lambda) {
  cf <- fit@coefficients
  des <- statmod_design(fit@spec)
  ll <- function(b) {
    c2 <- cf
    c2$mu <- b
    statmod_loglik_at(fit@spec, c2, des)
  }
  b <- cf$mu
  g <- numDeriv::grad(ll, b)
  act <- which(abs(b[-1L]) > 1e-10) + 1L
  zer <- setdiff(seq_along(b)[-1L], act)
  list(intercept = abs(g[1L]),
       active = if (length(act)) max(abs(g[act] - lambda * sign(b[act]))) else 0,
       zero = if (length(zer)) max(abs(g[zer])) / lambda else 0,
       b = b)
}

test_that("a kinked block beside a structural term reaches its KKT point", {
  skip_on_cran()
  # the data of the measurement in NEWS, drawn in the same order
  n <- 400L
  set.seed(11)
  z1 <- stats::rnorm(n)
  z2 <- stats::rnorm(n)
  z3 <- stats::rnorm(n)
  # a score-driven level, drawn inside the loop because the recursion is
  # driven by its own residual
  f <- numeric(n)
  y <- numeric(n)
  for (t in seq_len(n)) {
    y[t] <- 1 + 0.8 * z1[t] + f[t] + stats::rnorm(1, 0, 0.5)
    if (t < n) f[t + 1L] <- 0.3 * (y[t] - 1 - 0.8 * z1[t] - f[t]) + 0.7 * f[t]
  }
  dg <- data.frame(y = y, z1 = z1, z2 = z2, z3 = z3)
  s <- numeric(n)
  s[1L] <- 1
  for (t in 2:n) {
    s[t] <- if (stats::runif(1) < 0.9) s[t - 1L] else 3 - s[t - 1L]
  }
  yr <- 1 + 0.8 * z1 + ifelse(s == 2, 2, 0) + stats::rnorm(n, 0, 0.5)
  yp <- stats::rpois(n, exp(0.5 + 0.5 * z1 + ifelse(s == 2, 1, 0)))
  dp <- data.frame(y = yp, z1 = z1, z2 = z2, z3 = z3)

  # beside a filter: the level of gas() enters the predictor of mu and no
  # column of the lasso carries it. Built on the whole predictor, the
  # working response left an active score 38 from its KKT value.
  fg <- statmod(y ~ lasso(~ z1 + z2 + z3, lambda = 15) + gas(1, 1),
                distributions7::gaussian1_distrib(), dg)
  k <- kkt_gap(fg, 15)
  expect_lt(k$intercept, 1e-3)
  expect_lt(k$active, 1e-3)
  expect_lt(k$zero, 1)

  # beside a mixture over latent states, with a Poisson response: the score
  # at the posterior-mean predictor is not the posterior-weighted score, and
  # the fit diverged to an intercept of -309
  fp <- statmod(y ~ lasso(~ z1 + z2 + z3, lambda = 5) + regime(2),
                distributions7::poisson_distrib(), dp)
  k <- kkt_gap(fp, 5)
  expect_lt(k$intercept, 1e-3)
  expect_lt(k$active, 1e-3)
  expect_lt(abs(k$b[2L] - 0.5), 0.1)
})
