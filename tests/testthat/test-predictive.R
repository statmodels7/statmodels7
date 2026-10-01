gauss_fit <- function(n = 40, seed = 1) {
  set.seed(seed)
  dd <- data.frame(x = runif(n))
  dd$y <- 1 + 2 * dd$x + rnorm(n, sd = 0.5)
  statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
}

test_that("the plug-in interval is the family's quantiles at the estimates", {
  fit <- gauss_fit()
  nw <- data.frame(x = c(0.2, 0.8))
  p <- predict(fit, "response", nw, interval = "prediction",
               predictive = "plugin")
  mu <- predict(fit, "mu", nw)
  s <- predict(fit, "sigma", nw)
  expect_equal(p$lower, stats::qnorm(0.025, mu, s), tolerance = 1e-8)
  expect_equal(p$upper, stats::qnorm(0.975, mu, s), tolerance = 1e-8)
  expect_equal(p$se, s, tolerance = 1e-10)
})

test_that("the averaged interval adds se0^2 and sigma-hat^2 e^{2v}", {
  fit <- gauss_fit()
  nw <- data.frame(x = 0.7)
  p <- predict(fit, "response", nw, interval = "prediction")
  se0 <- predict(fit, "mu", nw, se = TRUE)$se
  V <- vcov(fit)
  v <- V["sigma:(Intercept)", "sigma:(Intercept)"]
  s <- predict(fit, "sigma", nw)
  expect_equal(p$se, sqrt(se0^2 + s^2 * exp(2 * v)), tolerance = 1e-6)
  # and it is wider than the plug-in one, which leaves both terms out
  q <- predict(fit, "response", nw, interval = "prediction",
               predictive = "plugin")
  expect_gt(p$upper - p$lower, q$upper - q$lower)
})

test_that("a bootstrap replica on the fit's own response is the fit", {
  fit <- gauss_fit()
  sp <- fit@spec
  m <- fit@methods$smooth
  cfg <- inner_settings(m)
  b0 <- unlist(fit@coefficients, use.names = FALSE)
  for (rf in c("coefficients", "full")) {
    r <- bootstrap_refit(fit, sp, rf, m, cfg, b0)
    expect_equal(unlist(r$coefficients), unlist(fit@coefficients),
                 tolerance = 1e-5, label = rf)
  }
})

test_that("the bootstrap gives an interval, reproducibly, for any family", {
  set.seed(5)
  x <- runif(50)
  dd <- data.frame(x = x, y = rpois(50, exp(1 + x)))
  fp <- statmod(y ~ x, distributions7::poisson_distrib(), dd)
  set.seed(9)
  a <- predict(fp, "response", data.frame(x = 0.5), interval = "prediction",
               predictive = "bootstrap", n_boot = 15)
  set.seed(9)
  b <- predict(fp, "response", data.frame(x = 0.5), interval = "prediction",
               predictive = "bootstrap", n_boot = 15)
  expect_identical(a, b)
  expect_true(a$lower <= a$fit && a$fit <= a$upper)
  # integer ends: the support is discrete
  expect_equal(c(a$lower, a$upper), round(c(a$lower, a$upper)))
})

test_that("the bootstrap redraws the effects set aside and keeps the rest", {
  set.seed(6)
  gg <- data.frame(g = factor(rep(1:10, each = 6)), x = rnorm(60))
  gg$y <- 1 + gg$x + rnorm(10, sd = 1)[gg$g] + rnorm(60, sd = 0.3)
  fit <- statmod(y ~ x + random(~ 1 | g), distributions7::gaussian1_distrib(),
                 gg)
  m_fit <- tapply(fit@fitted$mu, gg$g, mean)
  # set aside: the group means of a replica are not the fit's
  draw <- bootstrap_simulator(fit, fit@spec, statmod_design(fit@spec),
                              random_modes(fit@spec, "zero"))
  set.seed(1)
  expect_gt(stats::sd(tapply(draw(), gg$g, mean) - m_fit), 0.3)
  # conditional: the groups keep their effects, and only the residual moves
  rc <- random_modes(fit@spec, "conditional")
  draw_c <- bootstrap_simulator(fit, fit@spec, statmod_design(fit@spec),
                                rc[rc$mode != "conditional", , drop = FALSE])
  set.seed(1)
  expect_lt(stats::sd(tapply(draw_c(), gg$g, mean) - m_fit), 0.3)
})

test_that("a prior that is not Gaussian is drawn for an observed group", {
  set.seed(8)
  gg <- data.frame(g = factor(rep(1:12, each = 5)), x = rnorm(60))
  gg$y <- 1 + gg$x + 0.8 * rt(12, df = 3)[gg$g] + rnorm(60, sd = 0.4)
  tp <- distributions7::fixed(distributions7::student_t1_distrib(), mu = 0)
  fit <- suppressWarnings(statmod(y ~ x + random(~ 1 | g, distrib = tp),
                                  distributions7::gaussian1_distrib(), gg))
  old <- data.frame(x = 0, g = "3")
  set.seed(2)
  b <- predict(fit, "response", old, interval = "prediction",
               predictive = "bootstrap", n_boot = 10)
  a <- predict(fit, "response", old, interval = "prediction")
  expect_true(is.finite(b$lower) && is.finite(b$upper))
  # the observed group keeps its effect, so the two intervals are of one size
  expect_lt(abs(log((b$upper - b$lower) / (a$upper - a$lower))), 0.5)
})

test_that("a new group under a Student t prior is a scale mixture", {
  set.seed(8)
  gg <- data.frame(g = factor(rep(1:30, each = 5)), x = rnorm(150))
  gg$y <- 1 + gg$x + 0.8 * rt(30, df = 1.5)[gg$g] + rnorm(150, sd = 0.4)
  tp <- distributions7::fixed(distributions7::student_t1_distrib(), mu = 0)
  fit <- suppressWarnings(statmod(y ~ x + random(~ 1 | g, distrib = tp),
                                  distributions7::gaussian1_distrib(), gg))
  nw <- data.frame(x = 0.5, g = "new")
  p <- predict(fit, "response", nw, random = "zero", interval = "prediction",
               predictive = "plugin")
  h <- fit@hyper$mu[[1]]
  m <- predict(fit, "mu", nw, random = "zero")
  s <- predict(fit, "sigma", nw, random = "zero")
  # the mixture's distribution function at the reported ends, by adaptive
  # quadrature over the mixing variable, which shares nothing with the rule
  Fy <- function(y) stats::integrate(function(w)
    stats::pnorm((y - m) / sqrt(s^2 + h[["sigma"]]^2 / w)) *
      stats::dgamma(w, h[["nu"]] / 2, h[["nu"]] / 2), 0, Inf,
    rel.tol = 1e-10)$value
  # the rule in log w is exact to 1e-07 here; what remains, 4e-05 in
  # probability, is the 20-node Gauss-Hermite grid over the predictor,
  # whose spread under a small w is far wider than sigma
  expect_lt(abs(Fy(p$lower) - 0.025), 1e-4)
  expect_lt(abs(Fy(p$upper) - 0.975), 1e-4)
})

test_that("a new group under another prior is averaged by Monte Carlo", {
  set.seed(8)
  gg <- data.frame(g = factor(rep(1:12, each = 5)), x = rnorm(60))
  gg$y <- 1 + gg$x + rnorm(12)[gg$g] + rnorm(60, sd = 0.4)
  lp <- distributions7::fixed(distributions7::laplace_distrib(), mu = 0)
  fit <- suppressWarnings(statmod(y ~ x + random(~ 1 | g, distrib = lp),
                                  distributions7::gaussian1_distrib(), gg))
  nw <- data.frame(x = 0.5, g = "new")
  set.seed(1)
  p <- predict(fit, "response", nw, random = "zero", interval = "prediction",
               predictive = "plugin")
  sb <- fit@hyper$mu[[1]][["sigma"]]
  m <- predict(fit, "mu", nw, random = "zero")
  s <- predict(fit, "sigma", nw, random = "zero")
  set.seed(2)
  k <- 2e5
  yb <- m + stats::rexp(k, 1 / sb) * sample(c(-1, 1), k, TRUE) +
    stats::rnorm(k, 0, s)
  expect_equal(c(p$lower, p$upper),
               unname(stats::quantile(yb, c(0.025, 0.975))), tolerance = 0.03)
})

test_that("a prior with no variance reports none", {
  sp <- gauss_fit()@spec
  sp <- statmodels7:::spec_at(gauss_fit(), data.frame(x = 0.5),
                              need_response = FALSE)
  comps <- list(list(eta = list(mu = 1, sigma = 0), C = NULL, weight = 1))
  attr(comps, "infinite_variance") <- TRUE
  r <- predictive_response(sp, comps, 0.95)
  expect_true(is.na(r$se))
  expect_true(is.finite(r$lower))
})

test_that("predictive is refused where it is not read", {
  fit <- gauss_fit()
  expect_error(predict(fit, "mu", predictive = "plugin"), "read only")
  expect_error(predict(fit, "response", data.frame(x = 0.5),
                       interval = "prediction", predictive = "bootstrap",
                       n_boot = 0), "n_boot")
})
