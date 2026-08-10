# What a fit reports about its own uncertainty.

set.seed(4)
n <- 400
dd <- data.frame(x = runif(n), z = runif(n),
                 g = factor(rep(c("a", "b"), length.out = n)))

test_that("an unpenalized gaussian fit reproduces the closed-form variance", {
  # the reference is the linear model's own, which shares no code with the
  # penalized information assembled here
  sim <- rstatmod(y ~ x + g, distributions7::gaussian1_distrib(), dd,
                  par = list(mu = c(1, 2, -0.5), sigma = log(0.4)))
  fit <- statmod(y ~ x + g, distributions7::gaussian1_distrib(), sim)
  V <- vcov(fit)

  lm_fit <- stats::lm(y ~ x + g, sim)
  s2 <- sum(stats::residuals(lm_fit)^2) / n   # the MLE, not the unbiased one
  V_lm <- s2 * solve(crossprod(stats::model.matrix(lm_fit)))

  idx <- grep("^mu:", rownames(V))
  expect_equal(unname(V[idx, idx]), unname(V_lm), tolerance = 1e-6)

  # and the scale block: eta = log(sigma), so the information is 2n and the
  # standard error 1/sqrt(2n) whatever sigma is
  expect_equal(unname(sqrt(V["sigma:(Intercept)", "sigma:(Intercept)"])),
               1 / sqrt(2 * n), tolerance = 1e-6)
})

test_that("with no penalty the two conventions are the same matrix", {
  sim <- rstatmod(y ~ x, distributions7::gaussian1_distrib(), dd,
                  par = list(mu = c(1, 2), sigma = log(0.4)))
  fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), sim)
  expect_equal(vcov(fit, "bayesian"), vcov(fit, "frequentist"),
               tolerance = 1e-10)
})

test_that("a ridge makes the bayesian variance the larger of the two", {
  sim <- rstatmod(y ~ ridge(~ x + z), distributions7::gaussian1_distrib(), dd,
                  par = list(mu = c(1, 2, -1), sigma = log(0.5)))
  fit <- statmod(y ~ ridge(~ x + z), distributions7::gaussian1_distrib(), sim)
  vb <- diag(vcov(fit, "bayesian"))
  vf <- diag(vcov(fit, "frequentist"))
  expect_true(all(vb >= vf - 1e-12))
  expect_true(any(vb > vf * 1.000001))
  # V_f = V_b H V_b is the sandwich the definition states, computed here from
  # the pieces rather than read back from the function under test
  H <- -hessian(fit, expected = TRUE)
  Vb <- vcov(fit, "bayesian")
  expect_equal(Vb %*% H %*% Vb, vcov(fit, "frequentist"), tolerance = 1e-8)
})

test_that("vcov agrees with a numerical Hessian of the objective", {
  # the assembled information is a sum of crossprods; numDeriv on the whole
  # penalized objective shares none of that arithmetic
  skip_if_not_installed("numDeriv")
  sim <- rstatmod(y ~ x | sigma ~ z, distributions7::gaussian1_distrib(), dd,
                  par = list(mu = c(1, 2), sigma = c(-1, 0.8)))
  fit <- statmod(y ~ x | sigma ~ z, distributions7::gaussian1_distrib(), sim)
  spec <- fit@spec
  design <- statmod_design(spec)
  obj <- statmod_objective(spec, fit@hyper, design, expected = FALSE)
  b <- unlist(fit@coefficients[spec@distrib@params], use.names = FALSE)
  Hn <- numDeriv::hessian(obj$fn, b)
  expect_equal(unname(vcov(fit, expected = FALSE)), solve(Hn),
               tolerance = 1e-5)
})

test_that("a confidence interval is symmetric and covers the truth", {
  sim <- rstatmod(y ~ x, distributions7::gaussian1_distrib(), dd,
                  par = list(mu = c(1, 2), sigma = log(0.4)))
  fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), sim)
  ci <- confint(fit)
  expect_named(ci, c("parameter", "term", "coefficient", "estimate", "se",
                     "lower", "upper"))
  # a coefficient of a linear predictor is free on the whole line, so the
  # interval is symmetric about the estimate on the scale it is built on
  expect_equal(ci$upper - ci$estimate, ci$estimate - ci$lower,
               tolerance = 1e-12)
  # 400 observations: the truth is inside
  expect_true(ci["mu:x", "lower"] < 2 && 2 < ci["mu:x", "upper"])
  expect_true(ci["mu:(Intercept)", "lower"] < 1 &&
                1 < ci["mu:(Intercept)", "upper"])
  # the width follows the level
  wide <- confint(fit, level = 0.99)
  expect_true(all(wide$upper - wide$lower > ci$upper - ci$lower))
})

test_that("confint selects by parameter, by term and by label", {
  sim <- rstatmod(y ~ x | sigma ~ z, distributions7::gaussian1_distrib(), dd,
                  par = list(mu = c(1, 2), sigma = c(-1, 0.8)))
  fit <- statmod(y ~ x | sigma ~ z, distributions7::gaussian1_distrib(), sim)
  expect_identical(unique(confint(fit, "sigma")$parameter), "sigma")
  expect_identical(nrow(confint(fit, "mu:x")), 1L)
  expect_error(confint(fit, "nonesuch"), "matched nothing")
  expect_error(confint(fit, level = 1), "strictly between")
})

test_that("a coefficient a lasso set to zero has no variance", {
  # at the kink the penalty is not twice differentiable, so there is no
  # curvature to read and NA is the reading rather than a number
  sim <- dd
  sim$y <- 1 + 2 * dd$x + stats::rnorm(n, sd = 0.4)
  sim$noise1 <- stats::runif(n)
  sim$noise2 <- stats::runif(n)
  fit <- statmod(y ~ x + lasso(~ noise1 + noise2),
                 distributions7::gaussian1_distrib(), sim,
                 hyper = list(mu = list(lasso = c(lambda = 200))))
  b <- fit@coefficients$mu
  V <- vcov(fit)
  zero <- which(b == 0)
  expect_gt(length(zero), 0)
  expect_true(all(is.na(diag(V)[zero])))
  # while the unpenalized ones are perfectly available
  expect_true(all(is.finite(diag(V)[grep("^mu:x$", rownames(V))])))
  s <- summary(fit)
  expect_true(any(grepl("kinked penalty", s@notes)))
})

test_that("a summary carries the tables, the criteria and the notes", {
  sim <- rstatmod(y ~ x | sigma ~ z, distributions7::gaussian1_distrib(), dd,
                  par = list(mu = c(1, 2), sigma = c(-1, 0.8)))
  fit <- statmod(y ~ x | sigma ~ z, distributions7::gaussian1_distrib(), sim)
  s <- summary(fit)
  expect_named(s@tables, c("mu", "sigma"))
  expect_identical(nrow(s@tables$mu), 2L)
  expect_true(all(c("estimate", "se", "statistic", "p_value", "lower",
                    "upper") %in% names(s@tables$mu)))
  # the criteria are built on the same log-likelihood and df logLik reports,
  # so AIC() and the summary cannot disagree
  expect_equal(s@aic, stats::AIC(fit), tolerance = 1e-12)
  expect_equal(s@bic, stats::BIC(fit), tolerance = 1e-12)
  expect_equal(s@loglik, as.numeric(stats::logLik(fit)), tolerance = 1e-12)
  # a slope of 2 on 400 observations is not in doubt
  expect_lt(s@tables$mu$p_value[2L], 1e-10)
  expect_output(print(s), "log-likelihood")
  expect_output(print(s), "Distribution")
})

test_that("a penalized term's rows say they are conditional", {
  sim <- rstatmod(y ~ ridge(~ x + z), distributions7::gaussian1_distrib(), dd,
                  par = list(mu = c(1, 2, -1), sigma = log(0.5)))
  fit <- statmod(y ~ ridge(~ x + z), distributions7::gaussian1_distrib(), sim)
  s <- summary(fit)
  expect_true(any(s@tables$mu$penalized))
  expect_true(any(grepl("smoothing parameter", s@notes)))
})

test_that("a fit whose information is singular says so rather than guessing", {
  # two identical columns: the data does not identify their difference, and a
  # pseudo-inverse would report a standard error for a direction it cannot see
  sim <- rstatmod(y ~ x, distributions7::gaussian1_distrib(), dd,
                  par = list(mu = c(1, 2), sigma = log(0.4)))
  sim$xcopy <- sim$x
  expect_error(
    vcov(statmod(y ~ x + xcopy, distributions7::gaussian1_distrib(), sim)),
    "not positive definite")
})
