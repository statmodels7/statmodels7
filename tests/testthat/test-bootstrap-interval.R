g1 <- distributions7::gaussian1_distrib()

quiet <- function(e) suppressWarnings(suppressMessages(e))

eng0 <- bootstrap_engine

# the limits that the interval objects must give, written from the stored
# replicates with quantile() alone
ref_limits <- function(fit, n, est, level, basic = FALSE, cols = NULL,
                       resample = parametric_resample()) {
  B <- coef(statmod_bootstrap(fit, resample = resample))
  B <- B[seq_len(n), , drop = FALSE]
  if (!is.null(cols)) B <- B[, cols, drop = FALSE]
  a <- 1 - level
  q <- apply(B, 2L, stats::quantile, probs = c(a / 2, 1 - a / 2),
             names = FALSE)
  if (basic) list(lower = 2 * est - q[2L, ], upper = 2 * est - q[1L, ])
  else list(lower = q[1L, ], upper = q[2L, ])
}

pois_fit <- function() {
  set.seed(2)
  d <- data.frame(x = rnorm(150))
  d$y <- rpois(150, exp(0.5 + 0.4 * d$x))
  statmod(y ~ x, distributions7::poisson_distrib(), d)
}

test_that("the percentile interval is the quantiles of the stored replicates", {
  fit <- pois_fit()
  set.seed(3)
  fit <- add_bootstrap(fit, n = 200)
  ci <- confint(fit, level = 0.9, interval = percentile_interval(n = 200))
  est <- unlist(fit@coefficients, use.names = FALSE)
  ref <- ref_limits(fit, 200, est, 0.9)
  expect_equal(ci$lower, unname(ref$lower), tolerance = 1e-12)
  expect_equal(ci$upper, unname(ref$upper), tolerance = 1e-12)
  expect_equal(ci$estimate, est)
  # the standard error is the bootstrap one
  B <- coef(statmod_bootstrap(fit, resample = parametric_resample()))
  expect_equal(ci$se, unname(apply(B, 2L, stats::sd)), tolerance = 1e-12)
  # readable = FALSE gives the same table
  ci2 <- confint(fit, level = 0.9, readable = FALSE,
                 interval = percentile_interval(n = 200))
  expect_equal(ci2$lower, ci$lower)
  # arguments NULL read the single stored entry
  ci3 <- confint(fit, level = 0.9, interval = percentile_interval())
  expect_equal(ci3$lower, ci$lower)
})

test_that("the basic interval reflects the percentile interval about the estimate", {
  fit <- pois_fit()
  set.seed(3)
  fit <- add_bootstrap(fit, n = 200)
  est <- unlist(fit@coefficients, use.names = FALSE)
  pc <- confint(fit, level = 0.9, interval = percentile_interval())
  bs <- confint(fit, level = 0.9, interval = basic_interval())
  expect_equal(bs$lower, 2 * est - pc$upper, tolerance = 1e-12)
  expect_equal(bs$upper, 2 * est - pc$lower, tolerance = 1e-12)
  ref <- ref_limits(fit, 200, est, 0.9, basic = TRUE)
  expect_equal(bs$lower, unname(ref$lower), tolerance = 1e-12)
  expect_equal(bs$upper, unname(ref$upper), tolerance = 1e-12)
  expect_true(all(bs$lower < bs$estimate & bs$estimate < bs$upper))
})

test_that("a Gaussian with a modeled scale gives the quantiles on the link scale", {
  set.seed(1)
  d <- data.frame(x = runif(120))
  d$y <- 1 + 2 * d$x + rnorm(120, sd = exp(-0.5 + d$x))
  fit <- statmod(y ~ x | sigma ~ x, g1, d)
  set.seed(2)
  fit <- add_bootstrap(fit, n = 200)
  ci <- confint(fit, level = 0.9, interval = percentile_interval())
  est <- unlist(fit@coefficients, use.names = FALSE)
  # the readable quantity of a linear predictor's coefficient is the
  # coefficient itself, on the link scale of its parameter
  expect_equal(ci$estimate, est)
  ref <- ref_limits(fit, 200, est, 0.9)
  expect_equal(ci$lower, unname(ref$lower), tolerance = 1e-12)
  expect_equal(ci$upper, unname(ref$upper), tolerance = 1e-12)
  expect_equal(rownames(ci), c("mu:(Intercept)", "mu:x", "sigma:(Intercept)",
                               "sigma:x"))
  # parm selects the rows
  cs <- confint(fit, parm = "sigma", level = 0.9,
                interval = basic_interval())
  expect_equal(rownames(cs), c("sigma:(Intercept)", "sigma:x"))
  bs <- ref_limits(fit, 200, est, 0.9, basic = TRUE)
  expect_equal(cs$lower, unname(bs$lower[3:4]), tolerance = 1e-12)
})

test_that("coordinates without replicates have missing limits", {
  set.seed(4)
  d <- data.frame(x = rnorm(60))
  d$x2 <- 2 * d$x
  d$y <- 1 + d$x + rnorm(60)
  fit <- quiet(statmod(y ~ x + x2, g1, d))
  expect_true(length(fit@aliased) > 0L)
  set.seed(5)
  fit <- quiet(add_bootstrap(fit, n = 40))
  ci <- quiet(confint(fit, level = 0.5, interval = percentile_interval()))
  al <- rownames(ci) %in% fit@aliased
  expect_true(all(is.na(ci$lower[al])) && all(is.na(ci$upper[al])))
  expect_true(all(is.finite(ci$lower[!al])))
  # a coordinate missing in some replicates, as a random effect under the
  # cluster scheme is, is missing in its limits and leaves the others alone
  fit2 <- pois_fit()
  set.seed(6)
  fit2 <- add_bootstrap(fit2, n = 40)
  env <- statmodels7:::boot_env(fit2)
  key <- ls(env)[[1L]]
  e <- get(key, envir = env)
  M <- e@coefficients
  M[3L, 2L] <- NA_real_
  e@coefficients <- M
  assign(key, e, envir = env)
  ci2 <- quiet(confint(fit2, level = 0.5, interval = percentile_interval()))
  expect_true(is.na(ci2$lower[[2L]]) && is.na(ci2$upper[[2L]]))
  expect_true(is.finite(ci2$lower[[1L]]))
})

test_that("n beyond the stored replicates computes only the difference", {
  fitA <- pois_fit()
  set.seed(5)
  fitA <- add_bootstrap(fitA, n = 30, resample = case_resample())
  set.seed(99)
  expect_message(
    ciA <- suppressWarnings(confint(
      fitA, level = 0.5,
      interval = percentile_interval(n = 50, resample = case_resample()))),
    "computed 20 replicates")
  fitB <- pois_fit()
  set.seed(5)
  fitB <- add_bootstrap(fitB, n = 50, resample = case_resample())
  cA <- coef(statmod_bootstrap(fitA, resample = case_resample()))
  cB <- coef(statmod_bootstrap(fitB, resample = case_resample()))
  expect_identical(cA, cB)
  ciB <- confint(fitB, level = 0.5,
                 interval = percentile_interval(n = 50,
                                                resample = case_resample()))
  expect_identical(ciA$lower, ciB$lower)
  expect_identical(ciA$upper, ciB$upper)
  # a smaller n reads the first n replicates and computes nothing
  expect_no_message(ci30 <- suppressWarnings(confint(
    fitA, level = 0.5,
    interval = percentile_interval(n = 30, resample = case_resample()))))
  ref <- ref_limits(fitA, 30, unlist(fitA@coefficients, use.names = FALSE),
                    level = 0.5, resample = case_resample())
  expect_equal(ci30$lower, unname(ref$lower), tolerance = 1e-12)
})


test_that("a tail with fewer than ten replicates warns", {
  fit <- pois_fit()
  set.seed(3)
  fit <- add_bootstrap(fit, n = 40)
  expect_warning(confint(fit, interval = percentile_interval()),
                 "40 replicates leave 1.0 in each tail")
  expect_no_warning(confint(fit, level = 0.5,
                            interval = percentile_interval()))
})

test_that("the argument interval is checked and the default is unchanged", {
  fit <- pois_fit()
  set.seed(3)
  fit <- add_bootstrap(fit, n = 20)
  expect_error(confint(fit, interval = "percentile"),
               "must be NULL or a bootstrap_interval")
  expect_error(confint(fit, interval = percentile_interval(),
                       vcov_type = frequentist_vcov()),
               "vcov_type")
  expect_error(confint(fit, readable = FALSE, interval = basic_interval(),
                       test = lr_test()), "cannot be combined")
  expect_error(percentile_interval(n = -1), "positive whole number")
  expect_error(basic_interval(refit = "x"), "refit")
  expect_error(bootstrap_interval(), "abstract")
  expect_true(S7::S7_inherits(basic_interval(), bootstrap_interval))
  # the default is the Wald interval of the Bayesian variance
  ci <- confint(fit, level = 0.9)
  V <- vcov(fit)
  est <- unlist(fit@coefficients, use.names = FALSE)
  z <- stats::qnorm(0.95)
  expect_equal(ci$lower, unname(est - z * sqrt(diag(V))), tolerance = 1e-12)
  expect_equal(ci$upper, unname(est + z * sqrt(diag(V))), tolerance = 1e-12)
  expect_identical(confint(fit, level = 0.9, interval = NULL), ci)
})

test_that("two entries that match make the choice ambiguous", {
  fit <- pois_fit()
  set.seed(3)
  fit <- add_bootstrap(fit, n = 20, resample = case_resample())
  fit <- add_bootstrap(fit, n = 20, resample = parametric_resample())
  expect_error(confint(fit, interval = percentile_interval()),
               "More than one stored bootstrap matches")
  expect_no_error(suppressWarnings(confint(
    fit, level = 0.5,
    interval = percentile_interval(resample = case_resample()))))
})

test_that("pred_bootstrap reuses the stored parametric replicates", {
  set.seed(7)
  d <- data.frame(x = runif(40))
  d$y <- 1 + d$x + rnorm(40, sd = 0.4)
  nd <- data.frame(x = c(0.2, 0.7))
  ivl <- function(fit, n) {
    predict(fit, what = "response", newdata = nd,
            interval = prediction_interval(
              predictive = pred_bootstrap(n = n)))
  }
  fit <- statmod(y ~ x, g1, d)
  set.seed(11)
  fit <- add_bootstrap(fit, n = 30)
  # the stored path never reaches the engine
  testthat::local_mocked_bindings(
    bootstrap_engine = function(...) stop("the engine was called"))
  p_reuse <- ivl(fit, 20)
  expect_error(ivl(statmod(y ~ x, g1, d), 20), "the engine was called")
  # more replicates than stored: the predictive bootstrap makes its own
  expect_error(ivl(fit, 40), "the engine was called")
  testthat::local_mocked_bindings(
    bootstrap_engine = eng0)
  # the replicates are those of the stream from the same seed: the first 20
  # replicates of a run of 30 are those of a run of 20
  fit2 <- statmod(y ~ x, g1, d)
  set.seed(11)
  p_own <- ivl(fit2, 20)
  expect_identical(p_reuse, p_own)
})

test_that("pred_bootstrap reuse is faster than the refits and keeps a model with random effects out", {
  set.seed(8)
  d <- data.frame(x = runif(40))
  d$y <- 1 + d$x + rnorm(40, sd = 0.4)
  nd <- data.frame(x = 0.5)
  ivl <- function(fit) {
    predict(fit, what = "response", newdata = nd,
            interval = prediction_interval(
              predictive = pred_bootstrap(n = 20)))
  }
  fit <- statmod(y ~ x, g1, d)
  set.seed(1)
  fit <- add_bootstrap(fit, n = 20)
  t_reuse <- system.time(ivl(fit))[["elapsed"]]
  t_own <- system.time(ivl(statmod(y ~ x, g1, d)))[["elapsed"]]
  expect_lt(t_reuse, t_own)
  # a model with a random-effect term is not reused
  gg <- data.frame(g = factor(rep(1:8, each = 5)), x = rnorm(40))
  gg$y <- 1 + gg$x + rnorm(8)[gg$g] + rnorm(40)
  fr <- statmod(y ~ x + random(~ 1 | g), g1, gg)
  set.seed(2)
  fr <- add_bootstrap(fr, n = 20)
  testthat::local_mocked_bindings(
    bootstrap_engine = function(...) stop("the engine was called"))
  expect_error(
    predict(fr, what = "response", newdata = data.frame(x = 0, g = "new"),
            random = "zero",
            interval = prediction_interval(
              predictive = pred_bootstrap(n = 20))),
    "the engine was called")
})

test_that("pred_bootstrap reuses a refit = \"full\" entry with its hyperparameters", {
  set.seed(9)
  d <- data.frame(x = runif(60))
  d$y <- sin(2 * pi * d$x) + rnorm(60, sd = 0.3)
  nd <- data.frame(x = c(0.25, 0.5))
  ivl <- function(fit) {
    predict(fit, what = "response", newdata = nd,
            interval = prediction_interval(
              predictive = pred_bootstrap(n = 6, refit = "full")))
  }
  fit <- quiet(statmod(y ~ s(x), g1, d))
  set.seed(12)
  fit <- quiet(add_bootstrap(fit, n = 6, refit = "full"))
  testthat::local_mocked_bindings(
    bootstrap_engine = function(...) stop("the engine was called"))
  p_reuse <- ivl(fit)
  testthat::local_mocked_bindings(
    bootstrap_engine = eng0)
  fit2 <- quiet(statmod(y ~ s(x), g1, d))
  set.seed(12)
  p_own <- quiet(ivl(fit2))
  expect_equal(p_reuse, p_own, tolerance = 1e-10)
})
