# predict(interval = "group") and predict(interval = "prediction").

ri_data <- function(m = 15, ni = 8, seed = 2, sd_b = 0.7, pois = FALSE) {
  set.seed(seed)
  g <- factor(rep(seq_len(m), each = ni))
  x <- stats::runif(m * ni)
  b <- stats::rnorm(m, 0, sd_b)[as.integer(g)]
  y <- if (pois) stats::rpois(m * ni, exp(0.5 + x + b)) else
    1 + 2 * x + b + stats::rnorm(m * ni)
  data.frame(g = g, x = x, y = y)
}

test_that("a new group's predictor adds its effect's variance", {
  d <- ri_data()
  f <- statmod(y ~ x + random(~ 1 | g), distributions7::gaussian1_distrib(), d)
  nd <- data.frame(x = c(0.2, 0.8), g = factor("new"))
  c0 <- predict(f, "link:mu", nd, random = "zero", se = TRUE)
  g1 <- predict(f, "link:mu", nd, random = "zero", interval = "group")
  key <- grep("random", names(f@spec@terms$mu), value = TRUE)
  sb <- f@hyper$mu[[key]][["sigma"]]
  expect_equal(g1$fit, c0$fit)
  expect_equal(g1$se, sqrt(c0$se^2 + sb^2), tolerance = 1e-10)
  z <- stats::qnorm(0.975)
  expect_equal(g1$upper - g1$lower, 2 * z * g1$se, tolerance = 1e-10)
  # on the parameter scale the interval is the predictor's carried through
  # the link, here the identity
  gp <- predict(f, "mu", nd, random = "zero", interval = "group")
  expect_equal(gp$lower, g1$lower, tolerance = 1e-10)
})

test_that("a random effect on sigma widens sigma's interval for a new group", {
  set.seed(5)
  m <- 20
  g <- factor(rep(seq_len(m), each = 10))
  x <- stats::runif(m * 10)
  bs <- stats::rnorm(m, 0, 0.4)[as.integer(g)]
  d <- data.frame(g = g, x = x,
                  y = 1 + x + exp(-0.3 + bs) * stats::rnorm(m * 10))
  f <- statmod(y ~ x | sigma ~ random(~ 1 | g),
               distributions7::gaussian1_distrib(), d)
  nd <- data.frame(x = 0.5, g = factor("new"))
  c0 <- predict(f, "link:sigma", nd, random = "zero", se = TRUE)
  g1 <- predict(f, "link:sigma", nd, random = "zero", interval = "group")
  key <- grep("random", names(f@spec@terms$sigma), value = TRUE)
  sb <- f@hyper$sigma[[key]][["sigma"]]
  expect_equal(g1$se, sqrt(c0$se^2 + sb^2), tolerance = 1e-10)
  # mu carries no effect, so its interval is the confidence one
  cm <- predict(f, "mu", nd, random = "zero", se = TRUE)
  gm <- predict(f, "mu", nd, random = "zero", interval = "group")
  expect_equal(gm, cm)
})

test_that("the prediction interval is the quantiles of the averaged family", {
  skip_on_cran()
  d <- ri_data()
  f <- statmod(y ~ x + random(~ 1 | g), distributions7::gaussian1_distrib(), d)
  nd <- data.frame(x = 0.5, g = factor("new"))
  pr <- predict(f, "response", nd, random = "zero", interval = "prediction")
  # a reference by simulation, sharing only the covariance of the estimates:
  # the predictors drawn jointly Gaussian, the new effect from its prior, and
  # y from the family
  cm <- predict(f, "link:mu", nd, random = "zero", se = TRUE)
  cs <- predict(f, "link:sigma", nd, random = "zero", se = TRUE)
  key <- grep("random", names(f@spec@terms$mu), value = TRUE)
  sb <- f@hyper$mu[[key]][["sigma"]]
  set.seed(11)
  S <- 4e5
  em <- cm$fit + sqrt(cm$se^2 + sb^2) * stats::rnorm(S)
  es <- cs$fit + cs$se * stats::rnorm(S)
  y <- em + exp(es) * stats::rnorm(S)
  q <- stats::quantile(y, c(0.025, 0.5, 0.975), names = FALSE)
  expect_equal(c(pr$lower, pr$fit, pr$upper), q, tolerance = 5e-3)
  expect_equal(pr$se, stats::sd(y), tolerance = 5e-3)
})

test_that("a discrete family's prediction interval lies on its support", {
  skip_on_cran()
  d <- ri_data(pois = TRUE)
  f <- statmod(y ~ x + random(~ 1 | g), distributions7::poisson_distrib(), d)
  nd <- data.frame(x = c(0.1, 0.9), g = factor("new"))
  pr <- predict(f, "response", nd, random = "zero", interval = "prediction")
  expect_equal(pr$lower, round(pr$lower))
  expect_equal(pr$upper, round(pr$upper))
  expect_equal(pr$fit, round(pr$fit))
  cm <- predict(f, "link:mu", nd, random = "zero", interval = "group")
  set.seed(3)
  S <- 2e5
  for (i in 1:2) {
    lam <- exp(cm$fit[i] + cm$se[i] * stats::rnorm(S))
    y <- stats::rpois(S, lam)
    # the smallest value at which the averaged distribution function reaches
    # the level, which is what quantile(type = 1) returns
    q <- stats::quantile(y, c(0.025, 0.5, 0.975), names = FALSE, type = 1)
    expect_equal(c(pr$lower[i], pr$fit[i], pr$upper[i]), q)
  }
})

test_that("a family with no location parameter is predicted through its cdf", {
  set.seed(8)
  m <- 15
  g <- factor(rep(seq_len(m), each = 8))
  x <- stats::runif(m * 8)
  mu <- exp(0.3 + x + stats::rnorm(m, 0, 0.5)[as.integer(g)])
  d <- data.frame(g = g, x = x, y = stats::rgamma(m * 8, 3, 3 / mu))
  f <- statmod(y ~ x + random(~ 1 | g), distributions7::gamma1_distrib(), d)
  pr <- predict(f, "response", data.frame(x = 0.5, g = factor("new")),
                random = "zero", interval = "prediction")
  expect_true(all(is.finite(unlist(pr))))
  expect_true(pr$lower > 0 && pr$lower < pr$fit && pr$fit < pr$upper)
})

test_that("effects a label ties across equations enter the prediction jointly", {
  skip_on_cran()
  set.seed(9)
  m <- 20
  g <- factor(rep(seq_len(m), each = 10))
  x <- stats::runif(m * 10)
  b <- matrix(stats::rnorm(2 * m), m) %*% chol(matrix(c(0.5, 0.2, 0.2, 0.3), 2))
  d <- data.frame(g = g, x = x)
  d$y <- 1 + x + b[as.integer(g), 1] +
    exp(-0.2 + b[as.integer(g), 2]) * stats::rnorm(m * 10)
  f <- statmod(y ~ x + random(~ 1 | u | g) | sigma ~ random(~ 1 | u | g),
               distributions7::gaussian1_distrib(), d)
  nd <- data.frame(x = 0.5, g = factor("new"))
  pr <- predict(f, "response", nd, random = "zero", interval = "prediction")
  expect_true(all(is.finite(unlist(pr))))
  sp <- f@spec
  rm <- random_modes(sp, "zero")
  snd <- spec_at(f, nd, need_response = FALSE)
  design <- statmod_design(snd, split(rm$key, rm$param))
  bl <- random_blocks(snd, design, f, rm)
  expect_length(bl, 1L)
  expect_identical(nrow(bl[[1]]$members), 2L)
  expect_true(abs(bl[[1]]$Sigma[1, 2]) > 0)
  # a member left with its own estimate is refused
  expect_error(random_blocks(snd, design, f, rm[rm$param == "mu", ]),
               "read with its own")
})

test_that("each interval is asked for with what it applies to", {
  d <- ri_data()
  f <- statmod(y ~ x + random(~ 1 | g), distributions7::gaussian1_distrib(), d)
  nd <- data.frame(x = 0.5, g = factor("new"))
  expect_error(predict(f, "response", nd, random = "zero"), "prediction")
  expect_error(predict(f, "mu", nd, random = "zero", interval = "prediction"),
               "what = \"response\"")
  expect_error(predict(f, "mu", d, interval = "group"), "set aside")
  # conditional effects enter with their own estimates: a group the fit saw
  pr <- predict(f, "response", d[1:3, ], interval = "prediction")
  expect_true(all(pr$lower < pr$fit & pr$fit < pr$upper))
})
