# A fit with modelterms7::scs(): the shape on a grid in each parameter, and
# the summary block of a smooth (the smoothing parameter, not the free
# coefficients).

scs_data <- function() {
  set.seed(31)
  x <- runif(200)
  f <- ifelse(x < 0.6, 1 - (1 - x / 0.6)^2, 1)
  data.frame(x = x, y = f + exp(-2.3 + 1.2 * f) * rnorm(200))
}

test_that("a fit with scs() in mu and sigma is monotone in both", {
  d <- scs_data()
  fit <- statmod(y ~ modelterms7::scs(x) | sigma ~ modelterms7::scs(x),
                 distrib = distributions7::gaussian1_distrib(), data = d)
  nd <- data.frame(x = seq(0, 1, length.out = 401))
  mu <- predict(fit, what = "mu", newdata = nd)
  sg <- predict(fit, what = "sigma", newdata = nd)
  expect_true(all(diff(mu) >= -1e-10))
  expect_true(all(diff(sg) >= -1e-10))
  expect_true(is.finite(as.numeric(logLik(fit))))
})

test_that("a convex scs() gives a convex fitted mean", {
  set.seed(32)
  d <- data.frame(x = runif(150))
  d$y <- (d$x - 0.4)^2 * 3 + rnorm(150, sd = 0.1)
  fit <- statmod(y ~ modelterms7::scs(x, modelterms7::shape_convex()),
                 distrib = distributions7::gaussian1_distrib(), data = d)
  mu <- predict(fit, what = "mu", newdata = data.frame(x = seq(0, 1, length.out = 401)))
  expect_true(all(diff(diff(mu)) >= -1e-10))
})

test_that("scs() is summarized as a smooth", {
  d <- scs_data()
  fit <- statmod(y ~ modelterms7::scs(x),
                 distrib = distributions7::gaussian1_distrib(), data = d)
  tm <- fit@spec@terms$mu[["modelterms7::scs(x)"]]
  if (is.null(tm)) tm <- Filter(function(t) S7::S7_inherits(t, modelterms7::ScsTerm),
                                fit@spec@terms$mu)[[1L]]
  expect_identical(term_block_kind(tm), "smooth")
  out <- paste(capture.output(print(summary(fit))), collapse = "\n")
  expect_match(out, "lambda [reml]", fixed = TRUE)
  expect_false(grepl("\\.g1", out))
})

test_that("a value at a point holds in a fit without an intercept and is inert with one", {
  set.seed(35)
  d <- data.frame(x = runif(150, 0.5, 3))
  d$y <- 2 * (1 - exp(-1.2 * d$x)) + rnorm(150, sd = 0.15)
  sm <- basis7::bspline_smooth(k = 8, lower = 0, upper = 3)
  g <- distributions7::gaussian1_distrib()
  f0 <- statmod(y ~ 0 + modelterms7::scs(
    x, modelterms7::shape_increasing_concave(value = c(at = 0, value = 0)), sm),
    distrib = g, data = d)
  expect_equal(predict(f0, what = "mu", newdata = data.frame(x = 0)), 0)
  mu <- predict(f0, what = "mu", newdata = data.frame(x = seq(0, 3, length.out = 301)))
  expect_true(all(diff(mu) >= -1e-10))
  expect_true(all(diff(diff(mu)) <= 1e-10))
  fits <- lapply(c(0, 5), function(v) statmod(
    y ~ modelterms7::scs(x, modelterms7::shape_increasing(value = c(at = 1, value = v)), sm),
    distrib = g, data = d))
  fc <- statmod(y ~ modelterms7::scs(x, modelterms7::shape_increasing(), sm),
                distrib = g, data = d)
  nd <- data.frame(x = c(0.7, 1.5, 2.8))
  for (f in fits) {
    expect_equal(as.numeric(logLik(f)), as.numeric(logLik(fc)), tolerance = 1e-8)
    expect_equal(predict(f, what = "mu", newdata = nd),
                 predict(fc, what = "mu", newdata = nd), tolerance = 1e-8)
  }
  expect_equal(coef(fits[[1L]])$mu[["(Intercept)"]] - coef(fits[[2L]])$mu[["(Intercept)"]],
               5, tolerance = 1e-8)
})

test_that("a slope at a point holds in a fitted convex mean", {
  set.seed(36)
  d <- data.frame(x = runif(150))
  d$y <- exp(2 * d$x) - 2 * d$x + rnorm(150, sd = 0.2)
  fit <- statmod(y ~ modelterms7::scs(
    x, modelterms7::shape_convex(slope = c(at = 0.5, value = 0.718))),
    distrib = distributions7::gaussian1_distrib(), data = d)
  h <- 1e-4
  mu <- function(x) predict(fit, what = "mu", newdata = data.frame(x = x))
  expect_equal((mu(0.5 + h) - mu(0.5 - h)) / (2 * h), 0.718, tolerance = 1e-6)
  grid <- seq(min(d$x), max(d$x), length.out = 401)
  expect_true(all(diff(diff(mu(grid))) >= -1e-10))
})
