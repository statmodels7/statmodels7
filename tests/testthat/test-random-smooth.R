rsm_data <- function(m = 12, ni = 20, seed = 2) {
  set.seed(seed)
  d <- data.frame(g = factor(rep(seq_len(m), each = ni)), x = runif(m * ni))
  a <- stats::rnorm(m, 0, 0.6)
  b1 <- stats::rnorm(m, 0, 0.5)
  w <- stats::rnorm(m, 0, 0.4)
  i <- as.integer(d$g)
  d$y <- sin(2 * pi * d$x) + a[i] + b1[i] * d$x +
    w[i] * sin(4 * pi * d$x) + stats::rnorm(nrow(d), 0, 0.25)
  d
}

test_that("a smooth random effect is the composition of its two parts", {
  d <- rsm_data()
  sm <- basis7::bspline_smooth(k = 6)
  f <- statmod(y ~ s(x) + random(~ s(x, smoother = sm) | g),
               distributions7::gaussian1_distrib(), data = d)
  # the same model written as a random level and slope on the smooth's own
  # linear column, and a per-group smooth without its null space
  d$lin <- as.matrix(modelterms7::term_matrix(
    modelterms7::term_build(modelterms7::s(x, smoother = sm), d)))[, 1L]
  g <- statmod(y ~ s(x) + random(~ 1 + lin | g) +
                 s(x, smoother = basis7::bspline_smooth(k = 6,
                                                        null_space = "drop"),
                   by = g),
               distributions7::gaussian1_distrib(), data = d)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(g)), tolerance = 1e-5)
  expect_equal(fitted(f), fitted(g), tolerance = 1e-3)
})

test_that("the summary names the null-space coordinates and the smooth's sd", {
  d <- rsm_data()
  f <- statmod(y ~ s(x) + random(~ s(x, smoother = basis7::bspline_smooth(
    k = 6)) | g), distributions7::gaussian1_distrib(), data = d)
  out <- paste(capture.output(print(summary(f))), collapse = "\n")
  expect_match(out, "sd[(Intercept)]", fixed = TRUE)
  expect_match(out, "sd[lin]", fixed = TRUE)
  expect_match(out, "cor[(Intercept), lin]", fixed = TRUE)
  expect_match(out, "smooth.sigma", fixed = TRUE)
})

test_that("prediction sets the random curves aside and integrates them", {
  d <- rsm_data()
  f <- statmod(y ~ s(x) + random(~ s(x, smoother = basis7::bspline_smooth(
    k = 6)) | g), distributions7::gaussian1_distrib(), data = d)
  nd <- d[1:4, ]
  n2 <- nd
  n2$g <- factor("7", levels = levels(d$g))
  n3 <- nd
  n3$g <- "new"
  z1 <- predict(f, what = "mu", newdata = nd, random = "zero")
  expect_equal(predict(f, what = "mu", newdata = n2, random = "zero"), z1)
  expect_equal(predict(f, what = "mu", newdata = n3, random = "zero"), z1)
  expect_false(isTRUE(all.equal(
    predict(f, what = "mu", newdata = nd), z1)))
  # under the identity link the average over new curves is the zero curve,
  # up to the Monte Carlo error of the draws; the variance adds theirs
  mg <- predict(f, what = "mu", newdata = n3, random = "marginal")
  expect_equal(mg, z1, tolerance = 0.05)
  v0 <- predict(f, what = "variance", newdata = n3, random = "zero")
  v1 <- predict(f, what = "variance", newdata = n3, random = "marginal")
  expect_true(all(v1 > v0))
})

test_that("a term with one prior per column is integrated as one Gaussian", {
  d <- rsm_data()
  f <- statmod(y ~ x + random(~ x | g, correlated = FALSE),
               distributions7::gaussian1_distrib(), data = d)
  n3 <- d[1:3, ]
  n3$g <- "new"
  v <- predict(f, what = "variance", newdata = n3, random = "marginal")
  hy <- f@hyper$mu
  s2 <- exp(2 * as.numeric(coef(f)$sigma))
  sd0 <- hy[[grep("::\\(Intercept\\)$", names(hy))]][["sigma"]]
  sd1 <- hy[[grep("::x$", names(hy))]][["sigma"]]
  expect_equal(v, s2 + sd0^2 + sd1^2 * n3$x^2, tolerance = 0.05,
               ignore_attr = TRUE)
})

test_that("an adaptive smooth random effect is fitted and simulated", {
  d <- rsm_data()
  fm <- y ~ s(x) + random(~ s(x, smoother = basis7::adaptive_smooth(
    k = 10, m = 3)) | g)
  f <- statmod(fm, distributions7::gaussian1_distrib(), data = d)
  out <- paste(capture.output(print(summary(f))), collapse = "\n")
  expect_match(out, "sigma3", fixed = TRUE)
  expect_true(is.finite(as.numeric(logLik(f))))
  set.seed(5)
  sim <- rstatmod(fm, distributions7::gaussian1_distrib(), data = d)
  expect_true(all(is.finite(sim$y)))
})
