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
