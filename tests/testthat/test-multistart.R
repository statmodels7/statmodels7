# A structural term that asks for several starts through
# modelterms7::term_starts() is fitted from each, and the best is kept. On
# MASS::geyser with three regimes the start at the quantiles reaches the
# maximum, -1053.39, and a start placed at a local maximum stays there,
# -1210.49; starting from the second and offering the first must end at the
# first.

geyser3 <- function() data.frame(waiting = geyser_waiting)
G <- distributions7::gaussian1_distrib()

# the quantile start, taken before any test replaces the method
quant_start <- S7::method(modelterms7::term_start, modelterms7::RegimeTerm)

# a start at the local maximum, with two transition probabilities near their
# edge. The displaced start this replaced (the third of
# term_starts(regime(k = 3, n_start = 3)), seed 102) reached the local
# maximum while the regime was alternated with the coefficients; fitted in
# one system, it reaches the maximum, -1053.39.
bad_start <- function(term, target) {
  z <- quant_start(term, target = target)
  z[c("gap2", "gap3", "alr1.1", "alr1.2", "alr2.1", "alr2.2", "alr3.1",
      "alr3.2")] <- c(3.718, 1.816, -2.451, 3.584, -5.527, 8, -5.048, 8)
  z
}

with_regime_methods <- function(start, starts, expr) {
  RT <- modelterms7::RegimeTerm
  g1 <- modelterms7::term_start
  g2 <- modelterms7::term_starts
  k1 <- S7::method(g1, RT)
  k2 <- S7::method(g2, RT)
  if (!is.null(start)) S7::method(g1, RT) <- start
  if (!is.null(starts)) S7::method(g2, RT) <- starts
  on.exit({
    S7::method(g1, RT) <- k1
    S7::method(g2, RT) <- k2
  }, add = TRUE)
  suppressMessages(suppressWarnings(expr))
}

test_that("a fit from the displaced start alone stops at the local maximum", {
  skip_on_cran()
  fit <- with_regime_methods(
    start = function(term, ..., target = NULL) bad_start(term, target),
    starts = NULL,
    statmod(waiting ~ regime(k = 3), G, geyser3()))
  expect_lt(as.numeric(logLik(fit)), -1200)
})

test_that("a second start reaches the maximum and is kept", {
  skip_on_cran()
  fit <- with_regime_methods(
    start = function(term, ..., target = NULL) bad_start(term, target),
    starts = function(term, ..., target = NULL) {
      list(bad_start(term, target), quant_start(term, target = target))
    },
    statmod(waiting ~ regime(k = 3), G, geyser3()))
  expect_equal(as.numeric(logLik(fit)), -1053.392, tolerance = 1e-5)
})

test_that("a worse second start leaves the first fit in place", {
  skip_on_cran()
  one <- statmod(waiting ~ regime(k = 3), G, geyser3())
  two <- with_regime_methods(
    start = NULL,
    starts = function(term, ..., target = NULL) {
      list(quant_start(term, target = target), bad_start(term, target))
    },
    statmod(waiting ~ regime(k = 3), G, geyser3()))
  expect_identical(logLik(two), logLik(one))
  expect_identical(two@coefficients, one@coefficients)
  expect_identical(two@structural, one@structural)
})

test_that("the starts do not move the caller's random number stream", {
  skip_on_cran()
  set.seed(42)
  fit <- statmod(waiting ~ regime(k = 2, n_start = 2), G, geyser3())
  u1 <- runif(1)
  set.seed(42)
  expect_identical(u1, runif(1))
  expect_equal(as.numeric(logLik(fit)), -1099.633, tolerance = 1e-5)
})
