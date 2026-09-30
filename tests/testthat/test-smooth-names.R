test_that("two smooths that would name their coefficients alike are rejected", {
  # up to statmodels7 0.170.0 this failed inside the fit on "duplicate
  # 'row.names' are not allowed", naming neither the terms nor the remedy
  set.seed(21)
  d <- data.frame(x = runif(100))
  d$y <- sin(2 * pi * d$x) + rnorm(100, sd = 0.3)
  expect_error(
    statmod(y ~ s(x, basis7::bspline_smooth(k = 10)) +
              s(x, basis7::pspline_smooth(k = 12)),
            distrib = gaussian1_distrib(), data = d),
    "same names.*label")
  # a label separates them
  fit <- statmod(y ~ s(x, basis7::bspline_smooth(k = 10)) +
                   s(x, basis7::pspline_smooth(k = 12), label = "p"),
                 distrib = gaussian1_distrib(), data = d,
                 outer_criterion = NULL)
  expect_false(anyDuplicated(names(coef(fit)$mu)) > 0L)
})

test_that("a smooth and its varying coefficient are fitted as in mgcv", {
  skip_if_not_installed("mgcv")
  set.seed(22)
  d <- data.frame(x = runif(300), z = rnorm(300))
  d$y <- sin(2 * pi * d$x) + (1 + 2 * d$x^2) * d$z + rnorm(300, sd = 0.5)
  fit <- statmod(y ~ s(x, basis7::bspline_smooth(k = 10)) +
                   s(x, basis7::bspline_smooth(k = 10), by = z),
                 distrib = gaussian1_distrib(), data = d)
  g <- mgcv::gam(y ~ s(x, bs = "bs", k = 10) + s(x, bs = "bs", k = 10, by = z),
                 data = d, method = "REML")
  expect_false(anyDuplicated(names(coef(fit)$mu)) > 0L)
  expect_lt(max(abs(fitted(fit) - stats::fitted(g))), 1e-4)
  expect_equal(sum(fit@edf$edf[fit@edf$parameter == "mu"]), sum(g$edf),
               tolerance = 1e-5)
})
