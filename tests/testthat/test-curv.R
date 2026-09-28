# SCAD and MCP are scaled by the curvature of the likelihood in each of their
# coordinates, read once at the starting coefficients (decision 4 of
# 2026-09-27): the penalty is sum_j c_j rho(b_j; lambda / c_j), so the knee
# sits at a lambda / c_j and the shape parameter means the same thing whatever
# the data's scale.

set.seed(11)
xc <- matrix(stats::rnorm(120 * 6, sd = rep(c(0.5, 1, 4), each = 240)),
             120, 6)
dc <- data.frame(y = 2 + xc[, 1] - xc[, 4] + stats::rnorm(120, sd = 1.5))
dc$x <- xc

curv_of_fit <- function(fit, fam) {
  tn <- grep(fam, names(fit@spec@terms$mu), fixed = TRUE, value = TRUE)
  expect_length(tn, 1L)
  modelterms7::term_penalty(fit@spec@terms$mu[[tn]])
}

test_that("a gaussian mean's curvature is sum x^2 over the starting variance", {
  # the reference is the information of a gaussian mean written out, read at
  # the intercept-only maximum the fit starts from; it shares no arithmetic
  # with statmod_information_at()
  s0 <- sqrt(mean((dc$y - mean(dc$y))^2))
  ref <- colSums(xc^2) / s0^2
  for (f in list(y ~ scad(x, lambda = 3, a = 3.7),
                 y ~ mcp(x, lambda = 3, gamma = 3))) {
    fit <- statmod(f, distributions7::gaussian1_distrib(), dc)
    pen <- curv_of_fit(fit, if (grepl("scad", deparse(f))) "scad" else "mcp")
    expect_equal(pen@curv, ref, tolerance = 1e-8)
  }
})

test_that("under standardize the curvature is on the standardized scale", {
  fit <- statmod(y ~ scad(x, lambda = 3, a = 3.7, standardize = TRUE),
                 distributions7::gaussian1_distrib(), dc)
  pen <- curv_of_fit(fit, "scad")
  d <- as.numeric(Matrix::diag(pen@map))
  s0 <- sqrt(mean((dc$y - mean(dc$y))^2))
  expect_equal(pen@curv, colSums(xc^2) / s0^2 / d^2, tolerance = 1e-8)
})

test_that("a lasso and a ridge are not touched", {
  for (f in list(y ~ lasso(x, lambda = 3), y ~ ridge(x, lambda = 3))) {
    fit <- statmod(f, distributions7::gaussian1_distrib(), dc)
    fam <- if (grepl("lasso", deparse(f))) "lasso" else "ridge"
    pen <- curv_of_fit(fit, fam)
    expect_false("curv" %in% S7::prop_names(pen))
  }
})

test_that("under standardize the fit does not depend on the units", {
  # a covariate rescaled by k has its coefficient divided by k, its spread
  # multiplied by k and its curvature by k^2, so on the standardized scale
  # both the slope and the knee are unchanged and so is the fit. Without
  # standardize the lasso part lambda |b| is not unit-free, as for a lasso.
  d2 <- dc
  d2$x <- xc * 10
  f <- y ~ scad(x, lambda = 3, a = 3.7, standardize = TRUE)
  a <- statmod(f, distributions7::gaussian1_distrib(), dc)
  b <- statmod(f, distributions7::gaussian1_distrib(), d2)
  expect_equal(fitted(a), fitted(b), tolerance = 1e-6)
  expect_equal(a@coefficients$mu[-1L], 10 * b@coefficients$mu[-1L],
               tolerance = 1e-6)
})
