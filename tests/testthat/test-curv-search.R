# The outer search reads a scaled SCAD or MCP at the curvature each
# alternation settles at, as the final refit does: the criterion recorded at
# the search's last point is the one the fit reports at the same
# hyperparameters. Reading the curvature written from the starting
# coefficients gave -879.101 against -878.793 on the SCAD below.

curv_search_data <- function() {
  set.seed(7)
  n <- 1000
  d <- data.frame(x = runif(n), z = runif(n), w = runif(n))
  d$y <- 2 * sin(6 * d$x) + 1.2 * cos(4 * d$z) + rnorm(n, 0, 0.5)
  d
}

test_that("the search criterion of a scaled SCAD matches the final one", {
  d <- curv_search_data()
  fit <- statmod(y ~ s(x) + scad(~ z + w, lambda = 6, standardize = TRUE),
                 distrib = distributions7::gaussian1_distrib(), data = d)
  h <- fit@history$outer
  last <- h$criterion[nrow(h)]
  expect_lt(abs(last - fit@criterion), 1e-5 * max(1, abs(fit@criterion)))
})

test_that("the search criterion of a scaled MCP matches the final one", {
  d <- curv_search_data()
  fit <- statmod(y ~ s(x) + mcp(~ z + w, lambda = 6, standardize = TRUE),
                 distrib = distributions7::gaussian1_distrib(), data = d)
  h <- fit@history$outer
  last <- h$criterion[nrow(h)]
  expect_lt(abs(last - fit@criterion), 1e-5 * max(1, abs(fit@criterion)))
})
