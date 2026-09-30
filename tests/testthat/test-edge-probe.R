# A hyperparameter at the edge of its chart can look stationary on the free
# scale without being at a maximum. On nlme::Machines with nested random
# effects, the default search stopped with the Worker standard deviation at
# 0.00025 and a REML criterion of -110.63, where lme() reaches -107.84.

machines_df <- function() {
  data.frame(
    Worker = factor(rep(rep(1:6, each = 3), 3)),
    Machine = factor(rep(c("A", "B", "C"), each = 18)),
    score = c(52, 52.8, 53.1, 51.8, 52.8, 53.1, 60, 60.2, 58.4, 51.1, 52.3,
              50.3, 50.9, 51.8, 51.4, 46.4, 44.8, 49.2, 62.1, 62.6, 64, 59.7,
              60, 59, 68.6, 65.8, 69.7, 63.2, 62.8, 62.2, 64.8, 65, 65.4, 43.7,
              44.2, 43, 67.5, 67.2, 66.9, 61.5, 61.7, 62.3, 70.8, 70.6, 71,
              64.1, 66.2, 64, 72.1, 72, 71.1, 62, 61.4, 60.5))
}

G <- distributions7::gaussian1_distrib()

test_that("nested random effects reach the maximum of the REML criterion", {
  skip_on_cran()
  m <- machines_df()
  fit <- statmod(score ~ Machine + random(~ 1 | Worker) +
                   random(~ 1 | Worker:Machine), distrib = G, data = m)
  # lme(score ~ Machine, random = ~ 1 | Worker/Machine): 4.781050 and
  # 3.729532, restricted log-likelihood -107.8438
  expect_equal(hyper(fit)$estimate, c(4.781050, 3.729532), tolerance = 1e-4)
  expect_equal(as.numeric(logLik(fit, type = "marginal")), -107.8438,
               tolerance = 1e-6)
  expect_identical(statmod_certificate(fit)$state, "converged")
})

test_that("without the probe the search stops at the edge", {
  skip_on_cran()
  # the negative control: the case above reaches the maximum through the
  # probe and not by its own path
  m <- machines_df()
  testthat::local_mocked_bindings(edge_probe_hyper = function(...) NULL)
  fit <- statmod(score ~ Machine + random(~ 1 | Worker) +
                   random(~ 1 | Worker:Machine), distrib = G, data = m)
  expect_lt(as.numeric(logLik(fit, type = "marginal")), -108)
  expect_lt(hyper(fit)$estimate[1], 0.01)
})

test_that("the probe returns the best point inside and nothing at an interior point", {
  # a criterion rising from the edge to 0, on the scale the search minimizes
  fn <- function(v) -v[[1]] + v[[1]]^2 / 20
  p <- edge_probe_hyper(c(-9, 0.3), c(0, 0.3), 1L, fn(c(-9, 0.3)), fn)
  expect_identical(p, c(0, 0.3))
  # a coordinate past the edge where the criterion falls towards the inside
  # is left where it is
  fn2 <- function(v) v[[1]]
  expect_null(edge_probe_hyper(c(-9, 0.3), c(0, 0.3), 1L, fn2(c(-9, 0.3)),
                               fn2))
  # no coordinate past the edge: no evaluation at all
  calls <- 0L
  fn3 <- function(v) { calls <<- calls + 1L; 0 }
  expect_null(edge_probe_hyper(c(-2, 9), c(0, 0), 1L, 0, fn3))
  expect_identical(calls, 0L)
})
