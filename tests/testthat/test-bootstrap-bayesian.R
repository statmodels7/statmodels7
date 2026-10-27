g1b <- distributions7::gaussian1_distrib()
test_that("the Bayesian weights are positive, sum to n and keep zero weights", {
  set.seed(1)
  d <- data.frame(x = runif(40))
  d$y <- 1 + d$x + rnorm(40, sd = 0.4)
  w0 <- c(rep(1, 35), rep(0, 5))
  f <- statmod(y ~ x, g1b, d, weights = w0)
  spec0 <- f@spec
  design0 <- statmod_design(spec0)
  draw <- statmodels7:::bootstrap_plan(
    f, bayesian_resample(), spec0, design0,
    statmodels7:::random_modes(spec0, "zero"), TRUE)
  expect_identical(attr(draw, "kind"), "weights")
  set.seed(2)
  w <- draw()$w
  expect_equal(length(w), 40L)
  expect_true(all(w[1:35] > 0))
  expect_true(all(w[36:40] == 0))
  expect_equal(sum(w), 35, tolerance = 1e-12)
  # a different draw each time, and the stream determines the draw
  set.seed(2)
  expect_identical(draw()$w, w)
  expect_false(identical(draw()$w, w))
})

test_that("the Bayesian bootstrap has the variance of the case bootstrap", {
  d <- data.frame(x = {set.seed(3); runif(100)})
  d$y <- 1 + 2 * d$x + rnorm(100, sd = exp(2 * d$x) / 4)
  f <- statmod(y ~ x, g1b, d)
  set.seed(4)
  f <- add_bootstrap(f, n = 300, resample = bayesian_resample())
  f <- add_bootstrap(f, n = 300, resample = case_resample())
  sb <- sqrt(diag(vcov(f, vcov_type = bootstrap_vcov(
    resample = bayesian_resample()))))
  sc <- sqrt(diag(vcov(f, vcov_type = bootstrap_vcov(
    resample = case_resample()))))
  # HC0 written out for the mean's coefficients
  X <- cbind(1, d$x)
  r <- as.numeric(d$y - X %*% f@coefficients$mu)
  A <- solve(crossprod(X))
  hc0 <- sqrt(diag(A %*% crossprod(X * r) %*% A))
  tol <- 4 / sqrt(2 * 300)
  expect_lt(abs(sb[["mu:x"]] / hc0[[2]] - 1), tol)
  expect_lt(abs(sb[["mu:(Intercept)"]] / hc0[[1]] - 1), tol)
  expect_lt(abs(sb[["mu:x"]] / sc[["mu:x"]] - 1), 2 * tol)
})

test_that("the Bayesian scheme keeps a rare level that the case scheme loses", {
  set.seed(5)
  d <- data.frame(g = factor(c(rep("a", 97), rep("b", 3))), x = rnorm(100))
  d$y <- 1 + 0.5 * d$x + 0.4 * (d$g == "b") + rnorm(100, sd = 0.5)
  f <- statmod(y ~ x + g, g1b, d)
  set.seed(11)
  eb <- statmodels7:::bootstrap_engine(f, 150, bayesian_resample(),
                                       "coefficients", stream = FALSE)
  ec <- statmodels7:::bootstrap_engine(f, 150, case_resample(),
                                       "coefficients", stream = FALSE)
  expect_true(all(eb$status == 0L))
  # a draw that leaves the level out fails and is replaced
  expect_true(any(ec$status > 0L))
})

test_that("the Bayesian scheme is rejected where the case scheme is", {
  set.seed(6)
  gg <- data.frame(g = factor(rep(1:8, each = 5)), x = rnorm(40))
  gg$y <- 1 + gg$x + rnorm(8)[gg$g] + rnorm(40, sd = 0.5)
  fr <- statmod(y ~ x + random(~ 1 | g), g1b, gg)
  expect_error(add_bootstrap(fr, n = 3, resample = bayesian_resample()),
               "cluster_resample")
  expect_error(add_bootstrap(fr, n = 3, resample = bayesian_resample()),
               "bayesian_resample")
  expect_error(statmodels7:::as_resample_method("bayes"), "bayesian_resample")
  expect_output(print(bayesian_resample()), "bayesian")
  expect_true(S7::S7_inherits(bayesian_resample(), resample_method))
})

test_that("a Bayesian entry is stored and chosen like the others", {
  set.seed(7)
  d <- data.frame(x = runif(50))
  d$y <- 1 + d$x + rnorm(50, sd = 0.4)
  f <- statmod(y ~ x, g1b, d)
  set.seed(1)
  f <- add_bootstrap(f, n = 20, resample = bayesian_resample())
  e <- statmod_bootstrap(f, resample = bayesian_resample())
  expect_equal(e@n, 20L)
  expect_equal(e@resample@name, "bayesian")
  # 20 then 10 more are the 30 of one run from the same seed
  set.seed(1)
  f2 <- add_bootstrap(f, n = 30, resample = bayesian_resample())
  set.seed(1)
  f3 <- add_bootstrap(statmod(y ~ x, g1b, d), n = 30,
                      resample = bayesian_resample())
  expect_equal(unname(coef(statmod_bootstrap(f2, resample = bayesian_resample()))),
               unname(coef(statmod_bootstrap(f3, resample = bayesian_resample()))),
               tolerance = 1e-10)
})
