skip_on_cran()

test_that("one coordinate with no exact outer gradient is searched by brent()", {
  expect_identical(outer_default_optimizer(FALSE, FALSE, FALSE, 1L)@name, "brent")
  expect_identical(outer_default_optimizer(FALSE, FALSE, FALSE, 2L)@name,
                   "nelder-mead")
  expect_identical(outer_default_optimizer(TRUE, FALSE, FALSE, 1L)@name,
                   optimizers7::lbfgs()@name)
})

test_that("a marginal break-point is certified on the differenced criterion", {
  set.seed(4)
  m <- 10; ni <- 10
  id <- factor(rep(1:m, each = ni))
  t <- rep(seq(0, 10, length.out = ni), m)
  psi <- 5 + rnorm(m, sd = 0.7)
  y <- rnorm(m, sd = 0.5)[id] + 1.5 * (t > psi[id]) + rnorm(m * ni, sd = 0.4)
  d <- data.frame(id = id, t = t, y = y)
  f <- statmod(y ~ random(~1 | id) + jump(t, psi ~ random(~1 | id), marginal = TRUE),
               distrib = distributions7::gaussian1_distrib(), data = d)
  expect_identical(f@methods$search@name, "brent")
  ce <- statmod_certificate(f)
  expect_identical(ce$curvature, "differenced criterion")
  expect_identical(ce$state, "converged")
  expect_lt(abs(ce$gradient), 0.05)
  expect_lt(ce$decrement, 1e-2)
  # the random intercept's block is sparse beside the structural term, and
  # the information of the mixture assembles it into a dense matrix
  expect_no_error(s <- summary(f))
  expect_true(is.finite(confint(f)["mu:jump.delta1", "se"]))
})
