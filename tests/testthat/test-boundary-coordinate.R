## A coordinate at a boundary is one coordinate, not the whole point.
##
## Where a parameter reaches the clamp its link keeps it strictly inside, the
## family's curvature there is NaN. Until 0.82.0 that single entry denied the
## whole system a solve and the whole criterion a Cholesky factor, so a fit
## stopped with every OTHER coordinate far from stationary and `statmod()`
## raised at its first evaluation. The coordinate is held instead.

test_that("pin_boundary leaves a finite matrix exactly alone", {
  K <- crossprod(matrix(c(2, 1, 0, 1, 3, 1, 0, 1, 4), 3))
  P <- pin_boundary(K)
  expect_equal(attr(P, "held"), 0L)
  expect_identical(as.numeric(P), as.numeric(K))
})

test_that("pin_boundary pins the frozen coordinate and NOT its neighbours", {
  ## the shape a boundary leaves: the coordinate's whole ROW is NaN, cross
  ## terms included. Testing columns rather than the diagonal marks the
  ## neighbours too, which is what held sigma along with nu and left the fit
  ## exactly where it had been.
  K <- diag(c(4, 5, 6))
  K[1, 2] <- K[2, 1] <- 1
  K[3, ] <- NaN
  K[, 3] <- NaN
  P <- pin_boundary(K)
  expect_equal(attr(P, "held"), 1L)
  expect_equal(P[3, 3], 1)
  expect_true(all(P[3, -3] == 0) && all(P[-3, 3] == 0))
  ## the block that was finite is untouched, cross term included
  expect_equal(P[1:2, 1:2], K[1:2, 1:2])
  expect_true(all(is.finite(P)))
})

test_that("iwls_solve holds a frozen coordinate and solves the rest", {
  A <- diag(c(4, 5, 6))
  u <- c(1, 2, 3)
  full <- iwls_solve(list(R = NULL, C = NULL, A = A), u, "chol")
  expect_equal(full$held, integer(0))
  expect_equal(full$delta, u / c(4, 5, 6))

  ## freeze the third: its row and column go, the first two are solved as
  ## though it were not there
  A2 <- A
  A2[3, ] <- NaN
  A2[, 3] <- NaN
  got <- iwls_solve(list(R = NULL, C = NULL, A = A2), u, "chol")
  expect_equal(got$held, 3L)
  expect_equal(got$delta[3], 0)
  expect_equal(got$delta[1:2], u[1:2] / c(4, 5))
  expect_equal(got$rank, 2L)

  ## and with every coordinate frozen there is nothing to solve
  none <- iwls_solve(list(R = NULL, C = NULL, A = matrix(NaN, 3, 3)), u, "chol")
  expect_equal(none$rank, 0L)
  expect_equal(none$delta, numeric(3))
})

test_that("a shape run to its clamp fits instead of stopping the run", {
  skip_on_cran()
  ## ⚠️ THE CLAMP AND NOT MERELY A LARGE VALUE. Two regimes have to be kept
  ## apart, and only the first is what this repair covers: at the clamp the
  ## curvature is NaN and the solve dies for every coordinate, which is
  ## repaired here; approaching it without reaching it the curvature is
  ## finite but negligible (3.5/nu^4 is 5e-52 at nu = 9e12), the step in that
  ## coordinate is astronomically long and the line search shrinks the WHOLE
  ## step to nothing, which is NOT repaired. The construction below lands in
  ## the first: measured, nu reaches double.xmax.
  set.seed(7)
  n <- 1000
  g <- factor(rep_len(seq_len(30), n))
  x <- runif(n)
  z <- runif(n)
  b <- rnorm(30, 0, 0.5)
  eta <- 2 * sin(6 * x) + 1.2 * cos(4 * z) + b[as.integer(g)]
  d <- data.frame(y = eta + rt(n, df = 5) * 0.5, x = x, z = z, g = g)

  fit <- statmod(y ~ s(x, k = 20) + random(~1 | g),
                 distributions7::student_t1_distrib(), d)
  expect_s7_class(fit, StatmodFit)
  expect_true(is.finite(as.numeric(logLik(fit))))
  ## nu is at the clamp, and the OTHER coordinates were fitted there rather
  ## than left where the failed solve had put them
  expect_gt(fit@fitted$nu[1], 1e300)
  expect_true(is.finite(fit@fitted$sigma[1]))
  expect_gt(fit@fitted$sigma[1], 0)
  ## A t at an unbounded nu is the gaussian, and the two PREDICTORS agree.
  ## Not the log-likelihoods: the held coordinate takes one off the Laplace
  ## dimension, which moves the criterion by log(2*pi)/2 and sends the two
  ## searches to different hyperparameters -- measured, -1491.8 against
  ## -1494.0. The predictor is what the two models share.
  gfit <- statmod(y ~ s(x, k = 20) + random(~1 | g),
                  distributions7::gaussian1_distrib(), d)
  expect_gt(stats::cor(fit@fitted$mu, gfit@fitted$mu), 0.999)
})
