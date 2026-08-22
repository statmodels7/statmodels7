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

test_that("iwls_escalate starts at the curvature's own scale", {
  pc <- list(R = NULL, C = NULL, A = diag(c(2, 2328, 0.2357)))
  expect_equal(iwls_scale(pc), 2328)
  expect_equal(iwls_escalate(0, pc), 2328e-8)
  expect_equal(iwls_escalate(2328e-8, pc), 2328e-6)
  ## the augmented route reads the same diagonal off the pieces
  R <- matrix(c(1, 0, 0, 0, 3, 0, 0, 0, 0.5), 3)
  pc2 <- list(R = R, C = matrix(0, 1, 3), A = NULL)
  expect_equal(iwls_scale(pc2), 9)
})

test_that("the damping shortens the flat coordinate and not the others", {
  ## the disparity this exists for: one diagonal four orders below its
  ## neighbours, which a scalar step length cannot treat differently
  A <- diag(c(2328, 2328, 0.2357))
  u <- c(1, 1, 1)
  plain <- iwls_solve(list(R = NULL, C = NULL, A = A), u, "chol")$delta
  damped <- iwls_solve(list(R = NULL, C = NULL, A = A), u, "chol",
                       damp = 100)$delta
  ## the flat coordinate is shortened by orders, the curved ones barely move
  expect_lt(abs(damped[3]) / abs(plain[3]), 0.01)
  expect_gt(abs(damped[1]) / abs(plain[1]), 0.9)
  ## and a zero damping is the plain step exactly
  expect_identical(iwls_solve(list(R = NULL, C = NULL, A = A), u, "chol",
                              damp = 0)$delta, plain)
})

test_that("a shape run to its clamp fits instead of stopping the run", {
  skip_on_cran()
  ## ⚠️ TWO REGIMES, and this is the first: AT the clamp the curvature is
  ## NaN and the solve dies for every coordinate, which the hold repairs.
  ## The second -- approaching the clamp without reaching it, where the
  ## curvature is finite but negligible -- is the case below.
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

test_that("a shape APPROACHING its clamp no longer deadlocks the step", {
  skip_on_cran()
  ## The second regime. Here nu passes through about 9e12, where the
  ## information is 3.5/nu^4 = 5e-52: the curvature is finite so nothing is
  ## held, the scoring step in that coordinate is astronomically long, and
  ## the line search shrinks the WHOLE step to keep it admissible -- measured
  ## before the damping, 1, 0.125, 1.5e-05, 1.5e-08, and the run stalled at a
  ## score of 2.9e-02 with the mean nowhere near stationary.
  set.seed(11)
  n <- 300
  x <- runif(n)
  g <- factor(rep_len(seq_len(10), n))
  b <- rnorm(10, 0, 0.4)
  d <- data.frame(y = 2 * sin(4 * x) + b[as.integer(g)] + rnorm(n, 0, 0.5),
                  x = x, g = g)

  fit <- statmod(y ~ x + random(~1 | g),
                 distributions7::student_t1_distrib(), d,
                 outer_criterion = NULL)
  expect_true(fit@converged)
  h <- fit@history$inner
  ## the damping really did fire, and it started from zero
  expect_true(any(h$damp > 0))
  expect_identical(h$damp[1], 0)
  ## and it got PAST the stall: the loop's own score -- the penalized one it
  ## drives on, which is not the likelihood's -- had floored at 2.4e-02
  ## before the damping and reaches 1e-05 after it. The split is at the FIRST
  ## damped iteration and not at `damp == 0`, the damping decaying back to
  ## zero once it has done its work.
  first <- which(h$damp > 0)[1]
  expect_true(!is.na(first) && first > 1)
  expect_gt(min(h$score[seq_len(first - 1L)]), 1e-3)
  expect_lt(min(h$score[first:nrow(h)]), 1e-4)
})
