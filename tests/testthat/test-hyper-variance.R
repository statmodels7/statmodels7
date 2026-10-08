## The variance of the hyperparameters a marginal criterion estimated.
##
## THE FULL INVERSE IS THE ANSWER WHENEVER IT IS ONE. Holding a coordinate is
## reached only where the negative of the outer Hessian cannot be inverted
## into a variance -- a diagonal entry that is zero, negative or not finite,
## which is what a hyperparameter at the edge of its range or a search
## stopped before a maximum leaves behind. Where that happens the coordinate
## is held and the rest inverted, which is the variance CONDITIONAL on it,
## and that equals the marginal one only where the coupling contributes
## nothing to the kept curvature -- so the Schur correction is computed and
## compared rather than assumed small.

test_that("a usable Hessian gives the FULL inverse and nothing else", {
  set.seed(1)
  for (p in 2:5) {
    B <- matrix(stats::rnorm(p * p), p)
    A <- crossprod(B) + diag(p)
    dimnames(A) <- list(letters[seq_len(p)], letters[seq_len(p)])
    ## identical(), not equal(): the full inverse must be RETURNED, not
    ## recomputed by another route that happens to agree
    expect_identical(hyper_variance(A), solve(A), info = paste("p =", p))
  }
})

test_that("a wide but honest variance is still the full inverse", {
  ## a hyperparameter the criterion barely pins down has an enormous
  ## variance and that is an answer, not a failure: nothing is held
  A <- matrix(c(1e-8, 0, 0, 2), 2, 2,
              dimnames = list(c("a", "b"), c("a", "b")))
  V <- hyper_variance(A)
  expect_identical(V, solve(A))
  expect_equal(V[["a", "a"]], 1e8)
})

test_that("a coordinate with the wrong curvature is held, the rest inverted", {
  ## the matrix measured on a Poisson-inverse gaussian with a random effect
  ## in each equation, stopped before a maximum: -H has eigenvalues 1.1017
  ## and -10.38, so the full inverse would report a NEGATIVE variance
  A <- matrix(c(1.101698, 0.024245, 0.024245, -10.379824), 2, 2,
              dimnames = list(c("mu", "alpha"), c("mu", "alpha")))
  expect_true(any(diag(solve(A)) < 0))          # the full route is unusable
  V <- hyper_variance(A)
  expect_false(is.null(V))
  expect_identical(dim(V), c(2L, 2L))
  expect_identical(dimnames(V), dimnames(A))
  expect_true(is.finite(V[["mu", "mu"]]))
  expect_true(is.na(V[["alpha", "alpha"]]))
  ## what the kept block gives is the inverse of its own curvature
  expect_equal(V[["mu", "mu"]], 1 / A[["mu", "mu"]])
})

test_that("the coupling is tested and not assumed negligible", {
  ## the same shape with the off-diagonal made large: holding the broken
  ## coordinate would move the kept one's curvature, so the conditional
  ## variance is NOT the marginal one and the whole matrix is refused
  A <- matrix(c(1, 0.5, 0.5, -0.1), 2, 2,
              dimnames = list(c("a", "b"), c("a", "b")))
  expect_null(hyper_variance(A))
  ## and the threshold is where the argument says it is: a correction of
  ## 1e-6 of the kept diagonal is taken, one of 1e-2 is not
  small <- A; small[1, 2] <- small[2, 1] <- sqrt(1e-6 * 0.1)
  expect_false(is.null(hyper_variance(small)))
  big <- A; big[1, 2] <- big[2, 1] <- sqrt(1e-2 * 0.1)
  expect_null(hyper_variance(big))
})

test_that("an indefinite matrix gives no variance, whatever its inverse's diagonal", {
  ## the shape measured on a zero-inflated negative binomial with a random
  ## effect on mu and on zi: every diagonal entry of -H negative and the
  ## inverse's diagonal positive, which the old test of the inverse accepted
  A <- matrix(c(-1, 3, 3, -8), 2, 2,
              dimnames = list(c("a", "b"), c("a", "b")))
  expect_true(all(diag(solve(A)) > 0))
  expect_null(hyper_variance(A))
  ## a positive diagonal does not make a matrix positive definite either
  ## a Householder reflection: diagonal 1/3, eigenvalues (1, 1, -1), its
  ## own inverse
  B <- diag(3) - 2 / 3 * matrix(1, 3, 3)
  dimnames(B) <- list(letters[1:3], letters[1:3])
  expect_lt(min(eigen(B, symmetric = TRUE)$values), 0)
  expect_true(all(diag(solve(B)) > 0))
  expect_null(hyper_variance(B))
  ## and the kept block is held to the same test: the third coordinate is
  ## held, and what is left is still indefinite
  C <- matrix(0, 3, 3, dimnames = list(letters[1:3], letters[1:3]))
  C[1:2, 1:2] <- matrix(c(1, 2, 2, 1), 2, 2)
  C[3, 3] <- -1
  expect_null(hyper_variance(C))
})

test_that("nothing usable gives NULL", {
  A <- diag(c(-1, -2))
  dimnames(A) <- list(c("a", "b"), c("a", "b"))
  expect_null(hyper_variance(A))
  B <- matrix(c(1, NA, NA, 2), 2, 2, dimnames = dimnames(A))
  expect_null(hyper_variance(B))
})

test_that("a fit that converges keeps the variance it had", {
  ## the negative control for the whole change: where the criterion reaches
  ## a maximum the full inverse is usable, so the hyperparameter's standard
  ## error and interval are what they were before a coordinate could be held
  skip_on_cran()
  set.seed(9)
  m <- 30
  ni <- 20
  n <- m * ni
  g <- factor(rep(seq_len(m), each = ni))
  u <- stats::rnorm(m, 0, 0.7)
  x <- stats::runif(n, -1, 1)
  dat <- data.frame(y = stats::rnorm(n, 1 + 0.8 * x + u[as.integer(g)], 1),
                    x = x, g = g)
  fit <- statmod(y ~ x + modelterms7::random(~1 | g), gaussian1_distrib(),
                 data = dat)
  sp <- fit@spec
  des <- statmod_design(sp)
  V <- statmod_hyper_vcov(sp, des, fit@coefficients, fit@hyper,
                          fit@methods$outer)
  expect_false(is.null(V))
  expect_true(all(is.finite(V)))
  tb <- summary(fit)@tables$mu[[2L]]$table
  r <- tb$role == "estimated"
  expect_true(any(r))
  expect_true(all(is.finite(tb$se[r])))
  expect_true(all(is.finite(tb$lower[r])) && all(is.finite(tb$upper[r])))
  expect_true(all(tb$lower[r] < tb$estimate[r]))
  expect_true(all(tb$upper[r] > tb$estimate[r]))
})

test_that("a coordinate at the edge is held without the Schur test", {
  ## the case measured on two correlated random intercepts: the angle's own
  ## curvature is 8e-7 and its coupling to the two standard deviations is
  ## above the criterion's resolution, so the whole matrix is indefinite
  A <- matrix(c(34.015, -16.941, 3.8598e-03,
                -16.941, 18.563, 5.1931e-03,
                3.8598e-03, 5.1931e-03, 7.9943e-07), 3, 3,
              dimnames = list(c("a", "b", "z"), c("a", "b", "z")))
  expect_null(hyper_variance(A))
  V <- hyper_variance(A, at_edge = 3L)
  expect_false(is.null(V))
  expect_true(all(is.na(V[3, ])) && all(is.na(V[, 3])))
  expect_equal(V[1:2, 1:2], solve(A[1:2, 1:2]), tolerance = 1e-12)
  ## a held coordinate is held even where the full matrix would invert, so
  ## the others read the same conditional variance at either point
  B <- A
  B[3, 1:2] <- B[1:2, 3] <- c(1e-5, -1e-5)
  B[3, 3] <- 1e-3
  expect_equal(hyper_variance(B, at_edge = 3L)[1:2, 1:2],
               solve(B[1:2, 1:2]), tolerance = 1e-12)
  ## a coordinate held for its curvature alone keeps the Schur test
  C <- A
  C[3, 3] <- -1
  C[3, 1:2] <- C[1:2, 3] <- c(1, 1)
  expect_null(hyper_variance(C))
  expect_false(is.null(hyper_variance(C, at_edge = 3L)))
})
