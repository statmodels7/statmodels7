# The rank of a penalized augmented matrix is read off the JACOBI-
# EQUILIBRATED diagonal of its triangular factor, which is the correction
# solve_pd() took in 0.70.0 and which had not been propagated here. Since
# R'R = A'A, scaling A's columns by their norms scales that diagonal by the
# same factors, so |R_jj| / ||a_j|| is the diagonal of a decomposition with
# unit column norms. Without it a matrix was called deficient for having
# columns of different SIZE -- a large smoothing parameter makes the penalty
# rows of its own block enormous beside an unpenalized one -- and the whole
# solve fell through to a dense QR. Measured on `s(x) + random(~1|g)`, 87 of
# 127 solves were rejected and the dense route found FULL RANK in all 127.
test_that("the sparse route keeps a matrix whose columns differ in size", {
  skip_if_not_installed("Matrix")
  set.seed(31)
  n <- 2000; m <- 40
  g <- factor(sample(m, n, TRUE))
  Z <- Matrix::sparse.model.matrix(~ 0 + g)
  x <- runif(n)
  R <- Matrix::cbind2(Matrix::Matrix(cbind(1, x), sparse = TRUE), Z)
  p <- ncol(R)
  # a penalty on the indicator block alone, and a large one: the columns of
  # the augmented matrix then span many orders of magnitude in norm while
  # the matrix stays of full rank
  d <- c(0, 0, rep(sqrt(1e12), m))
  C <- Matrix::Diagonal(x = d)
  u <- rnorm(p)
  out <- statmodels7:::sparse_augmented_solve(R, C, u, "qr")
  expect_false(is.null(out))
  expect_identical(out$rank, p)

  # and it is the same increment the dense pivoted route gives
  A <- rbind(as.matrix(R), as.matrix(C))
  qa <- qr(A)
  expect_identical(qa$rank, p)
  ref <- numeric(p)
  ref[qa$pivot] <- backsolve(qr.R(qa),
                             forwardsolve(t(qr.R(qa)), u[qa$pivot]))
  expect_equal(out$delta, ref, tolerance = 1e-6)
})

test_that("an exactly collinear column is still refused", {
  skip_if_not_installed("Matrix")
  set.seed(32)
  n <- 500
  x <- runif(n)
  M <- cbind(1, x, x)                        # the third column repeats the second
  R <- Matrix::Matrix(M, sparse = TRUE)
  C <- Matrix::Diagonal(x = rep(0, 3))
  expect_null(statmodels7:::sparse_augmented_solve(R, C, rnorm(3), "qr"))
})

test_that("a nearly collinear column is aliased on the sparse route as on the dense one", {
  # the third column's part orthogonal to the others is about 2e-10 of its
  # norm: above p eps, which the sparse route used to test against, and
  # below the solver's 1e-7, which the dense pivoted QR tests against
  set.seed(4)
  x <- runif(50)
  M <- cbind(1, x, x + 1e-10 * rnorm(50))
  u <- as.numeric(crossprod(M, rnorm(50)))
  C <- matrix(0, nrow = 3, ncol = 3)
  dense <- statmodels7:::augmented_solve(M, C, u, "qr")
  sparse <- statmodels7:::augmented_solve(Matrix::Matrix(M, sparse = TRUE),
                                          Matrix::Diagonal(x = rep(0, 3)), u, "qr")
  expect_identical(dense$dropped, 3L)
  expect_identical(sparse$dropped, dense$dropped)
  expect_identical(sparse$rank, dense$rank)
  # well separated columns stay on the sparse route at full rank
  M2 <- cbind(1, x, x^2)
  ok <- statmodels7:::sparse_augmented_solve(Matrix::Matrix(M2, sparse = TRUE),
                                             Matrix::Diagonal(x = rep(0, 3)),
                                             u, "qr")
  expect_identical(ok$rank, 3L)
})
