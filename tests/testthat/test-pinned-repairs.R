# A linear column aliased with the null space of s(x): the criterion and the
# hyperparameter correction to the degrees of freedom are read with the
# coordinate pinned, and the correction equals that of the model without the
# aliased column.

pinned_case <- function() {
  set.seed(4)
  n <- 300
  d <- data.frame(x = runif(n))
  d$y <- rpois(n, exp(0.5 + sin(3 * d$x)))
  fit <- suppressWarnings(statmod(y ~ x + s(x),
                                  distrib = distributions7::poisson_distrib(),
                                  data = d))
  ref <- statmod(y ~ s(x), distrib = distributions7::poisson_distrib(),
                 data = d)
  list(fit = fit, ref = ref)
}

test_that("the criterion pins a deficient coordinate as its gradient does", {
  cs <- pinned_case()
  fit <- cs$fit
  sp <- fit@spec
  dz <- statmod_design(sp)
  H <- statmod_information_at(sp, fit@coefficients, dz, FALSE)
  S <- statmod_penalty_at(sp, fit@coefficients, fit@hyper, dz, "hessian")
  fl <- deficient_coords(as_dense(H + S))
  expect_length(fl, 1L)
  dzp <- dz
  attr(dzp, "pinned") <- fl
  c0 <- statmod_marginal(sp, dz, fit@coefficients, fit@hyper,
                         fit@methods$outer)
  c1 <- statmod_marginal(sp, dzp, fit@coefficients, fit@hyper,
                         fit@methods$outer)
  expect_false(is.null(c0))
  expect_equal(c0$value, c1$value, tolerance = 1e-12)
  expect_identical(c0$q, c1$q)
})

test_that("the edf correction pins a deficient coordinate", {
  cs <- pinned_case()
  fit <- cs$fit
  ref <- cs$ref
  for (ex in c(FALSE, TRUE)) {
    e1 <- statmod_edf_correction(fit@spec, fit@coefficients, fit@hyper,
                                 statmod_design(fit@spec), fit@methods$outer,
                                 expected = ex, pinned = fit@methods$pinned)
    e2 <- statmod_edf_correction(ref@spec, ref@coefficients, ref@hyper,
                                 statmod_design(ref@spec), ref@methods$outer,
                                 expected = ex, pinned = ref@methods$pinned)
    expect_gt(e1$total, 0)
    expect_equal(e1$total, e2$total, tolerance = 1e-8)
  }
})

test_that("pinned_inverse() zeroes the pinned rows of the inverse", {
  set.seed(1)
  X <- matrix(rnorm(40), 10, 4)
  K <- crossprod(cbind(X, X[, 4]))
  V <- pinned_inverse(K, integer(0))
  j <- which(rowSums(V != 0) == 0)
  expect_length(j, 1L)
  expect_true(j %in% 4:5)
  expect_equal(V[-j, -j], solve(K[-j, -j]), tolerance = 1e-10)
  V2 <- pinned_inverse(K, 2L)
  expect_true(all(V2[2, ] == 0))
  expect_null(pinned_inverse(matrix(0, 3, 3)))
})

test_that("a pinned coordinate does not decide the Cholesky test", {
  set.seed(2)
  A <- crossprod(matrix(rnorm(60), 12, 5)) * 1e9
  K <- matrix(0, 6, 6)
  K[1:5, 1:5] <- A
  K[6, 6] <- 7
  Kp <- pin_boundary(K, 6L)
  expect_identical(pinned_free(Kp), 1:5)
  r <- pd_factor(Kp)
  expect_true(r$ok)
  expect_false(is.null(r$factor))
  expect_equal(r$logdet, as.numeric(determinant(A)$modulus), tolerance = 1e-12)
  # a matrix changed after pinning keeps the test on the whole matrix
  Kc <- Kp
  Kc[6, 1] <- Kc[1, 6] <- 1e-3
  expect_null(pinned_free(Kc))
  # the sparse route
  Ks <- Matrix::bdiag(Matrix::Matrix(A, sparse = TRUE), Matrix::Diagonal(400))
  Ks <- methods::as(Ks, "CsparseMatrix")
  attr(Ks, "held_at") <- 6:405
  expect_identical(pinned_free(Ks), 1:5)
  rs <- pd_factor(Ks)
  expect_true(rs$ok)
  expect_false(is.null(rs$factor))
  expect_equal(rs$logdet, as.numeric(determinant(A)$modulus), tolerance = 1e-10)
})
