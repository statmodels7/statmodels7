# The joint reading of the outer curvature: a direction along which several
# coordinates move with the criterion unchanged, while each alone is curved.

ridge_case <- function(eps = 1e-9) {
  # four coordinates; the last three are an intercept and two contrasts, and
  # the criterion depends on them only through the second and third levels,
  # so (-1, +1, +1) on them is flat
  L <- rbind(c(1, 0, 0, 0),
             c(0, 1, 1, 0),
             c(0, 1, 0, 1),
             c(0.3, 0, 0, 0))
  A <- crossprod(L * sqrt(c(5, 1, 1, 1)))
  A[2:4, 2:4] <- A[2:4, 2:4] + eps
  A
}

test_that("a flat direction is found where every coordinate is curved", {
  A <- ridge_case()
  expect_true(all(diag(A) > 0.5))
  g <- c(1e-6, 1e-7, -1e-7, 2e-7)
  r <- certificate_ridges(g, A, seq_len(4), flat = 2e-3, tol = 1e-2)
  expect_length(r$dirs, 1L)
  x <- r$dirs[[1L]]
  x <- x * sign(x[3])
  expect_equal(x, c(0, -1, 1, 1), tolerance = 1e-6)
  expect_length(r$settled, 1L)
  # a gradient along the ridge is not settled: the direction stays under
  # test, and the verdict then cannot be read on it
  r2 <- certificate_ridges(c(0, -1, 1, 1) * 1e-2, A, seq_len(4), 2e-3, 1e-2)
  expect_length(r2$dirs, 1L)
  expect_length(r2$settled, 0L)
  # a curvature with no dependence has no flat direction
  expect_length(certificate_ridges(g, diag(4), seq_len(4), 2e-3, 1e-2)$dirs, 0L)
  # the sign of an eigenvalue this small is rounding: -1e-9 is named alike
  rn <- certificate_ridges(g, ridge_case(-1e-9), seq_len(4), 2e-3, 1e-2)
  expect_length(rn$dirs, 1L)
})

test_that("the decrement off a ridge is the constrained maximum", {
  A <- ridge_case()
  set.seed(3)
  g <- stats::rnorm(4) * 1e-3
  g <- g - sum(g * c(0, -1, 1, 1)) / 3 * c(0, -1, 1, 1)
  d0 <- decrement_off(g, A, seq_len(4), list(c(0, -1, 1, 1)))
  # by brute force: the maximum of 2 g'x - x'Ax over x with s^2 * r ' x = 0
  s2 <- diag(A)
  cc <- s2 * c(0, -1, 1, 1)
  N <- qr.Q(qr(cbind(cc)), complete = TRUE)[, -1]
  z <- solve(crossprod(N, A %*% N), crossprod(N, g))
  expect_equal(d0, 0.5 * sum(crossprod(N, g) * z), tolerance = 1e-10)
  # and it is finite where the full matrix is too flat to read reliably
  expect_true(is.finite(d0))
  # no larger than the maximum over every coordinate where that one is read
  B <- A + diag(4) * 1e-2
  expect_lte(decrement_off(g, B, seq_len(4), list(c(0, -1, 1, 1))),
             joint_decrement(g, B) + 1e-15)
})

test_that("the coefficients that move along a ridge lose their rows", {
  rows <- c("mu:(Intercept)", "mu:loc2", "k:(Intercept)", "k:loc2", "c:(Intercept)")
  # the inner mu coefficients move with k's combination by 0.4, c does not
  base <- diag(5) * 0.01
  r <- c(0.4, -0.4, 1, -1, 0)
  V <- base + 1e6 * tcrossprod(r)
  dimnames(V) <- list(rows, rows)
  ridge <- list(structure(c(`k/(Intercept)` = -1, `k/loc2` = 1),
                          where = c(3L, 4L)))
  expect_identical(ridge_coef_rows(ridge, V, rows), 1:4)
  # where the variance is missing on the ridge, the named coefficients alone
  V[4, ] <- V[, 4] <- NA_real_
  expect_identical(ridge_coef_rows(ridge, V, rows), 3:4)
  expect_identical(ridge_coef_rows(ridge, NULL, rows), 3:4)
  # a ridge over hyperparameters alone names no coefficient
  expect_identical(
    ridge_coef_rows(list(structure(c(a = 1, b = -1), where = c(NA, NA))),
                    V, rows), integer(0))
})
