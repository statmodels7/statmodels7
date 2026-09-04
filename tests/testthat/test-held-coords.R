## Holding a coordinate the penalized information carries nothing about.
##
## A combination c'beta has a finite variance exactly where c is orthogonal
## to the null space of K, so a coefficient the null space does not touch
## keeps its variance whatever happens to the others.  What decides whether
## dropping a coordinate leaves the others alone is the Schur correction it
## removes, K_Aj K_jj^-1 K_jA: where that is negligible, holding IS the
## Moore-Penrose inverse restricted to what is kept, and the tests below pin
## the two against each other rather than trusting the derivation.
##
## Every case carries its opposite as well: a matrix that is invertible must
## hold nothing, and a flat direction spread over two coordinates must hold
## nothing either, since dropping one of them would report the variance
## CONDITIONAL on the other as if it were the marginal one.

pd_matrix <- function(p = 5, seed = 1) {
  set.seed(seed)
  X <- matrix(stats::rnorm(200 * p), 200, p)
  crossprod(X)
}

test_that("uninformative_coords finds an empty row and nothing else", {
  A <- pd_matrix()

  ## an invertible matrix holds nothing
  expect_identical(uninformative_coords(A), integer(0))

  ## a row and column that vanish
  B <- A
  B[3, ] <- 0
  B[, 3] <- 0
  expect_identical(uninformative_coords(B), 3L)

  ## two of them
  C <- A
  C[c(2, 4), ] <- 0
  C[, c(2, 4)] <- 0
  expect_identical(uninformative_coords(C), c(2L, 4L))

  ## a non-finite diagonal, which is what a parameter run out of its range
  ## leaves. Read BY ROW this marks every neighbour, the cross terms being
  ## non-finite too, so the rule is boundary_coords()'.
  E <- A
  E[2, ] <- NaN
  E[, 2] <- NaN
  expect_identical(uninformative_coords(E), 2L)
})

test_that("uninformative_coords holds nothing where holding would not be exact", {
  ## TWO COLLINEAR COLUMNS. The null vector is (e_2 - e_4)/sqrt(2) and each
  ## coordinate carries half of it, so neither can be dropped on its own:
  ## what is estimable there is the sum, not either coefficient.
  set.seed(2)
  X <- matrix(stats::rnorm(200 * 5), 200, 5)
  X[, 4] <- X[, 2]
  expect_identical(uninformative_coords(crossprod(X)), integer(0))

  ## A COORDINATE ON A TINY SCALE is not a flat one. Its row and its
  ## diagonal are both small, and the Schur correction sees through both:
  ## here it is of order one, so dropping the coordinate would move every
  ## other variance.
  A <- pd_matrix()
  G <- A
  G[3, ] <- G[3, ] * 1e-14
  G[, 3] <- G[, 3] * 1e-14
  expect_identical(uninformative_coords(G), integer(0))

  ## and neither is a coordinate on an enormous one
  H <- A
  H[3, 3] <- H[3, 3] * 1e14
  expect_identical(uninformative_coords(H), integer(0))

  ## A ROW THAT VANISHES WITH A POSITIVE DIAGONAL is estimable, at a
  ## standard error of 1e15. solve_pd() inverts that matrix, so the hold is
  ## never reached; it must not claim the coordinate either.
  K <- A
  K[3, ] <- 0
  K[, 3] <- 0
  K[3, 3] <- 1e-30
  expect_identical(uninformative_coords(K), integer(0))
})

test_that("holding is the pseudo-inverse where it fires", {
  ## The two share no arithmetic: one inverts a submatrix, the other
  ## reconstructs from the spectrum of the whole.
  A <- pd_matrix()
  for (drop in list(3L, c(2L, 4L))) {
    M <- A
    M[drop, ] <- 0
    M[, drop] <- 0
    j <- uninformative_coords(M)
    expect_identical(j, drop)
    k <- setdiff(seq_len(ncol(M)), j)
    e <- eigen(M, symmetric = TRUE)
    pos <- e$values > 1e-10 * max(abs(e$values))
    pinv <- e$vectors[, pos] %*% diag(1 / e$values[pos]) %*%
      t(e$vectors[, pos])
    expect_equal(solve(M[k, k, drop = FALSE]), pinv[k, k, drop = FALSE],
                 tolerance = 1e-12)
  }
})

test_that("a fit whose information is invertible is untouched", {
  ## The hold runs only where solve_pd() has already refused, so nothing
  ## that works today changes. Two shapes, one with a modelled scale.
  set.seed(11)
  n <- 300
  dd <- data.frame(x1 = stats::runif(n, -1, 1), x2 = stats::runif(n, -1, 1))
  dd$y <- stats::rnorm(n, 0.6 + 1.2 * dd$x1 - 0.7 * dd$x2, 0.5)

  f <- statmod(y ~ x1 + x2 | sigma ~ x1,
               distributions7::gaussian1_distrib(), dd)
  V <- expect_silent(vcov(f))
  expect_false(anyNA(diag(V)))

  dd$yc <- stats::rpois(n, exp(0.6 + 1.2 * dd$x1))
  g <- statmod(yc ~ x1 + x2, distributions7::poisson_distrib(), dd)
  expect_false(anyNA(diag(vcov(g))))
})

test_that("a design of deficient rank is still refused", {
  ## Two identical columns leave a flat direction that is a COMBINATION,
  ## and there is no coordinate to hold: the refusal stands rather than a
  ## conditional variance being reported as a marginal one.
  set.seed(3)
  n <- 200
  dd <- data.frame(x1 = stats::runif(n, -1, 1), x2 = stats::runif(n, -1, 1))
  dd$x3 <- dd$x1
  dd$y <- stats::rnorm(n, 1 + dd$x1 - dd$x2, 0.5)
  f <- statmod(y ~ x1 + x2 + x3, distributions7::gaussian1_distrib(), dd)
  expect_error(vcov(f), "not positive definite")
})
