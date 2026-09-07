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

test_that("a design of deficient rank is aliased, not refused", {
  ## Two identical columns leave a flat direction that is a COMBINATION, so
  ## there is no coordinate to HOLD in the sense of the tests above -- each
  ## carries half of it. What settles which column goes is the pivot that
  ## fitted the model, exactly as in lm() and glm(): the column it left out
  ## has no estimate, and the model whose variance is reported is the one
  ## without it. Refusing the whole matrix instead lost every standard error
  ## for the sake of one that does not exist.
  set.seed(3)
  n <- 200
  dd <- data.frame(x1 = stats::runif(n, -1, 1), x2 = stats::runif(n, -1, 1))
  dd$x3 <- dd$x1
  dd$y <- stats::rnorm(n, 1 + dd$x1 - dd$x2, 0.5)
  f <- statmod(y ~ x1 + x2 + x3, distributions7::gaussian1_distrib(), dd)

  ## the same column base R aliases, and the same estimates
  g <- stats::glm(y ~ x1 + x2 + x3, data = dd)
  expect_identical(f@aliased, "mu:x3")
  expect_identical(names(coef(g))[is.na(coef(g))], "x3")
  expect_equal(unname(coef(f)$mu[1:3]), unname(coef(g)[1:3]), tolerance = 1e-6)

  V <- vcov(f)
  expect_true(is.na(V["mu:x3", "mu:x3"]))
  expect_false(anyNA(diag(V)[c("mu:(Intercept)", "mu:x1", "mu:x2")]))
  ## and the count of what was spent drops with it, or every criterion built
  ## on the total carries a parameter the fit never estimated
  expect_equal(sum(f@edf$edf), 4)

  ## AND EVERY ROUTE NAMES THE SAME COLUMN. Only iwls() on a pivoting
  ## decomposition reports a dropped column as a by-product of the solve;
  ## the others are answered after the fact by deficient_coords(), on the
  ## information at the mode. The two mechanisms share no arithmetic, so
  ## their agreeing on the coordinate AND on the variance of what is kept
  ## is what says the after-the-fact route reports the same model.
  se <- function(f) sqrt(diag(vcov(f)))[c("mu:(Intercept)", "mu:x1", "mu:x2")]
  ref <- se(f)
  for (m in list(optimizers7::newton(), optimizers7::bfgs(),
                 optimizers7::lbfgs(), iwls(decomposition = "chol"),
                 iwls(decomposition = "svd"))) {
    g <- statmod(y ~ x1 + x2 + x3, distributions7::gaussian1_distrib(), dd,
                 inner_optimizer = m)
    expect_identical(g@aliased, "mu:x3", info = class(m)[1])
    expect_true(is.na(vcov(g)["mu:x3", "mu:x3"]), info = class(m)[1])
    expect_equal(se(g), ref, tolerance = 1e-5, info = class(m)[1])
    expect_equal(sum(g@edf$edf), 4, tolerance = 1e-6, info = class(m)[1])
  }
})


test_that("deficient_coords names what a pivot would, and nothing else", {
  ## The after-the-fact rank test, at unit level. Every case carries its
  ## opposite: what must be named, and what must not be, since the whole
  ## risk of the thing is naming a coordinate the data does identify.
  set.seed(1)
  X <- matrix(stats::rnorm(200 * 5), 200, 5)
  A <- crossprod(X)
  dup <- function(j, k) { Y <- X; Y[, j] <- Y[, k]; crossprod(Y) }

  ## nothing to name
  expect_identical(deficient_coords(A), integer(0))
  expect_identical(deficient_coords(matrix(3)), integer(0))

  ## an exact duplication: the LATER column goes, as the pivot leaves it
  expect_identical(deficient_coords(dup(4, 2)), 4L)
  Y <- X; Y[, 4] <- Y[, 2]; Y[, 5] <- Y[, 2]
  expect_identical(deficient_coords(crossprod(Y)), c(4L, 5L))

  ## BY THE DIAGONAL: a coordinate at a link's clamp makes its whole row
  ## non-finite, so a row test would name its neighbours; and neither a
  ## clamp nor an empty row is an ALIASING -- there the estimate stands and
  ## only the variance does not, which is uninformative_coords()' business
  C <- A; C[2, ] <- NaN; C[, 2] <- NaN
  expect_identical(deficient_coords(C), integer(0))
  D <- A; D[3, ] <- 0; D[, 3] <- 0
  expect_identical(deficient_coords(D), integer(0))
  ## and a clamp beside a real duplication must not hide it
  E <- dup(4, 2); E[5, ] <- NaN; E[, 5] <- NaN
  expect_identical(deficient_coords(E), 4L)

  ## SCALE IS NOT DEFICIENCY, which is what the equilibration buys
  G <- A; G[3, ] <- G[3, ] * 1e-14; G[, 3] <- G[, 3] * 1e-14
  expect_identical(deficient_coords(G), integer(0))
  H <- A; H[3, 3] <- H[3, 3] * 1e14
  expect_identical(deficient_coords(H), integer(0))
})


test_that("an ordinary correlated pair is not aliased", {
  ## THE GATE IS LOAD-BEARING. K is X'X up to the weights, so it squares the
  ## conditioning of the design and a pivot read at dqrdc2's own tolerance
  ## is twice as strict here as on the augmented system. Measured on two
  ## columns collinear to within 1e-4 -- an ordinary pair of correlated
  ## covariates -- the bare pivot names one of them while solve_pd() inverts
  ## the matrix without difficulty. Nothing is named where it inverts.
  set.seed(5)
  n <- 300
  X <- cbind(1, matrix(stats::runif(n * 3, -1, 1), n, 3))
  near <- function(eps) {
    Y <- X
    Y[, 4] <- Y[, 2] + eps * stats::runif(n, -1, 1)
    crossprod(Y)
  }
  expect_identical(deficient_coords(near(1e-4)), integer(0))
  ## and where it does NOT invert, the column is named rather than the whole
  ## matrix refused
  expect_identical(deficient_coords(near(1e-8)), 4L)
  expect_identical(deficient_coords(near(0)), 4L)
})


test_that("a PENALIZED column is identified and is never aliased", {
  ## The pivot runs on the augmented system, design and penalty factor
  ## together, so a column the design alone does not identify is identified
  ## there. This is the whole reason the rule needs no clause of its own, and
  ## it is the control that keeps the test above from passing vacuously: the
  ## same two identical columns, penalized, must alias NOTHING.
  set.seed(3)
  n <- 200
  dd <- data.frame(x1 = stats::runif(n, -1, 1), x2 = stats::runif(n, -1, 1))
  dd$x3 <- dd$x1
  dd$y <- stats::rnorm(n, 1 + dd$x1 - dd$x2, 0.5)

  f <- statmod(y ~ x2 + ridge(~ 0 + x1 + x3),
               distributions7::gaussian1_distrib(), dd)
  expect_length(f@aliased, 0L)
  expect_false(anyNA(diag(vcov(f))))
  ## and the penalty splits the effect between them rather than dropping one
  b <- coef(f)$mu
  expect_equal(unname(b[["ridge.x1"]]), unname(b[["ridge.x3"]]),
               tolerance = 1e-8)
})


test_that("a factor with two redundant levels aliases the same two as glm", {
  ## more than one alias, and from a factor rather than a duplicated numeric
  set.seed(5)
  dd <- data.frame(g = factor(rep(c("a", "b", "c"), each = 60)))
  dd$h <- factor(ifelse(dd$g == "a", "p", "q"))
  dd$k <- factor(ifelse(dd$g == "c", "q", "p"))
  dd$y <- stats::rnorm(180)
  f <- statmod(y ~ g + h + k, distributions7::gaussian1_distrib(), dd)
  g <- stats::glm(y ~ g + h + k, data = dd)
  expect_identical(f@aliased,
                   paste0("mu:", names(coef(g))[is.na(coef(g))]))
  expect_equal(sum(is.na(coef(f)$mu)), sum(is.na(coef(g))))
})
