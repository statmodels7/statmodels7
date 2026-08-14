# Standardization end to end. A separable penalty is standardized by a
# DIAGONAL MAP on the penalty and never by rescaling the design, so what has
# to be checked here is that the fit is the one a hand-standardized design
# gives, that it no longer depends on the units a column was measured in, and
# that a sparse block in the same equation survives the whole chain.

set.seed(21)
n <- 200
sd_dat <- data.frame(
  x1 = rnorm(n), x2 = rnorm(n) * 100, x3 = rnorm(n) / 100, x4 = rnorm(n))
sd_dat$y <- with(sd_dat, 1 + 1.5 * x1 + 0.02 * x2 + rnorm(n))
# the hyperparameter is held IN THE TERM now, so every fit below writes it
sd_cols <- c("x1", "x2", "x3", "x4")

test_that("a standardized term fits what a hand-standardized design fits", {
  f_std <- statmod(y ~ lasso(~ x1 + x2 + x3 + x4, standardize = TRUE,
                             lambda = 2),
                   distributions7::gaussian1_distrib(), sd_dat)
  s <- vapply(sd_dat[sd_cols], stats::sd, numeric(1))
  dz <- data.frame(y = sd_dat$y, z1 = sd_dat$x1 / s[1], z2 = sd_dat$x2 / s[2],
                   z3 = sd_dat$x3 / s[3], z4 = sd_dat$x4 / s[4])
  f_man <- statmod(y ~ lasso(~ z1 + z2 + z3 + z4, lambda = 2),
                   distributions7::gaussian1_distrib(), dz)
  # the coefficients are compared BY POSITION: the intercept is first and the
  # four slopes follow in the order the formula names them
  bx <- unname(coef(f_std)[[1L]])
  bz <- unname(coef(f_man)[[1L]])
  expect_length(bx, 5L)
  expect_equal(bz[-1L], unname(s) * bx[-1L])
  expect_equal(bz[1L], bx[1L])
  expect_equal(fitted(f_std)[[1L]], fitted(f_man)[[1L]])
})

test_that("standardizing makes the fit independent of a column's units", {
  # the test that says whether it buys anything: multiplying a column by a
  # thousand must leave a standardized fit where it was, and must NOT leave an
  # unstandardized one there
  scaled <- sd_dat
  scaled$x1 <- sd_dat$x1 * 1000
  fit <- function(dat, std) {
    statmod(y ~ lasso(~ x1 + x2 + x3 + x4, standardize = std, lambda = 2),
            distributions7::gaussian1_distrib(), dat)
  }
  a <- fit(sd_dat, TRUE)
  b <- fit(scaled, TRUE)
  expect_equal(fitted(b)[[1L]], fitted(a)[[1L]])
  # and the coefficient carries the rescaling exactly
  expect_equal(unname(coef(b)[[1L]])[2L] * 1000, unname(coef(a)[[1L]])[2L])

  a0 <- fit(sd_dat, FALSE)
  b0 <- fit(scaled, FALSE)
  expect_gt(max(abs(fitted(b0)[[1L]] - fitted(a0)[[1L]])), 1e-3)
})

test_that("a sparse block survives a standardized penalty in the same equation", {
  # a grouping indicator is built sparse, and standardization works on the
  # penalty rather than on the design, so nothing in the chain has an excuse
  # to densify it
  set.seed(3)
  m <- 60
  dr <- data.frame(g = factor(rep(seq_len(m), each = 6)))
  dr$v1 <- rnorm(nrow(dr))
  dr$v2 <- rnorm(nrow(dr)) * 50
  dr$y <- rnorm(nrow(dr), 1 + 0.5 * dr$v1)
  fml <- y ~ lasso(~ v1 + v2, standardize = TRUE) + random(~ 1 | g)
  spec <- statmod_spec(fml, distributions7::gaussian1_distrib(), dr)
  blocks <- lapply(spec@terms[["mu"]], modelterms7::term_matrix)
  kinds <- vapply(blocks, function(X) class(X)[1L], character(1))
  expect_true("dgCMatrix" %in% kinds)

  design <- statmod_design(spec)
  X <- design[["mu"]]$X
  expect_s4_class(X, "dgCMatrix")
  expect_lt(Matrix::nnzero(X) / prod(dim(X)), 0.1)

  # and the fit runs on it: the compiled coordinate descent takes an
  # arma::mat, so a sparse block reaching the .Call aborted the process
  # rather than raising, which is why this is asserted end to end
  # the lambda is HELD here, and that is the point rather than a detail: with
  # it estimated the criterion zeroes a column, and a coefficient at the kink
  # has no variance to report, so the finiteness below would fail for a
  # reason that has nothing to do with sparsity
  fml2 <- y ~ lasso(~ v1 + v2, standardize = TRUE, lambda = 1) +
    random(~ 1 | g)
  fit <- statmod(fml2, distributions7::gaussian1_distrib(), dr)
  expect_true(fit@converged)
  expect_s4_class(statmod_information_at(spec, fit@coefficients, design,
                                         FALSE, "bartlett"), "dgCMatrix")
  expect_true(all(is.finite(diag(vcov(fit)))))
})

test_that("the sparse kernel is the dense one, bit for bit", {
  # One algorithm instantiated twice over a column accessor, so the claim is
  # IDENTITY and not a tolerance: skipping a structural zero omits an
  # addition of zero, which is exact. A reordering of the arithmetic would
  # fail this, which is the point of asserting it this way.
  set.seed(17)
  n <- 800L; p <- 30L
  for (cv in c(FALSE, TRUE)) {
    for (dens in c(0.01, 0.2)) {
      X <- matrix(0, n, p)
      nz <- round(n * p * dens)
      X[sample(n * p, nz)] <- rnorm(nz)
      S <- coord_block(methods::as(methods::as(
        methods::as(X, "dMatrix"), "generalMatrix"), "CsparseMatrix"),
        seq_len(p))
      expect_s4_class(S, "dgCMatrix")
      z <- rnorm(n); w <- runif(n, 0.5, 2); b0 <- numeric(p)
      v <- as.numeric(crossprod(w, X^2))
      tab <- penalties7::penalty_prox_spec(
        penalties7::lasso_penalty(n_coef = p), list(lambda = 3), 1 / v)
      scr <- as.integer(seq_len(p) - 1L)
      a <- coord_call(X, z, w, b0, tab, scr, 1e-10, cv)
      b <- coord_call(S, z, w, b0, tab, scr, 1e-10, cv)
      lab <- sprintf("covariance=%s density=%g", cv, dens)
      expect_identical(a$beta, b$beta, label = lab)
      expect_identical(a$sweeps, b$sweeps, label = lab)
      # and the fit is not the empty one, or the identity is vacuous
      expect_gt(sum(a$beta != 0), 0)
    }
  }
  # the screened route computes a gradient the unscreened one does not, so
  # that pass needs its own comparison
  X <- matrix(0, n, p); nz <- round(n * p * 0.05)
  X[sample(n * p, nz)] <- rnorm(nz)
  S <- coord_block(methods::as(methods::as(
    methods::as(X, "dMatrix"), "generalMatrix"), "CsparseMatrix"), seq_len(p))
  z <- rnorm(n); w <- runif(n, 0.5, 2); b0 <- numeric(p)
  v <- as.numeric(crossprod(w, X^2))
  tab <- penalties7::penalty_prox_spec(
    penalties7::lasso_penalty(n_coef = p), list(lambda = 3), 1 / v)
  scr <- as.integer(0:14)
  a <- coord_call(X, z, w, b0, tab, scr, 1e-10, FALSE)
  b <- coord_call(S, z, w, b0, tab, scr, 1e-10, FALSE)
  expect_identical(a$beta, b$beta)
  expect_identical(a$grad, b$grad)
})

test_that("a sparse equation fits under every penalized term", {
  # the same composition for each kind: the smooth branch and the kinked one
  # reach the block by different routes
  set.seed(4)
  m <- 30
  dr <- data.frame(g = factor(rep(seq_len(m), each = 5)))
  dr$v1 <- rnorm(nrow(dr))
  dr$v2 <- rnorm(nrow(dr)) * 20
  dr$y <- rnorm(nrow(dr), 1 + 0.5 * dr$v1)
  # the hyperparameters are HELD, which is what this test is about: whether a
  # standardized penalty reaches a sparse block at all, not how its
  # hyperparameters are chosen. Left free they are selected along a path --
  # under the product, 125 fits per term for an answer nothing here reads.
  # The selection is covered in test-path.R.
  cases <- list(
    ridge = quote(ridge(~ v1 + v2, standardize = TRUE, lambda = 1)),
    lasso = quote(lasso(~ v1 + v2, standardize = TRUE, lambda = 1)),
    scad  = quote(scad(~ v1 + v2, standardize = TRUE, lambda = 1, a = 3.7)),
    mcp   = quote(mcp(~ v1 + v2, standardize = TRUE, lambda = 1, gamma = 3))
  )
  for (nm in names(cases)) {
    fml <- stats::as.formula(sprintf(
      "y ~ %s + random(~ 1 | g)",
      paste(deparse(cases[[nm]]), collapse = "")))
    fit <- statmod(fml, distributions7::gaussian1_distrib(), dr)
    expect_true(fit@converged, label = nm)
    expect_true(all(is.finite(unlist(coef(fit)))), label = nm)
  }
})
