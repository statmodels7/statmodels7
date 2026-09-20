# A penalty with a kink over a structural term's own parameters. Those
# parameters are not columns of any design and not entries of the stacked
# coefficient vector, so neither route sparse_fit() takes for a block of
# coefficients can address them: the coordinate descent reads the positions as
# columns and the proximal branch reads beta[index], which is empty. The fit
# died on "subscript out of bounds" in the published 0.124.0 and at 0.140.0.

gas_panel <- function(seed = 1L, m = 8L, ni = 40L,
                      om = c(0.9, 0.8, 0.85, 0, 0, 0, 0, 0)) {
  set.seed(seed)
  n <- m * ni
  id <- factor(rep(seq_len(m), each = ni))
  y <- numeric(n)
  for (g in seq_len(m)) {
    fl <- 0
    for (t in seq_len(ni)) {
      i <- (g - 1L) * ni + t
      mu <- om[[g]] + fl
      y[[i]] <- stats::rnorm(1, mu, 1)
      fl <- 0.4 * (y[[i]] - mu) + 0.6 * fl
    }
  }
  data.frame(y = y, id = id)
}

fit_at <- function(dd, lambda) {
  suppressWarnings(statmod(
    y ~ 0 + gas(p = 1, q = 1, by = id,
                omega ~ 0 + lasso(~ id, lambda = lambda)),
    distrib = distributions7::gaussian1_distrib(), data = dd,
    outer_criterion = NULL))
}

omega_of <- function(f) {
  z <- f@structural[[1L]]$parameter
  z[grepl("^omega", names(z))]
}

test_that("a lasso held over a filter's own parameters fits", {
  skip_on_cran()
  dd <- gas_panel()
  f <- fit_at(dd, 1)
  expect_true(is.finite(as.numeric(stats::logLik(f))))
  expect_identical(length(omega_of(f)), 8L)
})

test_that("it selects, and the survivors are the levels that are there", {
  skip_on_cran()
  dd <- gas_panel()
  z5 <- omega_of(fit_at(dd, 5))
  z20 <- omega_of(fit_at(dd, 20))
  z80 <- omega_of(fit_at(dd, 80))
  # more shrinkage, never fewer zeros
  expect_true(sum(z5 == 0) <= sum(z20 == 0))
  expect_true(sum(z20 == 0) <= sum(z80 == 0))
  expect_identical(sum(z80 == 0), 8L)
  # the three groups that carry a level are the last to go
  expect_true(all(abs(z20[1:3]) > abs(z20[6:8])))
})

test_that("the point is the penalized mode, by a route that shares no arithmetic", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  dd <- gas_panel()
  lam <- 5
  fm <- y ~ 0 + gas(p = 1, q = 1, by = id,
                    omega ~ 0 + lasso(~ id, lambda = lam))
  f <- suppressWarnings(statmod(fm, distributions7::gaussian1_distrib(), dd,
                                outer_criterion = NULL))
  # the smooth objective rebuilt as a function of the free structural
  # parameters and differenced: neither the filter's adjoint nor the proximal
  # operator is involved
  spec <- statmod_spec(fm, distributions7::gaussian1_distrib(), dd)
  des <- statmod_design(spec)
  sst <- statmod_structural_state(des)
  key <- names(sst$zeta)[[1L]]
  nm <- names(sst$zeta[[key]])
  free <- setdiff(nm, sst$held[[key]])
  cf <- f@coefficients
  zhat <- as.numeric(f@structural[[1L]]$unconstrained[free])
  negll <- function(z) {
    v <- sst$zeta[[key]]
    v[free] <- z
    sst$zeta[[key]] <- v
    sst$key <- NULL
    sst$value <- NULL
    -statmod_loglik_at(spec, cf, des)
  }
  g <- numDeriv::grad(negll, zhat)
  pen <- which(grepl("^omega", free))
  b <- zhat[pen]
  gp <- g[pen]
  away <- abs(b) > 1e-8
  expect_true(any(away))
  # stationarity where the coefficient is away from the kink
  expect_lt(max(abs(gp[away] + lam * sign(b[away]))), 1e-5)
  # and inside the interval the kink opens where it is at zero
  if (any(!away)) expect_lte(max(abs(gp[!away])), lam * (1 + 1e-6))
})

test_that("the block the layer builds is marked as structural", {
  skip_on_cran()
  dd <- gas_panel()
  spec <- statmod_spec(
    y ~ 0 + gas(p = 1, q = 1, by = id,
                omega ~ 0 + lasso(~ id, lambda = 5)),
    distributions7::gaussian1_distrib(), dd)
  des <- statmod_design(spec)
  key <- names(statmod_structural_state(des)$zeta)[[1L]]
  blocks <- statmod_blocks(spec, des)
  expect_identical(length(blocks$sparse), 1L)
  expect_true(isTRUE(blocks$sparse[[1L]]$structural))
  # `zterm` is the TERM's name, which the structural state is keyed by, where
  # `term` is the penalty's key and would carry `::entry` for a term
  # declaring more than one
  expect_identical(blocks$sparse[[1L]]$zterm, key)
  expect_null(blocks$sparse[[1L]]$index)
  # the positions are among the term's own parameters, and they are the
  # levels the lasso covers and not the loading or the persistence
  nm <- names(statmod_structural_state(des)$zeta[[key]])
  expect_true(all(grepl("^omega", nm[blocks$sparse[[1L]]$cols])))
})

test_that("a smooth penalty there is untouched by the new route", {
  skip_on_cran()
  # the control: a ridge over the same parameters has no kink, so it is not a
  # sparse block at all and the smooth structural fit carries it as before
  dd <- gas_panel()
  f <- suppressWarnings(statmod(
    y ~ 0 + gas(p = 1, q = 1, by = id,
                omega ~ 0 + ridge(~ id, lambda = 20)),
    distributions7::gaussian1_distrib(), dd, outer_criterion = NULL))
  expect_true(f@converged)
  expect_true(all(omega_of(f) != 0))
})
