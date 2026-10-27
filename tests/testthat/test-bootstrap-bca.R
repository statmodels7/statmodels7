g1c <- distributions7::gaussian1_distrib()

bca_data <- function(n = 40, seed = 4) {
  set.seed(seed)
  d <- data.frame(x = runif(n), id = rep(seq_len(n / 4), each = 4))
  d$y <- 1 + 2 * d$x + rnorm(n, sd = 0.2 + d$x)
  d
}

# The BCa limits written out from the definition: bias correction from the
# share of replicates below the estimate, acceleration from the jackknife.
bca_hand <- function(q, est, jk, alpha) {
  z0 <- qnorm(mean(q < est))
  u <- mean(jk) - jk
  a <- sum(u^3) / (6 * sum(u^2)^1.5)
  z <- qnorm(c(alpha / 2, 1 - alpha / 2))
  p <- pnorm(z0 + (z0 + z) / (1 - a * (z0 + z)))
  quantile(q, p, names = FALSE)
}

test_that("the BCa limits follow the definition", {
  set.seed(1)
  x <- rexp(15)
  q <- replicate(2000, mean(sample(x, replace = TRUE)))
  est <- mean(x)
  jk <- vapply(seq_along(x), function(i) mean(x[-i]), 0)
  got <- statmodels7:::boot_bca_limits(q, est, jk, 0.05)
  expect_equal(as.numeric(got), bca_hand(q, est, jk, 0.05), tolerance = 1e-12)
  # skewed data: the interval moves away from the percentile one
  expect_gt(abs(got[[2L]] - quantile(q, 0.975, names = FALSE)), 1e-4)
  # a symmetric set of replicates and of jackknife values gives the
  # percentile interval
  qs <- c(-3:3, -3:3)
  gs <- statmodels7:::boot_bca_limits(qs, 0, c(-1, 1, -2, 2), 0.1)
  expect_equal(as.numeric(gs),
               quantile(qs, c(0.05, 0.95), names = FALSE), tolerance = 1e-12)
  # an estimate outside the replicates has no bias correction: NA limits
  expect_true(all(is.na(statmodels7:::boot_bca_limits(q, max(q) + 1, jk, 0.05))))
  # ties with the estimate count by halves
  expect_equal(as.numeric(statmodels7:::boot_bca_limits(rep(1, 10), 1,
                                                        c(0, 1, 2), 0.05)),
               c(1, 1))
})

test_that("the refits without one observation equal the leave-one-out formula", {
  d <- bca_data()
  n <- nrow(d)
  f <- statmod(y ~ x, g1c, d, outer_criterion = ml())
  X <- cbind(1, d$x)
  A <- solve(crossprod(X))
  b <- drop(A %*% crossprod(X, d$y))
  r <- drop(d$y - X %*% b)
  h <- rowSums((X %*% A) * X)
  ref <- t(vapply(seq_len(n), function(i) b - drop(A %*% X[i, ]) * r[i] / (1 - h[i]),
                  numeric(2)))
  J <- statmodels7:::bootstrap_jackknife(f, case_resample(), "coefficients")
  expect_equal(dim(J), c(n, 3L))
  expect_equal(unname(J[, 1:2]), ref, tolerance = 1e-5)
  expect_identical(attr(J, "failed"), 0L)
  # a cluster is deleted whole
  Jc <- statmodels7:::bootstrap_jackknife(f, cluster_resample(~ id, data = d),
                                          "coefficients")
  groups <- split(seq_len(n), d$id)
  refc <- t(vapply(groups, function(idx) {
    Xg <- X[idx, , drop = FALSE]
    Hgg <- Xg %*% A %*% t(Xg)
    b - drop(A %*% t(Xg) %*% solve(diag(length(idx)) - Hgg, r[idx]))
  }, numeric(2)))
  expect_equal(dim(Jc), c(length(groups), 3L))
  expect_equal(unname(Jc[, 1:2]), unname(refc), tolerance = 1e-5)
  # the random number stream of the session is left where it was
  set.seed(5)
  s0 <- get(".Random.seed", envir = globalenv())
  statmodels7:::bootstrap_jackknife(f, case_resample(), "coefficients")
  expect_identical(get(".Random.seed", envir = globalenv()), s0)
})

test_that("a weight of zero removes the observation from the units", {
  d <- bca_data(24)
  w <- rep(1, 24)
  w[c(2, 9)] <- 0
  f <- statmod(y ~ x, g1c, d, weights = w, outer_criterion = ml())
  J <- statmodels7:::bootstrap_jackknife(f, case_resample(), "coefficients")
  expect_equal(nrow(J), 22L)
})

test_that("confint with bca_interval equals the hand interval from the stored replicates", {
  d <- bca_data(60, seed = 8)
  f <- statmod(y ~ x, g1c, d, outer_criterion = ml())
  set.seed(2)
  f <- add_bootstrap(f, n = 400, resample = case_resample())
  expect_message(
    ci <- confint(f, interval = bca_interval(n = 400,
                                             resample = case_resample())),
    "refits without one observation")
  X <- cbind(1, d$x)
  A <- solve(crossprod(X))
  b <- drop(A %*% crossprod(X, d$y))
  r <- drop(d$y - X %*% b)
  h <- rowSums((X %*% A) * X)
  JK <- t(vapply(seq_len(60), function(i) b - drop(A %*% X[i, ]) * r[i] / (1 - h[i]),
                 numeric(2)))
  Bm <- coef(statmod_bootstrap(f, resample = case_resample()))
  for (k in 1:2) {
    nm <- c("mu:(Intercept)", "mu:x")[k]
    hand <- bca_hand(Bm[, nm], b[k], JK[, k], 0.05)
    row <- ci[rownames(ci) == nm, ]
    expect_equal(c(row$lower, row$upper), hand, tolerance = 1e-4)
  }
  # the percentile interval of the same entry is different
  cp <- confint(f, interval = percentile_interval(n = 400,
                                                 resample = case_resample()))
  expect_gt(max(abs(ci$lower - cp$lower), abs(ci$upper - cp$upper)), 1e-6)
})

test_that("the refits are stored in the entry and survive an extension", {
  d <- bca_data(32, seed = 9)
  f <- statmod(y ~ x, g1c, d, outer_criterion = ml())
  set.seed(3)
  f <- add_bootstrap(f, n = 120, resample = case_resample())
  orig <- statmodels7:::bootstrap_jackknife
  calls <- 0L
  testthat::local_mocked_bindings(bootstrap_jackknife = function(...) {
    calls <<- calls + 1L
    orig(...)
  })
  bi <- bca_interval(n = 120, resample = case_resample())
  suppressWarnings(suppressMessages(confint(f, interval = bi)))
  expect_no_message(suppressWarnings(confint(f, interval = bi, level = 0.9)))
  expect_equal(calls, 1L)
  ent <- statmod_bootstrap(f, resample = case_resample())
  expect_false(is.null(ent@jackknife))
  f <- add_bootstrap(f, n = 200, resample = case_resample())
  expect_false(is.null(statmod_bootstrap(f, resample = case_resample())@jackknife))
  expect_no_message(suppressWarnings(confint(
    f, interval = bca_interval(n = 200, resample = case_resample()))))
  expect_equal(calls, 1L)
})

test_that("a default BCa interval draws case replicates", {
  d <- bca_data(24, seed = 10)
  f <- statmod(y ~ x, g1c, d, outer_criterion = ml())
  ci <- suppressWarnings(suppressMessages(
    confint(f, interval = bca_interval(n = 200))))
  expect_equal(statmod_bootstrap(f)$resample, "case")
  expect_true(all(is.finite(ci$lower)))
})

test_that("the Bayesian and the cluster schemes are accepted", {
  d <- bca_data(24, seed = 11)
  f <- statmod(y ~ x, g1c, d, outer_criterion = ml())
  cb <- suppressWarnings(suppressMessages(confint(f, interval = bca_interval(
    n = 200, resample = bayesian_resample()))))
  expect_true(all(is.finite(cb$lower)))
  cc <- suppressWarnings(suppressMessages(confint(f, interval = bca_interval(
    n = 200, resample = cluster_resample(~ id, data = d)))))
  expect_true(all(is.finite(cc$lower)))
})

test_that("the BCa interval is rejected where the acceleration is not defined", {
  d <- bca_data(24, seed = 12)
  f <- statmod(y ~ x, g1c, d, outer_criterion = ml())
  expect_error(confint(f, interval = bca_interval(
    n = 50, resample = parametric_resample())), "percentile_interval")
  f <- suppressMessages(add_bootstrap(f, n = 20,
                                     resample = parametric_resample()))
  expect_error(confint(f, interval = bca_interval(
    resample = parametric_resample())), "case_resample")
  gg <- data.frame(g = factor(rep(1:8, each = 5)), x = rnorm(40))
  gg$y <- 1 + gg$x + rnorm(8)[gg$g] + rnorm(40, sd = 0.5)
  fr <- statmod(y ~ x + random(~ 1 | g), g1c, gg)
  expect_error(confint(fr, interval = bca_interval(
    n = 20, resample = cluster_resample(~ g, data = gg))), "random")
  expect_error(confint(f, interval = bca_interval(), vcov_type = bayesian_vcov()),
               "vcov_type")
  expect_error(confint(f, interval = bca_interval(), test = lr_test()),
               "cannot be combined")
  expect_output(print(bca_interval(n = 300)), "bca bootstrap interval")
})
