fsw <- function(sw = NULL, information = NULL) {
  frequentist_vcov(information = information, sandwich = sw)
}

sandwich_hand <- function(V, U, factor = 1) {
  factor * V %*% crossprod(U) %*% V
}

test_that("the default frequentist variance does not change", {
  set.seed(1)
  dd <- data.frame(x = runif(80))
  dd$y <- 1 + 2 * dd$x + rnorm(80, sd = 0.4)
  f <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
  expect_null(frequentist_vcov()@sandwich)
  expect_identical(vcov(f, vcov_type = frequentist_vcov()),
                   vcov(f, vcov_type = frequentist_vcov(sandwich = NULL)))
  # V_f = Vb H Vb written out
  H <- statmod_information_at(f@spec, f@coefficients, statmod_design(f@spec),
                              FALSE, NULL)
  S <- statmod_penalty_at(f@spec, f@coefficients, f@hyper,
                          statmod_design(f@spec), "hessian")
  Vb <- solve(as.matrix(H + S))
  expect_equal(unname(vcov(f, vcov_type = frequentist_vcov(), readable = FALSE)),
               unname(Vb %*% as.matrix(H) %*% Vb), tolerance = 1e-12)
  expect_error(frequentist_vcov(sandwich = "HC0"), "sandwich")
  expect_error(bayesian_vcov(sandwich = hc0_sandwich()), "unused argument")
  expect_error(unconditional_vcov(sandwich = hc0_sandwich()),
               "unused argument")
})

test_that("HC0 and HC1 of a Gaussian linear model equal the hand formulas", {
  set.seed(2)
  n <- 200
  x <- runif(n)
  dd <- data.frame(x = x, y = 1 + 2 * x + rnorm(n, sd = 0.1 + 1.5 * x))
  f <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
  X <- cbind(1, x)
  r <- dd$y - drop(X %*% f@coefficients$mu)
  A <- solve(crossprod(X))
  hc0 <- A %*% crossprod(X * r) %*% A
  v0 <- vcov(f, vcov_type = fsw(hc0_sandwich()), param = "mu", readable = FALSE)
  v1 <- vcov(f, vcov_type = fsw(hc1_sandwich()), param = "mu", readable = FALSE)
  expect_equal(unname(v0), unname(hc0), tolerance = 1e-10)
  expect_equal(unname(v1), unname(hc0 * n / (n - 3)), tolerance = 1e-10)
  # the full matrix is symmetric and positive semi-definite
  V <- vcov(f, vcov_type = fsw(hc1_sandwich()), readable = FALSE)
  expect_equal(V, t(V))
  expect_gte(min(eigen(V, symmetric = TRUE)$values), 0)
  # the other methods that take a vcov_type
  ci <- confint(f, vcov_type = fsw(hc1_sandwich()))
  expect_true(all(is.finite(as.matrix(ci[, sapply(ci, is.numeric)]))))
  s <- summary(f, vcov_type = fsw(hc1_sandwich()))
  expect_no_error(capture.output(print(s)))
  pr <- predict(f, what = "mu", newdata = data.frame(x = c(0.2, 0.8)),
                se = TRUE, vcov_type = fsw(hc1_sandwich()))
  Xn <- cbind(1, c(0.2, 0.8))
  expect_equal(unname(pr$se),
               sqrt(rowSums((Xn %*% (hc0 * n / (n - 3))) * Xn)),
               tolerance = 1e-8)
})

test_that("Poisson and logistic models equal J^-1 (sum u u') J^-1", {
  set.seed(3)
  n <- 300
  x <- rnorm(n)
  dp <- data.frame(x = x, y = rpois(n, exp(0.5 + 0.4 * x) *
                                      rgamma(n, 2, 2)))
  fp <- statmod(y ~ x, distributions7::poisson_distrib(), dp)
  X <- cbind(1, x)
  mu <- exp(drop(X %*% fp@coefficients$mu))
  J <- crossprod(X * sqrt(mu))
  Ji <- solve(J)
  hc0 <- Ji %*% crossprod(X * (dp$y - mu)) %*% Ji
  v <- vcov(fp, vcov_type = fsw(hc0_sandwich()), readable = FALSE)
  expect_equal(unname(v), unname(hc0), tolerance = 1e-9)

  pr <- plogis(-0.3 + 0.8 * x)
  db <- data.frame(x = x, y = rbinom(n, 1, pr))
  fb <- statmod(y ~ x, distributions7::bernoulli_distrib(), db)
  p1 <- plogis(drop(X %*% fb@coefficients$mu))
  Jb <- solve(crossprod(X * sqrt(p1 * (1 - p1))))
  hb <- Jb %*% crossprod(X * (db$y - p1)) %*% Jb
  vb <- vcov(fb, vcov_type = fsw(hc0_sandwich()), readable = FALSE)
  expect_equal(unname(vb), unname(hb), tolerance = 1e-9)
  v1 <- vcov(fb, vcov_type = fsw(hc1_sandwich()), readable = FALSE)
  expect_equal(unname(v1), unname(hb * n / (n - 2)), tolerance = 1e-9)
})

test_that("clustered HC0 and HC1 are CR0 and CR1", {
  set.seed(4)
  G <- 25
  n <- 100
  id <- rep(seq_len(G), each = 4)
  x <- rnorm(n)
  y <- 1 + x + rnorm(G)[id] + rnorm(n)
  dd <- data.frame(x = x, y = y, id = id, row = seq_len(n))
  f <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
  X <- cbind(1, x)
  r <- y - drop(X %*% f@coefficients$mu)
  A <- solve(crossprod(X))
  Ug <- rowsum(X * r, id)
  cr0 <- A %*% crossprod(Ug) %*% A
  v0 <- vcov(f, vcov_type = fsw(hc0_sandwich(cluster = ~ id)),
             param = "mu", readable = FALSE)
  v1 <- vcov(f, vcov_type = fsw(hc1_sandwich(cluster = ~ id)),
             param = "mu", readable = FALSE)
  expect_equal(unname(v0), unname(cr0), tolerance = 1e-10)
  expect_equal(unname(v1), unname(cr0 * G / (G - 1) * (n - 1) / (n - 3)),
               tolerance = 1e-10)
  # a cluster for each row is the HC estimator, to the bit of the reduction
  for (sw in list(hc0_sandwich, hc1_sandwich)) {
    a <- vcov(f, vcov_type = fsw(sw(cluster = ~ row)), readable = FALSE)
    b <- vcov(f, vcov_type = fsw(sw()), readable = FALSE)
    expect_equal(a, b, tolerance = 1e-12)
  }
  # an explicit data frame replaces the data of the call
  v2 <- vcov(f, vcov_type = fsw(hc1_sandwich(cluster = ~ id, data = dd)),
             param = "mu", readable = FALSE)
  expect_identical(v1, v2)
})

test_that("the location-scale Gaussian sandwich equals the hand one", {
  set.seed(5)
  n <- 400
  x <- runif(n)
  y <- 1 + 2 * x + exp(-0.5 + 0.8 * x) * rt(n, 5) / sqrt(5 / 3)
  dd <- data.frame(x = x, y = y)
  f <- statmod(y ~ x | sigma ~ x, distributions7::gaussian1_distrib(), dd,
               outer_criterion = ml())
  mu <- predict(f, what = "mu")
  s <- predict(f, what = "sigma")
  X <- cbind(1, x)
  r <- (y - mu) / s
  U <- cbind(X * (r / s), X * (r^2 - 1))
  V <- vcov(f, vcov_type = fsw(), readable = FALSE)
  Vb <- vcov(f, vcov_type = bayesian_vcov(), readable = FALSE)
  # the information is the observed one, so the model variance is Vb itself
  Vs <- vcov(f, vcov_type = fsw(hc0_sandwich()), readable = FALSE)
  expect_equal(unname(Vs), unname(sandwich_hand(Vb, U)), tolerance = 1e-8)
  # on the readable scale the sigma coefficients are on the log scale already
  Vr <- vcov(f, vcov_type = fsw(hc0_sandwich()))
  expect_equal(unname(sqrt(diag(Vr))), unname(sqrt(diag(Vs))), tolerance = 1e-8)
})

test_that("a model with a smooth and a random effect runs", {
  set.seed(6)
  G <- 15
  n <- 150
  id <- rep(seq_len(G), each = 10)
  x <- runif(n)
  y <- sin(2 * pi * x) + rnorm(G, sd = 0.5)[id] + rnorm(n, sd = 0.3)
  dd <- data.frame(x = x, y = y, id = factor(id))
  f <- suppressWarnings(statmod(
    y ~ s(x, basis7::bspline_smooth(k = 8)) + random(~ 1 | id),
    distributions7::gaussian1_distrib(), dd))
  for (sw in list(hc0_sandwich(), hc1_sandwich(),
                  hc1_sandwich(cluster = ~ id))) {
    V <- vcov(f, vcov_type = fsw(sw), readable = FALSE)
    expect_equal(V, t(V), tolerance = 1e-12)
    ok <- !is.na(diag(V))
    expect_gte(min(eigen(V[ok, ok], symmetric = TRUE)$values), -1e-10)
    expect_true(all(diag(V)[ok] >= 0))
  }
  # the model-meat path is the frequentist matrix, to the bit
  expect_identical(vcov(f, vcov_type = fsw(NULL)),
                   vcov(f, vcov_type = frequentist_vcov()))
  # a cluster that splits a random-effect group is rejected
  dd$half <- rep(rep(1:2, each = 5), G)
  f2 <- suppressWarnings(statmod(
    y ~ s(x, basis7::bspline_smooth(k = 8)) + random(~ 1 | id),
    distributions7::gaussian1_distrib(), dd))
  expect_error(vcov(f2, vcov_type = fsw(hc1_sandwich(cluster = ~ half))),
               "more than one cluster")
})

test_that("unsupported fits and bad clusters are rejected", {
  set.seed(7)
  dd <- data.frame(x = runif(60), id = rep(1:12, each = 5))
  dd$y <- 1 + dd$x + rnorm(60)
  f <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
  expect_error(vcov(f, vcov_type = fsw(hc1_sandwich(cluster = ~ nope))),
               "not a column")
  dd2 <- dd
  dd2$id[3] <- NA
  expect_error(vcov(f, vcov_type = fsw(hc1_sandwich(cluster = ~ id, data = dd2))),
               "missing values")
  expect_error(hc1_sandwich(cluster = id ~ x), "one-sided")
  expect_error(hc1_sandwich(cluster = ~ id + x), "one-sided")
  # prior weights other than one
  fw <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd,
                weights = rep(c(1, 2), 30))
  expect_error(vcov(fw, vcov_type = fsw(hc1_sandwich())), "prior weights")
  # a structural term
  fg <- suppressWarnings(statmod(y ~ x + gas(p = 1, q = 1),
                                 distributions7::gaussian1_distrib(), dd))
  expect_error(vcov(fg, vcov_type = fsw(hc1_sandwich())), "structural")
  # a term whose block moves
  set.seed(8)
  x <- sort(runif(120))
  ds <- data.frame(x = x, y = ifelse(x < 0.5, x, 0.5 + 3 * (x - 0.5)) +
                     rnorm(120, sd = 0.05))
  fs <- suppressWarnings(statmod(y ~ seg(x, psi ~ 1),
                                 distributions7::gaussian1_distrib(), ds))
  expect_error(vcov(fs, vcov_type = fsw(hc1_sandwich())), "moves with")
})

test_that("a lasso coordinate at zero is missing, as in the other conventions", {
  set.seed(9)
  n <- 120
  dd <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  dd$y <- 1 + dd$x1 + rnorm(n)
  f <- suppressWarnings(statmod(
    y ~ lasso(~ x1 + x2 + x3, lambda = 50), distributions7::gaussian1_distrib(), dd))
  Vf <- vcov(f, vcov_type = frequentist_vcov(), readable = FALSE)
  Vs <- vcov(f, vcov_type = fsw(hc1_sandwich()), readable = FALSE)
  expect_identical(is.na(diag(Vf)), is.na(diag(Vs)))
  expect_true(anyNA(diag(Vs)))
})
