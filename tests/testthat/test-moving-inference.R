# Leverage, influence and the predictive bootstrap beside a block that moves
# with its coefficients: the hat matrix is the tangent-plane leverage, read
# on the Jacobian at the estimate that the fit stores.

mi_data <- function(seed = 1, n = 150) {
  set.seed(seed)
  d <- data.frame(x = runif(n))
  d$y <- 2 * (1 - exp(-3 * d$x)) + rnorm(n, sd = 0.2)
  d
}

test_that("the leverage of nl() is the tangent-plane leverage", {
  d <- mi_data()
  fit <- statmod(y ~ 0 + modelterms7::nl(~ a * (1 - exp(-r * x)), start = list(a = 1, r = 1)),
                 distrib = distributions7::gaussian1_distrib(), data = d)
  cf <- coef(fit)$mu
  a <- cf[[which(endsWith(names(cf), ".a"))]]
  r <- cf[[which(endsWith(names(cf), ".r"))]]
  mu <- function(a, r) a * (1 - exp(-r * d$x))
  h <- 1e-6
  J <- cbind((mu(a + h, r) - mu(a - h, r)) / (2 * h), (mu(a, r + h) - mu(a, r - h)) / (2 * h))
  lv <- statmod_leverage(fit)
  expect_equal(lv$leverage[, "mu"], rowSums((J %*% solve(crossprod(J))) * J), tolerance = 1e-8)
  # two coefficients in mu, the intercept of sigma in the other equation
  expect_equal(sum(lv$leverage[, "mu"]), 2, tolerance = 1e-8)
  expect_equal(sum(lv$trace), 3, tolerance = 1e-8)
})

test_that("the leverage of scs() reads the Jacobian of its contribution", {
  d <- mi_data()
  fit <- statmod(y ~ modelterms7::scs(x, modelterms7::shape_increasing_concave()),
                 distrib = distributions7::gaussian1_distrib(), data = d)
  tm <- Filter(function(t) S7::S7_inherits(t, modelterms7::ScsTerm), fit@spec@terms$mu)[[1L]]
  cf <- coef(fit)$mu[modelterms7::term_coef_names(tm)]
  p <- length(cf); h <- 1e-6
  Jt <- sapply(seq_len(p), function(j) {
    e <- replace(numeric(p), j, h)
    (modelterms7::term_value(tm, coef = cf + e) - modelterms7::term_value(tm, coef = cf - e)) / (2 * h)
  })
  Z <- cbind(1, Jt)
  s2 <- exp(2 * coef(fit)$sigma[[1L]])
  hy <- unlist(fit@hyper)
  S <- matrix(0, p + 1L, p + 1L)
  S[-1L, -1L] <- as.matrix(penalties7::penalty_hessian(modelterms7::term_penalty(tm), cf,
                                                       list(lambda = hy[[grep("lambda", names(hy))]])))
  K <- crossprod(Z) / s2 + S
  lv <- statmod_leverage(fit)
  expect_equal(lv$leverage[, "mu"], rowSums((Z %*% solve(K)) * Z) / s2, tolerance = 1e-7)
  # the influence graph reads the same hat matrix
  inf <- influence_deltas(fit, "mu")
  expect_equal(inf$leverage, lv$leverage[, "mu"], tolerance = 1e-8)
})

test_that("the predictive bootstrap refits a block that moves", {
  d <- mi_data(n = 100)
  fit <- statmod(y ~ 0 + modelterms7::nl(~ a * (1 - exp(-r * x)), start = list(a = 1, r = 1)),
                 distrib = distributions7::gaussian1_distrib(), data = d)
  nd <- data.frame(x = c(0.1, 0.5, 0.9))
  set.seed(2)
  pb <- predict(fit, what = "response", newdata = nd,
                interval = prediction_interval(predictive = pred_bootstrap(n = 30)))
  pa <- predict(fit, what = "response", newdata = nd, interval = prediction_interval())
  expect_true(all(is.finite(c(pb$lower, pb$upper))))
  expect_lt(max(abs(c(pb$lower - pa$lower, pb$upper - pa$upper))), 0.15)
})

test_that("the leverage of seg() counts its break-point as a coefficient", {
  set.seed(3)
  d <- data.frame(x = runif(200))
  d$y <- 1 + 2 * pmax(d$x - 0.5, 0) + rnorm(200, sd = 0.2)
  fit <- suppressWarnings(statmod(y ~ modelterms7::seg(x),
                                  distrib = distributions7::gaussian1_distrib(), data = d))
  lv <- statmod_leverage(fit)
  expect_equal(sum(lv$leverage[, "mu"]), length(coef(fit)$mu), tolerance = 1e-6)
})
