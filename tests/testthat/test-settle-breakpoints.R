skip_on_cran()

sim_jseg <- function(s, n = 200) {
  set.seed(s)
  x <- sort(runif(n, 0, 10))
  y <- 1 + 0.3 * x + ifelse(x > 6, 1.5 - 0.6 * (x - 6), 0) + rnorm(n, sd = 0.5)
  data.frame(x = x, y = y)
}

test_that("a sharp jseg is settled at the profile minimum and certified", {
  # seed 1: the working iteration stopped one interval away from the minimum
  d <- sim_jseg(1)
  f <- statmod(y ~ jseg(x), distrib = distributions7::gaussian1_distrib(),
               data = d)
  prof <- function(q) {
    sum(stats::lm.fit(cbind(1, d$x, d$x > q, pmax(d$x - q, 0)), d$y)$residuals^2)
  }
  u <- sort(unique(d$x))
  m <- (u[-1] + u[-length(u)]) / 2
  m <- m[m > stats::quantile(d$x, 0.05) & m < stats::quantile(d$x, 0.95)]
  best <- min(vapply(m, prof, 1))
  rss <- sum((d$y - fitted(f, what = "mu"))^2)
  expect_equal(rss, best, tolerance = 1e-8)
  expect_identical(statmod_certificate(f)$state, "converged")
  expect_length(f@aliased, 0L)
})

test_that("the coefficients of a held jseg are least squares at the position", {
  d <- sim_jseg(3)
  f <- statmod(y ~ jseg(x), distrib = distributions7::gaussian1_distrib(),
               data = d)
  p <- coef(f)$mu[["jseg.psi1"]]
  l <- stats::lm(y ~ x + I(x > p) + pmax(x - p, 0), data = d)
  cf <- coef(f)$mu
  expect_equal(unname(cf[c("(Intercept)", "jseg.beta", "jseg.delta1",
                           "jseg.gamma1")]), unname(coef(l)), tolerance = 1e-8)
  # sigma by REML, as lm() reads it, and the standard errors with it
  expect_equal(exp(coef(f)$sigma[[1]]), summary(l)$sigma, tolerance = 1e-6)
  ci <- confint(f)
  se_l <- sqrt(diag(stats::vcov(l)))
  expect_equal(unname(ci[c("mu:jseg.beta", "mu:jseg.delta1", "mu:jseg.gamma1"), "se"]),
               unname(se_l[2:4]), tolerance = 1e-5)
  # the position carries no standard error and counts one degree of freedom
  expect_true(is.na(ci["mu:jseg.psi1", "se"]))
  expect_equal(sum(f@edf$edf), 6)
  # and vcov() holds it without a warning
  expect_no_warning(V <- vcov(f))
  expect_true(all(is.na(V["mu:jseg.psi1", ])))
})

test_that("a jump beside a smooth is certified after the criterion is re-read", {
  set.seed(3)
  n <- 300
  x <- runif(n, 0, 10)
  z <- runif(n)
  y <- 1 + 1.5 * (x > 6) + sin(2 * pi * z) + rnorm(n, sd = 0.5)
  f <- statmod(y ~ jump(x) + s(z), distrib = distributions7::gaussian1_distrib(),
               data = data.frame(x, y, z))
  expect_identical(statmod_certificate(f)$state, "converged")
})

test_that("a sharp jseg in sigma's equation is settled at the profile's minimum", {
  skip_if_not_installed("MASS")
  data(mcycle, package = "MASS")
  f <- statmod(accel ~ s(times, bspline_smooth(k = 20)) | sigma ~ jseg(times),
               distrib = distributions7::gaussian1_distrib(), data = mcycle)
  p <- coef(f)$sigma[["jseg.psi1"]]
  # the global minimum of the REML profile over every interval is at 14.7,
  # between the observed times 14.6 and 14.8
  expect_equal(p, 14.7, tolerance = 1e-10)
  expect_identical(statmod_certificate(f)$state, "converged")
  # held, the fit is the model with the break-point's columns written out
  d <- mcycle
  d$h <- pmax(d$times - p, 0)
  d$s <- as.numeric(d$times > p)
  g <- statmod(accel ~ s(times, bspline_smooth(k = 20)) | sigma ~ times + h + s,
               distrib = distributions7::gaussian1_distrib(), data = d)
  expect_equal(f@criterion, g@criterion, tolerance = 1e-6)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(g)), tolerance = 1e-6)
})

test_that("a sharp jump in sigma's equation is certified", {
  skip_if_not_installed("MASS")
  data(mcycle, package = "MASS")
  f <- statmod(accel ~ s(times, bspline_smooth(k = 20)) | sigma ~ jump(times),
               distrib = distributions7::gaussian1_distrib(), data = mcycle)
  expect_identical(statmod_certificate(f)$state, "converged")
  expect_true(all(is.finite(f@coefficients$sigma)))
})

test_that("a jseg developed over groups is settled at the per-group minimum", {
  set.seed(1)
  g <- factor(rep(c("a", "b", "c"), each = 30))
  x <- runif(90, 0, 10)
  gi <- as.integer(g)
  pg <- c(3, 5, 7)[gi]
  y <- 1 + 0.3 * x + ifelse(x > pg, 1.5 - 0.6 * (x - pg), 0) +
    rnorm(90, sd = 0.5)
  d <- data.frame(x = x, g = g, y = y)
  f <- statmod(y ~ jseg(x, psi ~ 0 + g), distrib = distributions7::gaussian1_distrib(),
               data = d)
  rss <- function(q) {
    q <- q[gi]
    sum(stats::.lm.fit(cbind(1, x, pmax(x - q, 0), x > q), y)$residuals^2)
  }
  cand <- lapply(1:3, function(j) {
    u <- sort(unique(x[gi == j]))
    m <- (u[-1] + u[-length(u)]) / 2
    lim <- stats::quantile(x[gi == j], c(0.05, 0.95))
    m[m > lim[1] & m < lim[2]]
  })
  best <- Inf
  for (a in cand[[1]]) for (b in cand[[2]]) {
    best <- min(best, vapply(cand[[3]], function(v) rss(c(a, b, v)), 1))
  }
  # the polish also minimizes inside the interval, so it can only be lower
  expect_lte(sum((y - fitted(f, what = "mu"))^2), best + 1e-8)
  expect_identical(statmod_certificate(f)$state, "converged")
  cf <- coef(f)$mu
  expect_true(all(c("jseg.psi1.ga", "jseg.psi1.gb", "jseg.psi1.gc") %in% names(cf)))
  ci <- confint(f)
  expect_true(all(is.na(ci[c("mu:jseg.psi1.ga", "mu:jseg.psi1.gb"), "se"])))
})
