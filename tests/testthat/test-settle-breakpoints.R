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
