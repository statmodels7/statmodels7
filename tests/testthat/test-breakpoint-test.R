skip_on_cran()

test_that("the test of a change of slope agrees with Davies' bound by hand", {
  skip_if_not_installed("MASS")
  data(GAGurine, package = "MASS")
  set.seed(1)
  d <- GAGurine[sample(nrow(GAGurine), 40), ]
  d$l <- log(d$GAG)
  f <- statmod(l ~ seg(Age), distrib = distributions7::gaussian1_distrib(),
               data = d)
  r <- statmod_breakpoint_test(f, k = 10)
  expect_s3_class(r, "htest")
  # the score statistic of the change of slope at each position, written out
  # for a gaussian mean: the null fit is least squares on the line, sigma at
  # its maximum likelihood value
  X0 <- cbind(1, d$Age)
  e <- stats::lm.fit(X0, d$l)$residuals
  s2 <- mean(e^2)
  S <- vapply(r$process$psi, function(q) {
    z <- pmax(d$Age - q, 0)
    zr <- stats::lm.fit(X0, z)$residuals
    sum(z * e)^2 / (s2 * sum(zr^2))
  }, 1)
  expect_equal(r$process$score, S, tolerance = 1e-6)
  M <- max(S)
  V <- sum(abs(diff(sqrt(S))))
  p <- stats::pchisq(M, 1, lower.tail = FALSE) + V * exp(-M / 2) / sqrt(2 * pi)
  expect_equal(r$p.value, min(1, p), tolerance = 1e-6)
})

test_that("a sharp step is tested by the bootstrap over every interval", {
  nile <- data.frame(year = 1871:1970, flow = as.numeric(Nile))
  f <- statmod(flow ~ jump(year), distrib = distributions7::gaussian1_distrib(),
               data = nile)
  r <- statmod_breakpoint_test(f, n_boot = 19, seed = 1)
  expect_equal(r$p.value, 1 / 20)
  expect_true(floor(r$estimate) == 1898)
  later <- nile[nile$year > 1900, ]
  f2 <- statmod(flow ~ jump(year), distrib = distributions7::gaussian1_distrib(),
                data = later)
  r2 <- statmod_breakpoint_test(f2, n_boot = 19, seed = 1)
  expect_gt(r2$p.value, 0.2)
})

test_that("the test reaches a smoothed step inside nl()", {
  set.seed(2)
  dd <- data.frame(day = rep(1:30, each = 8),
                   x = rep(seq(0.25, 4, length.out = 8), 30))
  r_true <- exp(log(0.6) + 0.5 * (dd$day > 17))
  dd$y <- rpois(240, 80 * exp(-r_true * dd$x))
  f <- statmod(y ~ 0 + modelterms7::nl(~ a * exp(-r * x),
                 r ~ modelterms7::jump(day, smoothed = numericals7::smooth_probit()),
                 links = list(a = linkfunctions7::log_link(),
                              r = linkfunctions7::log_link()),
                 start = list(a = 50, r = 1)),
               distrib = distributions7::poisson_distrib(
                 link_mu = linkfunctions7::identity_link()),
               data = dd)
  r <- statmod_breakpoint_test(f)
  expect_lt(r$p.value, 1e-10)
  expect_true(r$estimate > 15 && r$estimate < 20)
})

test_that("the test rejects what it cannot read", {
  set.seed(3)
  dd <- data.frame(x = runif(60, 0, 10), g = factor(rep(1:2, 30)))
  dd$y <- dd$x + rnorm(60)
  f <- statmod(y ~ x, distrib = distributions7::gaussian1_distrib(), data = dd)
  expect_error(statmod_breakpoint_test(f), "no seg")
  f2 <- suppressWarnings(statmod(y ~ seg(x, npsi = 2), distrib = distributions7::gaussian1_distrib(),
                data = dd))
  expect_error(statmod_breakpoint_test(f2), "one break-point")
})
