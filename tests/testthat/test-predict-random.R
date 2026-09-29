# predict(random =): a group the fit never saw, read as the typical group
# ("zero") or averaged over the prior ("marginal"). Every reference below
# shares no arithmetic with the Gauss-Hermite grid: a closed form, or
# stats::integrate().

set.seed(3)
m <- 20; ni <- 8
g <- factor(rep(sprintf("g%02d", seq_len(m)), each = ni))
x <- stats::rnorm(m * ni)
b0 <- stats::rnorm(m, 0, 0.7)
new <- data.frame(g = c("zz", "g01"), x = c(0.5, -1))
sd_of <- function(fit) {
  h <- fit@hyper$mu
  h[[grep("random", names(h), fixed = TRUE)]][["sigma"]]
}

test_that("a new level names the argument that predicts it", {
  dg <- data.frame(g, x, y = 1 + 0.5 * x + b0[g] + stats::rnorm(m * ni))
  fg <- statmod(y ~ x + random(~ 1 | g), distributions7::gaussian1_distrib(),
                dg)
  expect_error(predict(fg, "mu", new), "random = \"zero\"")
})

test_that("identity link: zero and marginal agree, the variance adds sigma_b^2", {
  dg <- data.frame(g, x, y = 1 + 0.5 * x + b0[g] + stats::rnorm(m * ni, sd = 0.5))
  fg <- statmod(y ~ x + random(~ 1 | g), distributions7::gaussian1_distrib(),
                dg)
  z <- predict(fg, "mu", new, random = "zero")
  expect_equal(predict(fg, "mu", new, random = "marginal"), z,
               tolerance = 1e-12)
  expect_equal(predict(fg, "variance", new, random = "marginal") -
                 predict(fg, "variance", new, random = "zero"),
               rep(sd_of(fg)^2, 2), tolerance = 1e-12)
  # "zero" sets the effect aside at every row, a group the fit saw included,
  # and at the fitting rows it is the fixed part
  expect_equal(z[[2L]], sum(fg@coefficients$mu[1:2] * c(1, -1)),
               tolerance = 1e-12)
  f0 <- predict(fg, "link:mu", random = "zero")
  expect_equal(as.numeric(f0),
               as.numeric(cbind(1, dg$x) %*% fg@coefficients$mu[1:2]),
               tolerance = 1e-12)
})

test_that("log link: the marginal mean is exp(eta0 + sigma_b^2 / 2)", {
  dp <- data.frame(g, x, y = stats::rpois(m * ni, exp(0.5 + 0.3 * x + b0[g])))
  fp <- statmod(y ~ x + random(~ 1 | g), distributions7::poisson_distrib(), dp)
  e0 <- predict(fp, "link:mu", new, random = "zero")
  expect_equal(predict(fp, "mu", new, random = "marginal"),
               exp(e0 + sd_of(fp)^2 / 2), tolerance = 1e-12)
  s <- predict(fp, "mu", new, random = "marginal", se = TRUE)
  expect_true(all(s$lower <= s$fit & s$fit <= s$upper))
  expect_true(all(s$se > 0))
})

test_that("logit link: the marginal mean is the integral, pulled to one half", {
  db <- data.frame(g, x, y = stats::rbinom(m * ni, 1,
                                           stats::plogis(-0.2 + 0.8 * x +
                                                           1.2 * b0[g])))
  fb <- statmod(y ~ x + random(~ 1 | g), distributions7::bernoulli_distrib(),
                db)
  e0 <- predict(fb, "link:mu", new, random = "zero")
  s <- sd_of(fb)
  skip_if(s < 0.1, "the prior scale ran to zero on this sample")
  ref <- vapply(e0, function(e) stats::integrate(function(b)
    stats::plogis(e + b) * stats::dnorm(b, 0, s), -Inf, Inf,
    rel.tol = 1e-12)$value, numeric(1))
  mb <- predict(fb, "mu", new, random = "marginal")
  expect_equal(mb, ref, tolerance = 1e-9)
  expect_true(all(abs(mb - 0.5) < abs(stats::plogis(e0) - 0.5)))
})

test_that("a correlated random slope integrates exp(eta0 + z'Sz/2)", {
  b1 <- stats::rnorm(m, 0, 0.4)
  dq <- data.frame(g, x, y = stats::rpois(m * ni, exp(0.3 + 0.2 * x + b0[g] +
                                                       b1[g] * x)))
  fq <- statmod(y ~ x + random(~ 1 + x | g), distributions7::poisson_distrib(),
                dq)
  key <- grep("random", names(fq@spec@terms$mu), value = TRUE)
  S <- solve(as.matrix(penalties7::penalty_hessian(
    modelterms7::term_penalty(fq@spec@terms$mu[[key]]),
    rep(0, 2L * m), as.list(fq@hyper$mu[[key]])))[1:2, 1:2])
  e0 <- predict(fq, "link:mu", new, random = "zero")
  Z <- cbind(1, new$x)
  expect_equal(predict(fq, "mu", new, random = "marginal"),
               exp(e0 + rowSums((Z %*% S) * Z) / 2), tolerance = 1e-10)
})

test_that("modes are chosen term by term and unknown names are refused", {
  dg <- data.frame(g, x, y = 1 + b0[g] + stats::rnorm(m * ni))
  fg <- statmod(y ~ x + random(~ 1 | g), distributions7::gaussian1_distrib(),
                dg)
  expect_equal(predict(fg, "mu", new, random = c("random(~1|g)" = "zero")),
               predict(fg, "mu", new, random = "zero"))
  expect_error(predict(fg, "mu", new, random = c("random(~1|h)" = "zero")),
               "not a random-effect term")
  expect_error(predict(fg, "mu", new, random = "typical"), "must be")
  expect_error(predict(fg, "skewness", new, random = "marginal"),
               "no marginal prediction")
})

test_that("a Student t prior is integrated only under a bounded link", {
  tprior <- distributions7::fixed(distributions7::student_t1_distrib(),
                                  mu = 0)
  dp <- data.frame(g, x, y = stats::rpois(m * ni, exp(0.5 + b0[g])))
  fp <- statmod(y ~ x + random(~ 1 | g, distrib = tprior),
                distributions7::poisson_distrib(), dp)
  expect_error(predict(fp, "mu", new, random = "marginal"), "not bounded")
  db <- data.frame(g, x, y = stats::rbinom(m * ni, 1, stats::plogis(b0[g])))
  fb <- statmod(y ~ x + random(~ 1 | g, distrib = tprior),
                distributions7::bernoulli_distrib(), db)
  set.seed(99); a <- stats::runif(1)
  t1 <- predict(fb, "mu", new, random = "marginal")
  set.seed(99); t2 <- predict(fb, "mu", new, random = "marginal")
  a2 <- stats::runif(1)
  # reproducible, and the caller's random stream is left where it was
  expect_identical(t1, t2)
  expect_identical(a, a2)
  expect_true(all(t1 > 0 & t1 < 1))
})
