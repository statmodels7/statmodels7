# predict(interval = "group"): the interval and the standard deviation of a
# new group's parameter, theta* = h^-1(eta*), under every kind of link and
# prior. Every reference is an integral the package does not compute.

grp_data <- function(seed = 2, m = 20, ni = 8, kind = c("pois", "bern", "gauss"),
                     b = NULL) {
  kind <- match.arg(kind)
  set.seed(seed)
  g <- factor(rep(seq_len(m), each = ni))
  x <- stats::runif(m * ni)
  if (is.null(b)) b <- stats::rnorm(m, 0, 0.7)
  e <- b[as.integer(g)]
  y <- switch(kind,
              pois = stats::rpois(m * ni, exp(0.5 + x + e)),
              bern = stats::rbinom(m * ni, 1, stats::plogis(0.3 + x + e)),
              gauss = 1 + 2 * x + e + stats::rnorm(m * ni))
  data.frame(g = g, x = x, y = y)
}

# the standard deviation of h^-1(eta) for eta ~ N(m, s^2), by integrate()
sd_normal <- function(hinv, m, s) {
  e1 <- stats::integrate(function(x) hinv(x) * stats::dnorm(x, m, s),
                         m - 12 * s, m + 12 * s, rel.tol = 1e-12)$value
  e2 <- stats::integrate(function(x) hinv(x)^2 * stats::dnorm(x, m, s),
                         m - 12 * s, m + 12 * s, rel.tol = 1e-12)$value
  sqrt(e2 - e1^2)
}

test_that("under a log link the standard deviation is the lognormal one", {
  d <- grp_data()
  f <- statmod(y ~ x + random(~ 1 | g), distributions7::poisson_distrib(), d)
  nd <- data.frame(x = c(0.1, 0.9), g = factor("new"))
  gl <- predict(f, "link:mu", nd, random = "zero", interval = "group")
  gp <- predict(f, "mu", nd, random = "zero", interval = "group")
  z <- stats::qnorm(0.975)
  expect_equal(gp$lower, exp(gl$fit - z * gl$se), tolerance = 1e-12)
  expect_equal(gp$upper, exp(gl$fit + z * gl$se), tolerance = 1e-12)
  ref <- vapply(1:2, function(i) sd_normal(exp, gl$fit[i], gl$se[i]),
                numeric(1))
  expect_equal(gp$se, ref, tolerance = 1e-8)
  # the delta method is a different, smaller number here
  expect_true(all(gp$se > exp(gl$fit) * gl$se * 1.05))
})

test_that("under a logit the standard deviation is a Gauss-Hermite integral", {
  d <- grp_data(kind = "bern", m = 30, b = stats::rnorm(30, 0, 1.2))
  f <- statmod(y ~ x + random(~ 1 | g), distributions7::bernoulli_distrib(), d)
  nd <- data.frame(x = c(0, 0.5, 1), g = factor("new"))
  gl <- predict(f, "link:mu", nd, random = "zero", interval = "group")
  gp <- predict(f, "mu", nd, random = "zero", interval = "group")
  ref <- vapply(1:3, function(i) sd_normal(stats::plogis, gl$fit[i],
                                           gl$se[i]), numeric(1))
  expect_equal(gp$se, ref, tolerance = 1e-7)
  # random = "marginal" moves the fit to the population average and leaves
  # the interval of a new group's value where it is
  gm <- predict(f, "mu", nd, random = "marginal", interval = "group")
  expect_equal(gm[, c("se", "lower", "upper")],
               gp[, c("se", "lower", "upper")])
  expect_true(all(gm$fit < gp$fit | gm$fit > gp$fit))
})

test_that("a Student t prior gives the quantiles of the scale mixture", {
  b <- c(stats::rnorm(17, 0, 0.3), 2.6, -2.2, 3.1)
  d <- grp_data(kind = "gauss", m = 20, ni = 20, seed = 4, b = b)
  tp <- distributions7::fixed(distributions7::student_t1_distrib(), mu = 0)
  f <- statmod(y ~ x + random(~ 1 | g, distrib = tp),
               distributions7::gaussian1_distrib(), d)
  h <- f@hyper$mu[[grep("random", names(f@hyper$mu), value = TRUE)]]
  nd <- data.frame(x = 0.5, g = factor("new"))
  gp <- predict(f, "mu", nd, random = "zero", interval = "group")
  c0 <- predict(f, "mu", nd, random = "zero", se = TRUE)
  G <- function(e) stats::integrate(function(v) {
    w <- exp(v)
    stats::pnorm((e - c0$fit) / sqrt(c0$se^2 + h[["sigma"]]^2 / w)) *
      stats::dgamma(w, h[["nu"]] / 2, h[["nu"]] / 2) * w
  }, -60, 15, rel.tol = 1e-11, subdivisions = 2000L)$value
  expect_equal(c(G(gp$lower), G(gp$upper)), c(0.025, 0.975), tolerance = 1e-6)
  # the variance of a t with nu <= 2 does not exist
  expect_true(h[["nu"]] <= 2)
  expect_true(is.na(gp$se))
})

test_that("another univariate prior gives its own quantiles, seed-free", {
  b <- c(stats::rnorm(17, 0, 0.3), 2.6, -2.2, 3.1)
  d <- grp_data(kind = "gauss", m = 20, ni = 20, seed = 4, b = b)
  lp <- distributions7::fixed(distributions7::logistic_distrib(), mu = 0)
  f <- statmod(y ~ x + random(~ 1 | g, distrib = lp),
               distributions7::gaussian1_distrib(), d)
  s <- f@hyper$mu[[grep("random", names(f@hyper$mu), value = TRUE)]][["sigma"]]
  nd <- data.frame(x = 0.5, g = factor("new"))
  set.seed(1)
  g1 <- predict(f, "mu", nd, random = "zero", interval = "group")
  set.seed(2)
  g2 <- predict(f, "mu", nd, random = "zero", interval = "group")
  expect_identical(g1, g2)
  c0 <- predict(f, "mu", nd, random = "zero", se = TRUE)
  G <- function(e) stats::integrate(function(b) {
    stats::pnorm((e - c0$fit - b) / c0$se) * stats::dlogis(b, 0, s)
  }, -60 * s, 60 * s, rel.tol = 1e-12)$value
  expect_equal(c(G(g1$lower), G(g1$upper)), c(0.025, 0.975), tolerance = 1e-4)
  expect_equal(g1$se, sqrt(c0$se^2 + pi^2 * s^2 / 3), tolerance = 1e-10)
})

test_that("a prior with no moment generating function has no sd under a log link", {
  set.seed(2)
  d <- grp_data(kind = "pois", m = 40, b = 0.8 * stats::rt(40, 3))
  tp <- distributions7::fixed(distributions7::student_t1_distrib(), mu = 0)
  f <- statmod(y ~ x + random(~ 1 | g, distrib = tp),
               distributions7::poisson_distrib(), d)
  gp <- predict(f, "mu", data.frame(x = 0.5, g = factor("new")),
                random = "zero", interval = "group")
  expect_true(is.na(gp$se))
  expect_true(is.finite(gp$lower) && is.finite(gp$upper) &&
                gp$lower < gp$fit && gp$fit < gp$upper)
})

test_that("a predictor whose domain is not the whole line has no sd", {
  sq <- linkfunctions7::sqrt_link()
  expect_true(is.na(group_sd_gaussian(sq, 2, 0.5)))
  expect_identical(link_kind(linkfunctions7::inverse_link()), "restricted")
  expect_identical(link_kind(linkfunctions7::identity_link()), "identity")
  expect_identical(link_kind(linkfunctions7::log_link()), "log")
  expect_identical(link_kind(linkfunctions7::logit_link()), "bounded")
  expect_identical(link_kind(linkfunctions7::softplus_link()), "other")
})

test_that("a mixture quantile is the point where the distribution reaches p", {
  M <- matrix(c(0, 1, 3), 1)
  S <- matrix(c(1, 0.5, 0), 1)
  W <- c(0.5, 0.3, 0.2)
  q <- mixture_quantile(M, S, W, 0.4)
  Fq <- 0.5 * stats::pnorm(q) + 0.3 * stats::pnorm((q - 1) / 0.5)
  expect_equal(Fq, 0.4, tolerance = 1e-10)
  # a point mass: the smallest value with F >= p
  expect_equal(mixture_quantile(matrix(c(0, 3), 1), matrix(0, 1, 2), c(0.5, 0.5),
                                0.75), 3, tolerance = 1e-10)
})
