# predict(random = "marginal") under links other than the log and the logit,
# and under one univariate prior that is not Gaussian. The references are
# closed forms or integrate(), which share no code with the package.

ml_data <- function(seed = 2, m = 40, ni = 8, b = NULL) {
  set.seed(seed)
  g <- factor(rep(seq_len(m), each = ni))
  x <- stats::runif(m * ni)
  if (is.null(b)) b <- stats::rnorm(m, 0, 0.8)
  data.frame(g = g, x = x, b = b[as.integer(g)])
}
ml_new <- data.frame(x = c(0.2, 0.8), g = factor("new"))

marg_by_integrate <- function(f, hinv, dens) {
  e0 <- predict(f, "link:mu", ml_new, random = "zero")
  vapply(e0, function(e) stats::integrate(function(u) hinv(e + u) * dens(u),
                                          -Inf, Inf, rel.tol = 1e-12)$value,
         numeric(1))
}

test_that("the probit and the square root have their closed forms", {
  d <- ml_data()
  d$y <- stats::rbinom(nrow(d), 1, stats::pnorm(-0.3 + d$x + d$b))
  f <- statmod(y ~ x + random(~ 1 | g),
               distributions7::bernoulli_distrib(
                 link_mu = linkfunctions7::probit_link()), d)
  tau <- hyper(f)$estimate
  e0 <- predict(f, "link:mu", ml_new, random = "zero")
  expect_equal(predict(f, "mu", ml_new, random = "marginal"),
               stats::pnorm(e0 / sqrt(1 + tau^2)), tolerance = 1e-12)
  d$y <- stats::rpois(nrow(d), (2 + d$x + 0.3 * d$b)^2)
  f <- statmod(y ~ x + random(~ 1 | g),
               distributions7::poisson_distrib(
                 link_mu = linkfunctions7::sqrt_link()), d,
               outer_criterion = ml())
  tau <- hyper(f)$estimate
  e0 <- predict(f, "link:mu", ml_new, random = "zero")
  expect_equal(predict(f, "mu", ml_new, random = "marginal"), e0^2 + tau^2,
               tolerance = 1e-12)
})

test_that("a cloglog average is the integral over the prior", {
  d <- ml_data()
  d$y <- stats::rbinom(nrow(d), 1, 1 - exp(-exp(-0.8 + d$x + d$b)))
  f <- statmod(y ~ x + random(~ 1 | g),
               distributions7::bernoulli_distrib(
                 link_mu = linkfunctions7::cloglog_link()), d)
  tau <- hyper(f)$estimate
  ref <- marg_by_integrate(f, function(e) 1 - exp(-exp(e)),
                           function(u) stats::dnorm(u, 0, tau))
  expect_equal(predict(f, "mu", ml_new, random = "marginal"), ref,
               tolerance = 1e-12)
})

test_that("one univariate prior that is not Gaussian is integrated exactly", {
  set.seed(4)
  b <- c(stats::rnorm(37, 0, 0.4), 4, -3.5, 3.8)
  d <- ml_data(seed = 4, m = 40, ni = 12, b = b)
  d$y <- stats::rbinom(nrow(d), 1, stats::plogis(0.2 + d$x + d$b))
  tp <- distributions7::fixed(distributions7::student_t1_distrib(), mu = 0)
  f <- statmod(y ~ x + random(~ 1 | g, distrib = tp),
               distributions7::bernoulli_distrib(), d)
  h <- f@hyper$mu[[grep("random", names(f@hyper$mu), value = TRUE)]]
  dens <- function(u) stats::dt(u / h[["sigma"]], h[["nu"]]) / h[["sigma"]]
  ref <- marg_by_integrate(f, stats::plogis, dens)
  set.seed(1)
  p1 <- predict(f, "mu", ml_new, random = "marginal")
  set.seed(2)
  p2 <- predict(f, "mu", ml_new, random = "marginal")
  expect_identical(p1, p2)
  expect_equal(p1, ref, tolerance = 1e-7)
  # the interval's ends are the average at the ends of the predictor's
  # interval, and the standard error the averaged derivative times its se
  ps <- predict(f, "mu", ml_new, random = "marginal", se = TRUE)
  pl <- predict(f, "link:mu", ml_new, random = "zero", se = TRUE)
  at <- function(e) vapply(e, function(ei) stats::integrate(function(u)
    stats::plogis(ei + u) * dens(u), -Inf, Inf, rel.tol = 1e-12)$value, 0)
  expect_equal(ps$lower, at(pl$lower), tolerance = 1e-7)
  expect_equal(ps$upper, at(pl$upper), tolerance = 1e-7)
  dd <- vapply(pl$fit, function(ei) stats::integrate(function(u)
    stats::dlogis(ei + u) * dens(u), -Inf, Inf, rel.tol = 1e-12)$value, 0)
  expect_equal(ps$se, dd * pl$se, tolerance = 1e-7)
  # a logistic prior, and a moment of the response
  lp <- distributions7::fixed(distributions7::logistic_distrib(), mu = 0)
  fl <- statmod(y ~ x + random(~ 1 | g, distrib = lp),
                distributions7::bernoulli_distrib(), d)
  s <- fl@hyper$mu[[grep("random", names(fl@hyper$mu), value = TRUE)]][["sigma"]]
  ref <- marg_by_integrate(fl, stats::plogis,
                           function(u) stats::dlogis(u, 0, s))
  expect_equal(predict(fl, "mean", ml_new, random = "marginal"), ref,
               tolerance = 1e-7)
})

test_that("a predictor that leaves its domain gives NA, not a number", {
  inv <- linkfunctions7::inverse_link()
  nodes <- list(b = list(matrix(c(-1, 1), 2, 1)), w = c(0.5, 0.5))
  Zs <- list(matrix(1, 2, 1))
  # s = 1: at eta0 = 0.1 the predictor is negative with probability 0.46,
  # at eta0 = 8 with probability 6e-16
  expect_warning(bad <- marginal_out_of_domain(inv, c(0.1, 8), Zs, nodes, 1L),
                 "leaves the domain")
  expect_identical(bad, c(TRUE, FALSE))
  expect_identical(marginal_out_of_domain(linkfunctions7::logit_link(),
                                          c(0.1, 8), Zs, nodes, 1L),
                   c(FALSE, FALSE))
})

test_that("predict() reports NA where the predictor leaves its domain", {
  set.seed(3)
  g <- factor(rep(1:30, each = 6))
  d <- data.frame(g = g, x = stats::runif(180))
  d$y <- stats::rgamma(180, shape = 5,
                       rate = 5 * (1.5 + d$x + 0.15 * stats::rnorm(30)[g]))
  f <- statmod(y ~ x + random(~ 1 | g),
               distributions7::gamma1_distrib(
                 link_mu = linkfunctions7::inverse_link()), d)
  nd <- data.frame(x = c(0, 1), g = factor("new"))
  e0 <- predict(f, "link:mu", nd, random = "zero")
  # the fitted prior is narrow, so the average exists; a prior as wide as
  # the predictor itself makes it leave (0, Inf) with probability 0.16
  expect_true(all(is.finite(predict(f, "mu", nd, random = "marginal"))))
  rp <- random_prior
  local_mocked_bindings(random_prior = function(...) {
    pr <- rp(...)
    pr$chol <- pr$chol / pr$chol[1, 1] * mean(e0)
    pr
  })
  expect_warning(pm <- predict(f, "mu", nd, random = "marginal"),
                 "leaves the domain")
  expect_true(all(is.na(pm)))
})
