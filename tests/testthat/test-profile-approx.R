# profile_continuation() and profile_taylor(): the approximate methods of
# statmod_profile(). The references are profile_grid() (the refit at every
# point) and, for the quadratic expansion, the Wald profile of glm().

poisson_approx <- function() {
  set.seed(11)
  dd <- data.frame(x1 = rnorm(150), x2 = runif(150))
  dd$y <- rpois(150, exp(0.3 + 0.5 * dd$x1 - 0.4 * dd$x2))
  list(dd = dd,
       fit = statmod(y ~ x1 + x2, distributions7::poisson_distrib(), dd),
       glm = stats::glm(y ~ x1 + x2, stats::poisson, dd))
}

err <- function(a, b) max(abs(a$profile_loglik - b$profile_loglik))

test_that("continuation to convergence gives the rows of the grid", {
  p <- poisson_approx()
  v <- profile_range(width = 3, n = 7)
  g <- statmod_profile(p$fit, "mu:x1", v)
  for (o in 1:2) {
    c1 <- statmod_profile(p$fit, "mu:x1", v, method = profile_continuation(order = o))
    expect_equal(c1$profile_loglik, g$profile_loglik, tolerance = 1e-9)
    expect_equal(c1$profile_pen_loglik, g$profile_pen_loglik, tolerance = 1e-9)
    expect_true(all(c1$converged))
  }
  expect_identical(attr(c1, "method")@name, "continuation")
  # two held coefficients
  v2 <- list("mu:x1" = c(0.3, 0.5, 0.7), "mu:x2" = c(-0.6, -0.3))
  g2 <- statmod_profile(p$fit, c("mu:x1", "mu:x2"), v2)
  c2 <- statmod_profile(p$fit, c("mu:x1", "mu:x2"), v2, method = profile_continuation())
  expect_equal(c2$profile_loglik, g2$profile_loglik, tolerance = 1e-9)
})

test_that("a truncated corrector is a lower bound of the grid", {
  p <- poisson_approx()
  v <- profile_range(width = 3, n = 7)
  g <- statmod_profile(p$fit, "mu:x1", v)
  c1 <- statmod_profile(p$fit, "mu:x1", v,
                        method = profile_continuation(order = 1, maxit = 1))
  expect_true(all(c1$profile_pen_loglik <= g$profile_pen_loglik + 1e-9))
  expect_lt(err(c1, g), 1e-6)
})

test_that("the expansion approaches the grid as the order grows", {
  p <- poisson_approx()
  v <- profile_range(width = 2, n = 5)
  g <- statmod_profile(p$fit, "mu:x1", v)
  e <- vapply(2:5, function(o) {
    t <- statmod_profile(p$fit, "mu:x1", v, method = profile_taylor(order = o))
    expect_true(all(t$profile_pen_loglik <= g$profile_pen_loglik + 1e-9))
    expect_true(all(is.na(t$converged)))
    expect_true(all(is.na(t$outer_criterion)))
    err(t, g)
  }, numeric(1))
  expect_true(all(diff(e) < 0))
  expect_lt(e[4], 1e-8)
})

test_that("order 2 agrees with the quadratic profile near the estimate", {
  p <- poisson_approx()
  se <- sqrt(stats::vcov(p$glm)["x1", "x1"])
  b0 <- unname(stats::coef(p$glm)["x1"])
  b <- b0 + se * c(-0.3, -0.1, 0.1, 0.3)
  t <- statmod_profile(p$fit, "mu:x1", list(b), method = profile_taylor(order = 2))
  ref <- as.numeric(stats::logLik(p$glm)) - (b - b0)^2 / (2 * se^2)
  # the first neglected term is cubic in the step
  expect_equal(t$profile_loglik, ref, tolerance = 5e-4)
  expect_gt(max(abs(t$profile_loglik - ref)), 0)
})

test_that("a penalized smooth is followed by both methods", {
  set.seed(3)
  dd <- data.frame(z = runif(250), x = rnorm(250))
  dd$y <- rpois(250, exp(1 + sin(2 * pi * dd$z) + 0.4 * dd$x))
  fit <- statmod(y ~ s(z) + x, distributions7::poisson_distrib(), dd)
  v <- profile_range(width = 2, n = 5)
  g <- statmod_profile(fit, "mu:x", v)
  expect_equal(statmod_profile(fit, "mu:x", v, method = profile_continuation())$profile_loglik,
               g$profile_loglik, tolerance = 1e-9)
  t <- statmod_profile(fit, "mu:x", v, method = profile_taylor(order = 4))
  expect_lt(err(t, g), 1e-5)
  expect_true(all(t$profile_pen_loglik <= g$profile_pen_loglik + 1e-9))
  # the penalty column is that of the exact point
  expect_equal(t$profile_penalty, -(-t$profile_pen_loglik + t$profile_loglik))
})

test_that("a nonlinear term, whose block moves, is followed", {
  set.seed(2)
  dd <- data.frame(x = rep(seq(0.25, 4, length.out = 8), 25),
                   z = rep(runif(25), each = 8))
  dd$y <- 50 * exp(-exp(log(0.6) + 0.5 * dd$z) * dd$x) + stats::rnorm(200)
  fit <- statmod(
    y ~ 0 + modelterms7::nl(~ a * exp(-r * x), r ~ 1 + z,
      links = list(a = linkfunctions7::log_link(),
                   r = linkfunctions7::log_link()),
      start = list(a = 40, r = 1)),
    distrib = distributions7::gaussian1_distrib(
      link_mu = linkfunctions7::identity_link()),
    data = dd, outer_criterion = NULL)
  v <- profile_range(width = 2, n = 5)
  g <- statmod_profile(fit, "mu:nl.r.z", v)
  expect_equal(statmod_profile(fit, "mu:nl.r.z", v,
                               method = profile_continuation(order = 2))$profile_loglik,
               g$profile_loglik, tolerance = 1e-9)
  e <- vapply(c(2, 3, 5), function(o) {
    t <- statmod_profile(fit, "mu:nl.r.z", v, method = profile_taylor(order = o))
    expect_true(all(t$profile_pen_loglik <= g$profile_pen_loglik + 1e-9))
    err(t, g)
  }, numeric(1))
  expect_true(all(diff(e) < 0))
})

test_that("coordinates that reml estimates follow the nested problem", {
  set.seed(1)
  dd <- data.frame(x = runif(150), z = runif(150))
  dd$y <- 1 + 2 * dd$x + stats::rnorm(150, sd = exp(-0.9 - 0.2 * dd$z))
  fit <- statmod(y ~ x | sigma ~ z, distributions7::gaussian1_distrib(), dd)
  expect_identical(fit@methods$outer@kind, "reml")
  v <- profile_range(width = 2, n = 5)
  g <- statmod_profile(fit, "mu:x", v)
  for (o in 1:2) {
    c1 <- statmod_profile(fit, "mu:x", v, method = profile_continuation(order = o))
    expect_equal(c1$profile_loglik, g$profile_loglik, tolerance = 1e-6)
    expect_equal(c1$outer_criterion, g$outer_criterion, tolerance = 1e-6)
  }
  e <- vapply(2:5, function(o) {
    err(statmod_profile(fit, "mu:x", v, method = profile_taylor(order = o)), g)
  }, numeric(1))
  expect_true(all(diff(e) < 0))
  expect_lt(e[4], 1e-4)
})

test_that("the exact derivative of the nested problem agrees with the refits", {
  set.seed(1)
  dd <- data.frame(x = runif(150), z = runif(150))
  dd$y <- 1 + 2 * dd$x + stats::rnorm(150, sd = exp(-0.9 - 0.2 * dd$z))
  fit <- statmod(y ~ x | sigma ~ z, distributions7::gaussian1_distrib(), dd)
  ctx <- profile_context(fit, "mu:x", "held")
  base <- profile_point(ctx, ctx$estimate)
  base$val <- ctx$estimate
  N <- profile_nested_state(ctx, base)
  se <- profile_setup(ctx, base$par)$se
  # the central difference of the refits at two values, the reference that
  # shares no arithmetic with the Hessian of the criterion
  eps <- 0.1 * se
  zp <- profile_point(ctx, ctx$estimate + eps)
  zm <- profile_point(ctx, ctx$estimate - eps)
  fd <- (zp$par[N$gam] - zm$par[N$gam]) / (2 * eps)
  exact <- as.numeric(-solve(N$Vuu, N$Vuk))
  expect_equal(exact, fd, tolerance = 1e-3)
})

test_that("a halved step reaches the same point", {
  p <- poisson_approx()
  v <- profile_range(width = 3, n = 3)
  g <- statmod_profile(p$fit, "mu:x1", v)
  calls <- new.env()
  calls$n <- 0L
  calls$failed <- FALSE
  real <- profile_point
  testthat::local_mocked_bindings(
    profile_point = function(ctx, val, start = NULL, maxit = Inf) {
      calls$n <- calls$n + 1L
      if (!calls$failed && !is.null(start) && val > ctx$estimate) {
        calls$failed <- TRUE
        stop("a corrector that fails")
      }
      real(ctx, val, start, maxit)
    })
  c1 <- statmod_profile(p$fit, "mu:x1", v, method = profile_continuation())
  expect_true(calls$failed)
  expect_identical(nrow(c1), 3L)
  expect_equal(c1$profile_loglik, g$profile_loglik, tolerance = 1e-9)
  # the rows are the requested values only
  expect_equal(c1[["mu:x1"]], g[["mu:x1"]])
})

test_that("models without smooth coefficients are rejected with a pointer to the grid", {
  p <- poisson_approx()
  set.seed(7)
  dd <- data.frame(x1 = rnorm(80), x2 = rnorm(80), x3 = rnorm(80))
  dd$y <- 1 + dd$x1 + stats::rnorm(80)
  fl <- statmod(y ~ x3 + lasso(~ 0 + x1 + x2, lambda = 0.5),
                distributions7::gaussian1_distrib(), dd)
  expect_error(statmod_profile(fl, "mu:x3", method = profile_taylor()),
               "profile_grid")
  set.seed(3)
  n <- 400
  ds <- data.frame(w = runif(n, 0, 10))
  ds$y <- 0.8 * pmax(ds$w - 5, 0) + stats::rnorm(n, sd = 0.5)
  fs <- statmod(y ~ seg(w), distributions7::gaussian1_distrib(), ds)
  expect_error(statmod_profile(fs, "mu:(Intercept)", method = profile_continuation()),
               "without smoothing")
})

test_that("a smoothed break-point term is accepted", {
  set.seed(3)
  n <- 400
  ds <- data.frame(w = runif(n, 0, 10))
  ds$y <- 0.8 * pmax(ds$w - 5, 0) + stats::rnorm(n, sd = 0.5)
  fs <- statmod(y ~ seg(w, smoothed = numericals7::smooth_probit()),
                distributions7::gaussian1_distrib(), ds, outer_criterion = NULL)
  v <- profile_range(width = 0.3, n = 3)
  g <- statmod_profile(fs, "mu:(Intercept)", v)
  t <- statmod_profile(fs, "mu:(Intercept)", v, method = profile_taylor(order = 3))
  expect_lt(err(t, g), 1e-4)
  c1 <- statmod_profile(fs, "mu:(Intercept)", v, method = profile_continuation())
  expect_equal(c1$profile_loglik, g$profile_loglik, tolerance = 1e-8)
})

test_that("the constructors validate their arguments", {
  expect_error(profile_taylor(order = 1), "2 to 5")
  expect_error(profile_taylor(order = 6), "2 to 5")
  expect_error(profile_continuation(order = 3), "1 or 2")
  expect_error(profile_continuation(maxit = 0), "maxit")
  expect_s3_class(S7::S7_class(profile_taylor()), "S7_class")
  expect_true(S7::S7_inherits(profile_continuation(), profile_method))
})

test_that("the Hessian of a moving block is the difference of the gradient", {
  set.seed(2)
  dd <- data.frame(x = rep(seq(0.25, 4, length.out = 8), 25),
                   z = rep(runif(25), each = 8))
  dd$y <- 50 * exp(-exp(log(0.6) + 0.5 * dd$z) * dd$x) + stats::rnorm(200)
  fit <- statmod(
    y ~ 0 + modelterms7::nl(~ a * exp(-r * x), r ~ 1 + z,
      links = list(a = linkfunctions7::log_link(),
                   r = linkfunctions7::log_link()),
      start = list(a = 40, r = 1)),
    distrib = distributions7::gaussian1_distrib(
      link_mu = linkfunctions7::identity_link()),
    data = dd, outer_criterion = NULL)
  ctx <- profile_context(fit, "mu:nl.r.z", "held")
  S <- profile_setup(ctx, ctx$beta_fit)
  b <- ctx$beta_fit
  # the reference: central differences of the analytic gradient
  fd <- vapply(seq_along(b), function(j) {
    h <- 1e-5 * max(1, abs(b[j]))
    e <- numeric(length(b))
    e[j] <- h
    (S$obj$gr(b + e) - S$obj$gr(b - e)) / (2 * h)
  }, numeric(length(b)))
  fd <- (fd + t(fd)) / 2
  expect_lt(max(abs(S$H - fd)) / max(abs(fd)), 1e-7)
})
