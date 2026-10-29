# statmod_profile(): the profile likelihood of coefficients on a grid.
# References that share no arithmetic with the function: stats::profile.glm
# for one coefficient, glm() with an offset for two, and statmod_restrict()
# for a single point.

poisson_fit <- function() {
  set.seed(11)
  dd <- data.frame(x1 = rnorm(60), x2 = runif(60))
  dd$y <- rpois(60, exp(0.3 + 0.5 * dd$x1 - 0.4 * dd$x2))
  list(dd = dd,
       fit = statmod(y ~ x1 + x2, distributions7::poisson_distrib(), dd),
       glm = stats::glm(y ~ x1 + x2, stats::poisson, dd))
}

test_that("one coefficient reproduces stats::profile.glm", {
  p <- poisson_fit()
  pg <- stats::profile(p$glm, which = "x1", alpha = 0.001)
  b <- pg$x1$par.vals[, "x1"]
  ref <- as.numeric(stats::logLik(p$glm)) - pg$x1$z^2 / 2
  prof <- statmod_profile(p$fit, which = "mu:x1", values = list(b))
  expect_equal(prof$profile_loglik, ref, tolerance = 1e-9)
  expect_equal(prof$profile_pen_loglik, prof$profile_loglik)
  expect_equal(prof$profile_penalty, rep(0, length(b)))
  expect_true(all(is.na(prof$outer_criterion)))
  expect_true(all(prof$converged))
})

test_that("two coefficients reproduce glm() refits with an offset", {
  p <- poisson_fit()
  a <- c(0.3, 0.5, 0.7)
  b <- c(-0.6, -0.4)
  prof <- statmod_profile(p$fit, which = c("mu:x1", "mu:x2"),
                          values = list("mu:x2" = b, "mu:x1" = a))
  expect_identical(names(prof)[1:2], c("mu:x1", "mu:x2"))
  expect_identical(nrow(prof), 6L)
  # the first coefficient varies fastest
  expect_equal(prof[["mu:x1"]], rep(a, 2))
  expect_equal(prof[["mu:x2"]], rep(b, each = 3))
  ref <- vapply(seq_len(nrow(prof)), function(i) {
    g <- stats::glm(y ~ 1 + offset(prof[i, 1] * x1 + prof[i, 2] * x2),
                    stats::poisson, p$dd)
    as.numeric(stats::logLik(g))
  }, numeric(1))
  expect_equal(prof$profile_loglik, ref, tolerance = 1e-9)
})

test_that("a single point is the restricted fit", {
  p <- poisson_fit()
  r <- statmod_restrict(p$fit, param = "mu", coefname = "x1", value = 0.3)
  prof <- statmod_profile(p$fit, which = "mu:x1", values = list(0.3))
  expect_equal(prof$profile_loglik, r$loglik, tolerance = 1e-10)
  expect_equal(prof$profile_pen_loglik, -r$objective, tolerance = 1e-10)
})

test_that("the three starts agree, and the maximum is the fit", {
  p <- poisson_fit()
  v <- profile_range(width = 2, n = 5)
  a <- statmod_profile(p$fit, "mu:x1", v, method = profile_grid("neighbor"))
  b <- statmod_profile(p$fit, "mu:x1", v, method = profile_grid("fit"))
  d <- statmod_profile(p$fit, "mu:x1", v, method = profile_grid("both"))
  expect_equal(a$profile_loglik, b$profile_loglik, tolerance = 1e-9)
  expect_equal(a$profile_loglik, d$profile_loglik, tolerance = 1e-9)
  mx <- attr(a, "maximum")
  expect_identical(nrow(mx), 1L)
  expect_identical(names(mx), names(a))
  expect_equal(mx$profile_loglik, p$fit@loglik)
  expect_equal(mx$profile_pen_loglik, -p$fit@objective)
  # n odd: the middle row is the estimate, and it is the largest
  expect_equal(a[["mu:x1"]][3], mx[["mu:x1"]])
  expect_equal(a$profile_loglik[3], mx$profile_loglik, tolerance = 1e-8)
  expect_identical(which.max(a$profile_loglik), 3L)
  expect_true(S7::S7_inherits(attr(a, "method"), profile_method))
  expect_identical(attr(a, "hyper"), "held")
})

test_that("a model without hyperparameters is the same held or refitted", {
  p <- poisson_fit()
  v <- profile_range(width = 1, n = 3)
  a <- statmod_profile(p$fit, "mu:x2", v, hyper = "held")
  b <- statmod_profile(p$fit, "mu:x2", v, hyper = "refit")
  expect_equal(a$profile_loglik, b$profile_loglik)
  expect_equal(a$profile_pen_loglik, b$profile_pen_loglik)
})

test_that("a smooth: the penalty column makes the sum the objective", {
  set.seed(4)
  dd <- data.frame(x = runif(120), z = rnorm(120))
  dd$y <- rpois(120, exp(0.5 + sin(2 * pi * dd$x) + 0.3 * dd$z))
  fit <- statmod(y ~ s(x) + z, distributions7::poisson_distrib(), dd)
  prof <- statmod_profile(fit, "mu:z", profile_range(width = 1, n = 3))
  expect_true(all(prof$profile_penalty != 0))
  expect_equal(prof$profile_pen_loglik,
               prof$profile_loglik + prof$profile_penalty)
  # the hyperparameters are held, so the point is the restricted fit
  r <- statmod_restrict(fit, param = "mu", coefname = "z",
                        value = prof[["mu:z"]][1])
  expect_equal(prof$profile_pen_loglik[1], -r$objective, tolerance = 1e-9)
  # and the unrestricted maximum has the fit's own numbers
  mx <- attr(prof, "maximum")
  expect_equal(mx$profile_pen_loglik, -fit@objective)
  expect_equal(mx$profile_penalty, -(fit@objective + fit@loglik))
  # with the hyperparameters refitted the outer search runs at every point
  pr <- statmod_profile(fit, "mu:z", profile_range(width = 1, n = 3),
                        hyper = "refit")
  expect_true(all(is.finite(pr$outer_criterion)))
  expect_true(all(is.finite(pr$profile_pen_loglik)))
})

test_that("coefficients the outer criterion estimates are estimated by it", {
  set.seed(1)
  dd <- data.frame(x = runif(80), z = runif(80))
  dd$y <- 1 + 2 * dd$x + stats::rnorm(80, sd = exp(-0.9 - 0.2 * dd$z))
  fit <- statmod(y ~ x | sigma ~ z, distributions7::gaussian1_distrib(), dd)
  expect_identical(fit@methods$outer@kind, "reml")
  v <- profile_range(width = 2, n = 5)
  prof <- statmod_profile(fit, "mu:x", v)
  expect_true(all(is.finite(prof$outer_criterion)))
  # the criterion, not the joint mode: the restricted fit carries sigma's
  # coefficients to the maximum of the likelihood, and the two differ
  r <- statmod_restrict(fit, param = "mu", coefname = "x",
                        value = prof[["mu:x"]][1])
  expect_gt(abs(prof$profile_loglik[1] - r$loglik), 1e-4)
  # a coefficient of sigma is held and the criterion estimates the rest
  ps <- statmod_profile(fit, "sigma:z", profile_range(width = 2, n = 5))
  expect_true(all(is.finite(ps$profile_loglik)))
  expect_equal(which.max(ps$profile_loglik), 3L)
  # the start does not change the answer
  b <- statmod_profile(fit, "mu:x", v, method = profile_grid("fit"))
  expect_equal(prof$profile_loglik, b$profile_loglik, tolerance = 1e-5)
  mx <- attr(prof, "maximum")
  expect_equal(mx$outer_criterion, fit@criterion[[1L]])
})

test_that("the label of a coefficient of a nonlinear term is found", {
  set.seed(2)
  dd <- data.frame(x = rep(seq(0.25, 4, length.out = 8), 25),
                   z = rep(runif(25), each = 8))
  r <- exp(log(0.6) + 0.5 * dd$z)
  dd$y <- 50 * exp(-r * dd$x) + stats::rnorm(200, sd = 1)
  fit <- statmod(
    y ~ 0 + modelterms7::nl(~ a * exp(-r * x), r ~ 1 + z,
      links = list(a = linkfunctions7::log_link(),
                   r = linkfunctions7::log_link()),
      start = list(a = 40, r = 1)),
    distrib = distributions7::gaussian1_distrib(
      link_mu = linkfunctions7::identity_link()), data = dd)
  prof <- statmod_profile(fit, "mu:nl.r.z", profile_range(width = 1, n = 3))
  expect_identical(names(prof)[1L], "mu:nl.r.z")
  expect_identical(which.max(prof$profile_loglik), 2L)
})

test_that("values: a range, a list with names, a default, and the errors", {
  p <- poisson_fit()
  prof <- statmod_profile(
    p$fit, c("mu:x1", "mu:x2"),
    values = list("mu:x1" = c(0.4, 0.6), profile_range(width = 1, n = 3)))
  expect_identical(nrow(prof), 6L)
  expect_equal(prof[["mu:x1"]], rep(c(0.4, 0.6), 3))
  # a repeated value is dropped and the values are sorted
  q <- statmod_profile(p$fit, "mu:x1", list(c(0.6, 0.4, 0.6)))
  expect_equal(q[["mu:x1"]], c(0.4, 0.6))
  expect_error(statmod_profile(p$fit, "mu:x1", list(1, 2)), "more than one")
  expect_error(statmod_profile(p$fit, "mu:x1", list("mu:x2" = 1)),
               "not in 'which'")
  expect_error(statmod_profile(p$fit, "mu:x1", c(1, 2)), "wrap a numeric")
  expect_error(statmod_profile(p$fit, "mu:nope"), "Not a coefficient")
  expect_error(statmod_profile(p$fit, c("mu:x1", "mu:x1")), "repeat")
  expect_error(statmod_profile(p$fit, "mu:x1", cores = 0), "cores")
  expect_error(profile_range(width = -1), "positive")
  expect_error(profile_range(n = 1), "at least 2")
  expect_error(profile_grid("middle"))
  expect_error(profile_method(name = "x"))
})

test_that("a coefficient that cannot be held is rejected", {
  set.seed(5)
  dd <- data.frame(x1 = rnorm(80), x2 = rnorm(80))
  dd$y <- 1 + dd$x1 + stats::rnorm(80)
  fl <- statmod(y ~ lasso(~ 0 + x1 + x2, lambda = 0.5),
                distributions7::gaussian1_distrib(), dd)
  lab <- rownames(coef_labels(fl@spec, statmod_design(fl@spec)))
  expect_error(statmod_profile(fl, lab[2L]), "kinked penalty")
  dd$x3 <- dd$x1
  fa <- suppressWarnings(
    statmod(y ~ x1 + x3, distributions7::gaussian1_distrib(), dd))
  expect_error(statmod_profile(fa, fa@aliased[1L]), "aliased")
})

test_that("the positions of one break-point term stay ordered", {
  set.seed(5)
  dd <- data.frame(w = runif(300, 0, 10))
  dd$y <- 0.5 * dd$w + 1.2 * pmax(dd$w - 3, 0) - 1.6 * pmax(dd$w - 7, 0) +
    stats::rnorm(300, sd = 0.5)
  fit <- statmod(y ~ seg(w, npsi = 2), distributions7::gaussian1_distrib(), dd)
  est <- fit@coefficients$mu[c("seg.psi1", "seg.psi2")]
  prof <- statmod_profile(
    fit, c("mu:seg.psi1", "mu:seg.psi2"),
    values = list("mu:seg.psi1" = c(2.5, 3, 5, 7.5),
                  "mu:seg.psi2" = c(3, 6.5, 7)))
  # the crossed and the equal combinations are not evaluated
  expect_true(all(prof[["mu:seg.psi1"]] < prof[["mu:seg.psi2"]]))
  expect_identical(nrow(prof), 7L)
  expect_true(all(is.finite(prof$profile_loglik)))
  one <- statmod_profile(fit, "mu:seg.psi1", profile_range(width = 1, n = 3))
  expect_identical(which.max(one$profile_loglik), 2L)
  expect_error(
    statmod_profile(fit, c("mu:seg.psi1", "mu:seg.psi2"),
                    list("mu:seg.psi1" = 8, "mu:seg.psi2" = 2)),
    "increasing order")
})

test_that("two held positions reproduce lm() with the positions fixed", {
  set.seed(5)
  dd <- data.frame(w = runif(300, 0, 10))
  dd$y <- 0.5 * dd$w + 1.2 * pmax(dd$w - 3, 0) - 1.6 * pmax(dd$w - 7, 0) +
    stats::rnorm(300, sd = 0.5)
  fit <- statmod(y ~ seg(w, npsi = 2), distributions7::gaussian1_distrib(), dd,
                 outer_criterion = NULL)
  prof <- statmod_profile(
    fit, c("mu:seg.psi1", "mu:seg.psi2"),
    values = list("mu:seg.psi1" = c(2.8, 3.2), "mu:seg.psi2" = c(6.7, 7.3)))
  ref <- apply(prof[1:2], 1L, function(p) {
    as.numeric(stats::logLik(
      stats::lm(y ~ w + pmax(w - p[1L], 0) + pmax(w - p[2L], 0), dd)))
  })
  expect_equal(prof$profile_loglik, unname(ref), tolerance = 1e-9)
})

test_that("a point settled by the outer search is reported as converged", {
  skip_on_cran()
  # the model and the point where outer_fit() stops at the optimum with the
  # flag of its stopping rule FALSE and the decrement below the tolerance
  # (measured: the criterion there equals, to 1e-6, that of the neighboring
  # points whose flag is TRUE)
  set.seed(3)
  n <- 1000
  dd <- data.frame(z = runif(n), w = runif(n, 0, 10), x = runif(n))
  psi <- 5 + 1.5 * dd$x
  dd$y <- sin(2 * pi * dd$z) + 0.8 * pmax(dd$w - psi, 0) +
    stats::rnorm(n, sd = exp(-0.5 + 0.4 * dd$x))
  fit <- statmod(y ~ s(z) + seg(w, psi ~ x) | sigma ~ s(x),
                 distributions7::gaussian1_distrib(), dd)
  lab <- "mu:seg.psi1.x"
  ctx <- profile_context(fit, lab, "held")
  v <- ctx$estimate + sqrt(diag(stats::vcov(fit)))[[lab]] * 2 * (-1 / 3)
  spec <- profile_hold(ctx$spec, ctx$param, ctx$coefname, v)
  design <- statmod_design(spec)
  beta0 <- ctx$beta_fit
  beta0[ctx$at] <- v
  pins <- pin_held_breakpoints(spec, design, beta0, from = ctx$beta_fit)
  res <- outer_fit(spec, design, statmod_blocks(spec, design), ctx$hyper0,
                   ctx$inner, ctx$outer, NULL, pins$beta, ctx$emethod,
                   ctx$cfg$maxit, ctx$cfg$tol, verbosity(0), hold_hyper = TRUE)
  skip_if(isTRUE(res$converged) || !isTRUE(res$settled),
          "the search no longer stops short of its rule at this point")
  prof <- statmod_profile(fit, lab, list(v))
  expect_true(prof$converged)
  expect_equal(prof$profile_pen_loglik, -res$value, tolerance = 1e-8)
  # the maximum row keeps the flag of the stopping rule
  expect_identical(attr(prof, "maximum")$converged, isTRUE(fit@converged))
})

test_that("a model with a structural term is rejected", {
  skip_on_cran()
  set.seed(7)
  dd <- data.frame(t = 1:120)
  dd$y <- stats::arima.sim(list(ar = 0.5), 120)
  fit <- suppressWarnings(
    statmod(y ~ modelterms7::gas(p = 1, q = 0),
            distributions7::gaussian1_distrib(), dd))
  expect_error(statmod_profile(fit, "mu:(Intercept)"), "structural term")
})

test_that("cores gives the same table", {
  p <- poisson_fit()
  v <- profile_range(width = 1, n = 5)
  a <- statmod_profile(p$fit, "mu:x1", v)
  b <- statmod_profile(p$fit, "mu:x1", v, cores = 2)
  expect_equal(a, b)
})
