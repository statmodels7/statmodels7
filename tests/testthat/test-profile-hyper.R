# statmod_profile() with hyperparameters in `which`. The references are fits
# whose formula carries the value in the term (hyper = c(lambda = value)),
# which rebuilds the term, and the criterion read at their penalized mode.

hyper_data <- function() {
  set.seed(5)
  n <- 250
  dd <- data.frame(z = runif(n), w = runif(n), x = rnorm(n))
  dd$y <- rpois(n, exp(1 + sin(2 * pi * dd$z) + 0.5 * dd$w + 0.4 * dd$x))
  dd
}

# the criterion of a fit at its own mode, on the sign and scale of the fit
crit_at <- function(f, meth) {
  d <- statmod_design(f@spec)
  statmod_marginal(f@spec, d, f@coefficients, f@hyper, meth,
                   basis = integrated_basis(f@spec, d, meth@kind))$value
}

test_that("one held hyperparameter gives the criterion of a refit from the formula", {
  dd <- hyper_data()
  fam <- distributions7::poisson_distrib()
  fit <- statmod(y ~ s(z) + x, fam, dd)
  key <- "mu/s(z)/lambda"
  v <- profile_range(width = 2, n = 5)
  g <- statmod_profile(fit, key, v)
  expect_identical(names(g@points)[1L], key)
  ref <- vapply(g@points[[key]], function(e) {
    f <- statmod(y ~ s(z, hyper = c(lambda = exp(e))) + x, fam, dd)
    crit_at(f, fit@methods$outer)
  }, numeric(1))
  expect_equal(g@points$outer_criterion, ref, tolerance = 1e-9)
  # the grid is on the scale of the link, and holds the estimate
  est <- log(fit@hyper$mu[["s(z)"]][["lambda"]])
  expect_equal(g@points[[key]][3L], est)
  # at the estimate the criterion is the fit's
  expect_equal(g@points$outer_criterion[3L], fit@criterion[[1L]], tolerance = 1e-9)
  expect_equal(g@unrestricted_point$outer_criterion, fit@criterion[[1L]])
  expect_equal(g@unrestricted_point[[key]], est)
  # one hyperparameter: held and refit coincide
  r <- statmod_profile(fit, key, v, hyper = "refit")
  expect_equal(r@points$outer_criterion, g@points$outer_criterion, tolerance = 1e-9)
  expect_true(all(g@points$converged))
})

test_that("the other hyperparameters are held or searched again", {
  dd <- hyper_data()
  fam <- distributions7::poisson_distrib()
  fit <- statmod(y ~ s(z) + s(w) + x, fam, dd)
  key <- "mu/s(z)/lambda"
  v <- profile_range(width = 1.5, n = 3)
  ph <- statmod_profile(fit, key, v, hyper = "held")
  pr <- statmod_profile(fit, key, v, hyper = "refit")
  lw <- fit@hyper$mu[["s(w)"]][["lambda"]]
  refh <- vapply(ph@points[[key]], function(e) {
    f <- statmod(y ~ s(z, hyper = c(lambda = exp(e))) +
                   s(w, hyper = c(lambda = lw)) + x, fam, dd)
    crit_at(f, fit@methods$outer)
  }, numeric(1))
  expect_equal(ph@points$outer_criterion, refh, tolerance = 1e-9)
  refr <- vapply(pr@points[[key]], function(e) {
    f <- statmod(y ~ s(z, hyper = c(lambda = exp(e))) + s(w) + x, fam, dd)
    as.numeric(f@criterion[[1L]])
  }, numeric(1))
  expect_equal(pr@points$outer_criterion, refr, tolerance = 1e-9)
  # searching the others again cannot lower the criterion
  expect_true(all(pr@points$outer_criterion >= ph@points$outer_criterion - 1e-9))
})

test_that("a hyperparameter and a coefficient together", {
  dd <- hyper_data()
  fam <- distributions7::poisson_distrib()
  fit <- statmod(y ~ s(z) + s(w) + x, fam, dd)
  key <- "mu/s(z)/lambda"
  b <- coef(fit)$mu[["x"]]
  vals <- stats::setNames(list(c(0.5, 1.5), c(b - 0.1, b + 0.1)), c(key, "mu:x"))
  pm <- statmod_profile(fit, c(key, "mu:x"), vals, hyper = "refit")
  expect_identical(nrow(pm@points), 4L)
  # the reference carries the coefficient as an offset
  ref <- apply(pm@points[, 1:2], 1L, function(r) {
    d2 <- dd
    d2$off <- unname(r[2L]) * d2$x
    f <- statmod(y ~ s(z, hyper = c(lambda = exp(unname(r[1L])))) + s(w) +
                   offset(off), fam, d2)
    as.numeric(f@criterion[[1L]])
  })
  expect_equal(pm@points$outer_criterion, unname(ref), tolerance = 1e-9)
  # with the other hyperparameter held at its value
  lw <- fit@hyper$mu[["s(w)"]][["lambda"]]
  ph <- statmod_profile(fit, c(key, "mu:x"), vals, hyper = "held")
  refh <- apply(ph@points[, 1:2], 1L, function(r) {
    d2 <- dd
    d2$off <- unname(r[2L]) * d2$x
    f <- statmod(y ~ s(z, hyper = c(lambda = exp(unname(r[1L])))) +
                   s(w, hyper = c(lambda = lw)) + offset(off), fam, d2)
    crit_at(f, fit@methods$outer)
  })
  expect_equal(ph@points$outer_criterion, unname(refh), tolerance = 1e-9)
})

test_that("a hyperparameter at the edge takes the grid of the fallback", {
  set.seed(2)
  dd <- data.frame(z = runif(120), x = rnorm(120))
  dd$y <- 1 + 0.5 * dd$x + stats::rnorm(120, sd = 0.5)
  fit <- statmod(y ~ s(z) + x, distributions7::gaussian1_distrib(), dd)
  key <- "mu/s(z)/lambda"
  expect_gt(log(fit@hyper$mu[["s(z)"]][["lambda"]]), 8)
  g <- statmod_profile(fit, key, profile_range(width = 2, n = 5, fallback = 3))
  expect_equal(g@points[[key]], c(5, 6.5, 8))
  expect_true(all(is.finite(g@points$outer_criterion)))
  # a hyperparameter inside the interval takes the standard error
  set.seed(3)
  d2 <- data.frame(z = runif(250), x = rnorm(250))
  d2$y <- rpois(250, exp(1 + sin(2 * pi * d2$z) + 0.4 * d2$x))
  f2 <- statmod(y ~ s(z) + x, distributions7::poisson_distrib(), d2)
  g2 <- statmod_profile(f2, key, profile_range(width = 1, n = 3))
  expect_gt(diff(range(g2@points[[key]])), 0)
  expect_lt(diff(range(g2@points[[key]])), 4)
})

test_that("a hold on the specification reaches every reader of the index", {
  dd <- hyper_data()
  fit <- statmod(y ~ s(z) + s(w) + x, distributions7::poisson_distrib(), dd)
  sp <- fit@spec
  design <- statmod_design(sp)
  expect_identical(nrow(outer_hyper_index(sp, statmod_blocks(sp, design))), 2L)
  ctx <- profile_context(fit, "mu/s(z)/lambda", "held")
  sp2 <- profile_hold_hyper(sp, ctx$hyp, 0)
  keys <- statmod_penalty_keys(sp2)
  held <- Filter(function(u) length(u$fixed), keys)
  expect_length(held, 1L)
  expect_equal(held[[1L]]$fixed$lambda, 1)
  expect_identical(nrow(outer_hyper_index(sp2, statmod_blocks(sp2, design))), 1L)
  # the specification of the fit is unchanged
  expect_length(sp@held_hyper, 0L)
})

test_that("the errors name what cannot be held", {
  dd <- hyper_data()
  fam <- distributions7::poisson_distrib()
  fit <- statmod(y ~ s(z) + x, fam, dd)
  expect_error(statmod_profile(fit, "mu/s(w)/lambda"),
               "mu/s\\(z\\)/lambda")
  nofit <- statmod(y ~ s(z) + x, fam, dd, outer_criterion = NULL)
  expect_error(statmod_profile(nofit, "mu/s(z)/lambda"), "Not a coefficient")
  expect_error(profile_range(fallback = 0), "fallback")
})

test_that("the approximate methods hold a hyperparameter", {
  dd <- hyper_data()
  fam <- distributions7::poisson_distrib()
  fit <- statmod(y ~ s(z) + s(w) + x, fam, dd)
  key <- "mu/s(z)/lambda"
  v <- profile_range(width = 2, n = 5)
  # the others held: nothing is unknown, and the criterion is exact
  g <- statmod_profile(fit, key, v, hyper = "held")
  for (m in list(profile_continuation(), profile_taylor(order = 2))) {
    a <- statmod_profile(fit, key, v, hyper = "held", method = m)
    expect_equal(a@points$outer_criterion, g@points$outer_criterion, tolerance = 1e-9)
  }
  # the others searched again: the continuation follows the grid and the
  # expansion approaches it as the order grows
  gr <- statmod_profile(fit, key, v, hyper = "refit")
  cc <- statmod_profile(fit, key, v, hyper = "refit",
                        method = profile_continuation(order = 2))
  expect_equal(cc@points$outer_criterion, gr@points$outer_criterion, tolerance = 1e-9)
  e <- vapply(c(2, 3, 5), function(o) {
    t <- statmod_profile(fit, key, v, hyper = "refit",
                         method = profile_taylor(order = o))
    expect_true(all(is.na(t@points$converged)))
    max(abs(t@points$outer_criterion - gr@points$outer_criterion))
  }, numeric(1))
  expect_true(all(diff(e) < 0))
  expect_lt(e[3L], 1e-5)
  # with a coefficient beside it
  b <- coef(fit)$mu[["x"]]
  vals <- stats::setNames(list(c(0.5, 1.5), c(b - 0.1, b + 0.1)),
                          c(key, "mu:x"))
  gm <- statmod_profile(fit, c(key, "mu:x"), vals, hyper = "refit")
  cm <- statmod_profile(fit, c(key, "mu:x"), vals, hyper = "refit",
                        method = profile_continuation())
  expect_equal(cm@points$outer_criterion, gm@points$outer_criterion, tolerance = 1e-9)
  tm <- statmod_profile(fit, c(key, "mu:x"), vals, hyper = "refit",
                        method = profile_taylor(order = 3))
  expect_lt(max(abs(tm@points$outer_criterion - gm@points$outer_criterion)), 5e-3)
})

test_that("the hyperparameter of a kinked penalty is profiled by the grid", {
  set.seed(1)
  n <- 200
  dd <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n), x4 = rnorm(n))
  dd$y <- 1 + 0.8 * dd$x1 + 0.3 * dd$x2 + rnorm(n, sd = 0.5)
  fam <- distributions7::gaussian1_distrib()
  fit <- statmod(y ~ lasso(x1 + x2 + x3 + x4), fam, dd)
  key <- "mu/lasso(x1 + x2 + x3 + x4)/lambda"
  design <- statmod_design(fit@spec)
  tab <- profile_hyper_table(fit, fit@spec, statmod_blocks(fit@spec, design))
  expect_identical(tab$key, key)
  expect_true(tab$path)
  # a grid wide enough for the lasso to move the coefficients
  v <- c(-4, 4, 6, 7)
  g <- statmod_profile(fit, key, list(v))
  expect_true(all(is.na(g@points$outer_criterion)))
  ref <- vapply(g@points[[key]], function(e) {
    r <- statmod(y ~ lasso(x1 + x2 + x3 + x4, lambda = exp(e)), fam, dd)
    c(r@loglik, -r@objective)
  }, numeric(2))
  expect_equal(g@points$profile_loglik, ref[1L, ], tolerance = 1e-5)
  expect_equal(g@points$profile_pen_loglik, ref[2L, ], tolerance = 1e-5)
  expect_gt(diff(range(g@points$profile_loglik)), 1)
  # the approximate methods name the key and the reason
  for (m in list(profile_continuation(), profile_taylor())) {
    expect_error(statmod_profile(fit, key, list(v), method = m), "kink")
    expect_error(statmod_profile(fit, key, list(v), method = m), "lasso[(]x1")
  }
})

test_that("a hyperparameter that the term holds has no profile", {
  dd <- hyper_data()
  fam <- distributions7::poisson_distrib()
  fit <- statmod(y ~ s(z, hyper = c(lambda = 2)) + x, fam, dd)
  key <- "mu/s(z, hyper = c(lambda = 2))/lambda"
  expect_error(statmod_profile(fit, key), "holds at 2")
  expect_error(statmod_profile(fit, key), "lambda = 2[)][)]/lambda")
})

test_that("a coordinate of a covariance class is held with the others searched", {
  set.seed(4)
  n <- 300
  g <- factor(rep(1:30, each = 10))
  dd <- data.frame(x = rnorm(n), g = g)
  u <- rnorm(30)
  dd$y <- 1 + 0.5 * dd$x + u[g] + rnorm(n, sd = exp(0.2 * u[g]))
  fit <- statmod(y ~ x + random(~ 1 | u | g) | sigma ~ random(~ 1 | u | g),
                 distributions7::gaussian1_distrib(), dd)
  key <- "mu/u | g/sigma_log_sd1"
  ctx <- profile_context(fit, key, "refit")
  # the class unit carries the value held on the specification
  sp <- profile_hold_hyper(ctx$spec, ctx$hyp, 1)
  cl <- Filter(function(u) !is.null(u$class), statmod_penalty_keys(sp))
  expect_equal(cl[[1L]]$fixed$sigma_log_sd1, 1)
  val <- ctx$estimate + 1
  pr <- profile_point(ctx, val)
  expect_equal(pr$hyper$mu[["u | g"]][["sigma_log_sd1"]], val)
  ph <- profile_point(profile_context(fit, key, "held"), val)
  expect_gte(pr$criterion, ph$criterion - 1e-8)
  # at the estimate the profile reaches the criterion of the fit
  p0 <- profile_point(ctx, ctx$estimate)
  expect_equal(p0$criterion, as.numeric(fit@criterion[[1L]]), tolerance = 1e-5)
})
