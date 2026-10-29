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
  expect_identical(names(g)[1L], key)
  ref <- vapply(g[[key]], function(e) {
    f <- statmod(y ~ s(z, hyper = c(lambda = exp(e))) + x, fam, dd)
    crit_at(f, fit@methods$outer)
  }, numeric(1))
  expect_equal(g$outer_criterion, ref, tolerance = 1e-9)
  # the grid is on the scale of the link, and holds the estimate
  est <- log(fit@hyper$mu[["s(z)"]][["lambda"]])
  expect_equal(g[[key]][3L], est)
  # at the estimate the criterion is the fit's
  expect_equal(g$outer_criterion[3L], fit@criterion[[1L]], tolerance = 1e-9)
  expect_equal(attr(g, "maximum")$outer_criterion, fit@criterion[[1L]])
  expect_equal(attr(g, "maximum")[[key]], est)
  # one hyperparameter: held and refit coincide
  r <- statmod_profile(fit, key, v, hyper = "refit")
  expect_equal(r$outer_criterion, g$outer_criterion, tolerance = 1e-9)
  expect_true(all(g$converged))
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
  refh <- vapply(ph[[key]], function(e) {
    f <- statmod(y ~ s(z, hyper = c(lambda = exp(e))) +
                   s(w, hyper = c(lambda = lw)) + x, fam, dd)
    crit_at(f, fit@methods$outer)
  }, numeric(1))
  expect_equal(ph$outer_criterion, refh, tolerance = 1e-9)
  refr <- vapply(pr[[key]], function(e) {
    f <- statmod(y ~ s(z, hyper = c(lambda = exp(e))) + s(w) + x, fam, dd)
    as.numeric(f@criterion[[1L]])
  }, numeric(1))
  expect_equal(pr$outer_criterion, refr, tolerance = 1e-9)
  # searching the others again cannot lower the criterion
  expect_true(all(pr$outer_criterion >= ph$outer_criterion - 1e-9))
})

test_that("a hyperparameter and a coefficient together", {
  dd <- hyper_data()
  fam <- distributions7::poisson_distrib()
  fit <- statmod(y ~ s(z) + s(w) + x, fam, dd)
  key <- "mu/s(z)/lambda"
  b <- coef(fit)$mu[["x"]]
  vals <- stats::setNames(list(c(0.5, 1.5), c(b - 0.1, b + 0.1)), c(key, "mu:x"))
  pm <- statmod_profile(fit, c(key, "mu:x"), vals, hyper = "refit")
  expect_identical(nrow(pm), 4L)
  # the reference carries the coefficient as an offset
  ref <- apply(pm[, 1:2], 1L, function(r) {
    d2 <- dd
    d2$off <- unname(r[2L]) * d2$x
    f <- statmod(y ~ s(z, hyper = c(lambda = exp(unname(r[1L])))) + s(w) +
                   offset(off), fam, d2)
    as.numeric(f@criterion[[1L]])
  })
  expect_equal(pm$outer_criterion, unname(ref), tolerance = 1e-9)
  # with the other hyperparameter held at its value
  lw <- fit@hyper$mu[["s(w)"]][["lambda"]]
  ph <- statmod_profile(fit, c(key, "mu:x"), vals, hyper = "held")
  refh <- apply(ph[, 1:2], 1L, function(r) {
    d2 <- dd
    d2$off <- unname(r[2L]) * d2$x
    f <- statmod(y ~ s(z, hyper = c(lambda = exp(unname(r[1L])))) +
                   s(w, hyper = c(lambda = lw)) + offset(off), fam, d2)
    crit_at(f, fit@methods$outer)
  })
  expect_equal(ph$outer_criterion, unname(refh), tolerance = 1e-9)
})

test_that("a hyperparameter at the edge takes the grid of the fallback", {
  set.seed(2)
  dd <- data.frame(z = runif(120), x = rnorm(120))
  dd$y <- 1 + 0.5 * dd$x + stats::rnorm(120, sd = 0.5)
  fit <- statmod(y ~ s(z) + x, distributions7::gaussian1_distrib(), dd)
  key <- "mu/s(z)/lambda"
  expect_gt(log(fit@hyper$mu[["s(z)"]][["lambda"]]), 8)
  g <- statmod_profile(fit, key, profile_range(width = 2, n = 5, fallback = 3))
  expect_equal(g[[key]], c(5, 6.5, 8))
  expect_true(all(is.finite(g$outer_criterion)))
  # a hyperparameter inside the interval takes the standard error
  set.seed(3)
  d2 <- data.frame(z = runif(250), x = rnorm(250))
  d2$y <- rpois(250, exp(1 + sin(2 * pi * d2$z) + 0.4 * d2$x))
  f2 <- statmod(y ~ s(z) + x, distributions7::poisson_distrib(), d2)
  g2 <- statmod_profile(f2, key, profile_range(width = 1, n = 3))
  expect_gt(diff(range(g2[[key]])), 0)
  expect_lt(diff(range(g2[[key]])), 4)
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
    expect_equal(a$outer_criterion, g$outer_criterion, tolerance = 1e-9)
  }
  # the others searched again: the continuation follows the grid and the
  # expansion approaches it as the order grows
  gr <- statmod_profile(fit, key, v, hyper = "refit")
  cc <- statmod_profile(fit, key, v, hyper = "refit",
                        method = profile_continuation(order = 2))
  expect_equal(cc$outer_criterion, gr$outer_criterion, tolerance = 1e-9)
  e <- vapply(c(2, 3, 5), function(o) {
    t <- statmod_profile(fit, key, v, hyper = "refit",
                         method = profile_taylor(order = o))
    expect_true(all(is.na(t$converged)))
    max(abs(t$outer_criterion - gr$outer_criterion))
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
  expect_equal(cm$outer_criterion, gm$outer_criterion, tolerance = 1e-9)
  tm <- statmod_profile(fit, c(key, "mu:x"), vals, hyper = "refit",
                        method = profile_taylor(order = 3))
  expect_lt(max(abs(tm$outer_criterion - gm$outer_criterion)), 5e-3)
})
