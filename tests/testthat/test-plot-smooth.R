pt <- plot_test_fit()
fit <- pt$fit
dd <- pt$data
tm <- fit@spec@terms$mu[["s(x)"]]
ctx <- list(param = "mu", name = "s(x)", vars = "x", cols = NULL, data = dd)

# the derivative of the block of the Demmler-Reinsch smooth, in closed form:
# the free column is (x - mean) / sd and the rest is the derivative of the
# rotated basis
analytic_rows <- function(tm, x, order, sdx) {
  bp <- tm@build_info$core$build_info
  B <- basis7::basis_deriv(bp$basis, x, order = order)
  cbind(if (order == 1) 1 / sdx else 0, B)
}

test_that("the derivative rows are the closed-form derivative of the basis", {
  d <- plot_data(plot_derivative(), tm, fit, ctx)
  A <- analytic_rows(tm, d$curve$x, 1, stats::sd(dd$x))
  cb <- term_coef_block(fit, "mu", "s(x)", NULL, bayesian_vcov())
  expect_equal(d$curve$fit, drop(A %*% cb$coef), tolerance = 1e-6)
  d2 <- plot_data(plot_derivative(order = 2), tm, fit, ctx)
  A2 <- analytic_rows(tm, d2$curve$x, 2, stats::sd(dd$x))
  expect_equal(d2$curve$fit, drop(A2 %*% cb$coef), tolerance = 1e-5)
})

test_that("the derivative is the slope of the effect curve", {
  d <- plot_data(plot_derivative(n_grid = 100), tm, fit, ctx)
  e <- plot_data(plot_effect(n_grid = 100), tm, fit, ctx)$curve
  h <- e$x[2] - e$x[1]
  fd <- (e$fit[3:100] - e$fit[1:98]) / (2 * h)
  expect_lt(max(abs(fd - d$curve$fit[2:99])), 0.02 * diff(range(d$curve$fit)))
  # and a fine central difference of the prediction agrees closely
  x0 <- d$curve$x[c(10, 50, 90)]
  eps <- 1e-6
  nd <- function(x) data.frame(x = x)
  num <- (modelterms7::term_predict(tm, nd(x0 + eps)) -
            modelterms7::term_predict(tm, nd(x0 - eps))) / (2 * eps)
  cb <- term_coef_block(fit, "mu", "s(x)", NULL, bayesian_vcov())
  expect_equal(d$curve$fit[c(10, 50, 90)], drop(num %*% cb$coef),
               tolerance = 1e-5)
})

test_that("the standard error is the delta method on the derivative rows", {
  for (ord in 1:2) {
    d <- plot_data(plot_derivative(order = ord), tm, fit, ctx)
    A <- analytic_rows(tm, d$curve$x, ord, stats::sd(dd$x))
    cb <- term_coef_block(fit, "mu", "s(x)", NULL, bayesian_vcov())
    se <- sqrt(rowSums((A %*% cb$vcov) * A))
    expect_equal(d$curve$se, se, tolerance = 1e-5)
    z <- stats::qnorm(0.975)
    expect_equal(d$curve$upper - d$curve$fit, z * d$curve$se, tolerance = 1e-12)
  }
  d <- plot_data(plot_derivative(level = 0.8), tm, fit, ctx)
  expect_equal(d$curve$upper - d$curve$fit,
               stats::qnorm(0.9) * d$curve$se, tolerance = 1e-12)
})

test_that("the simultaneous critical value exceeds the pointwise one", {
  p <- plot_derivative(simultaneous = TRUE, n_sim = 2000, seed = 3)
  d <- plot_data(p, tm, fit, ctx)
  expect_gt(d$critical$used, d$critical$pointwise)
  expect_lt(d$critical$used, 6)
  d2 <- plot_data(p, tm, fit, ctx)
  expect_identical(d$critical$used, d2$critical$used)
  # the seed leaves the generator where it was
  set.seed(10)
  a <- stats::runif(1)
  set.seed(10)
  plot_data(p, tm, fit, ctx)
  expect_identical(stats::runif(1), a)
  # the coverage of the band over the grid is the level (re-draw from V)
  cb <- term_coef_block(fit, "mu", "s(x)", NULL, bayesian_vcov())
  A <- analytic_rows(tm, d$curve$x, 1, stats::sd(dd$x))
  set.seed(5)
  E <- mvn_draws(cb$vcov, 2000)
  m <- apply(abs(A %*% t(E)) / d$curve$se, 2, max)
  expect_equal(mean(m <= d$critical$used), 0.95, tolerance = 0.02)
})

test_that("the zones are the runs in which the band excludes zero", {
  d <- plot_data(plot_derivative(), tm, fit, ctx)
  sig <- d$curve$lower > 0 | d$curve$upper < 0
  inside <- vapply(d$curve$x, function(v) any(v >= d$zones$xmin &
                                                v <= d$zones$xmax), TRUE)
  expect_identical(inside, sig)
})

test_that("a smooth with a by factor has a curve and a critical value per level", {
  f <- statmod(y ~ g + s(x, by = g), distributions7::gaussian1_distrib(), dd)
  nm <- grep("^s\\(x", names(f@spec@terms$mu), value = TRUE)
  d <- plot_data(plot_derivative(simultaneous = TRUE, n_sim = 500, seed = 1),
                 f@spec@terms$mu[[nm]], f,
                 ctx = list(param = "mu", name = nm, vars = "x", cols = NULL,
                            data = dd))
  expect_identical(levels(d$curve$group), levels(dd$g))
  expect_equal(nrow(d$critical), nlevels(dd$g))
  expect_true(all(d$critical$used > d$critical$pointwise))
})

test_that("the derivative graph rejects what it cannot draw", {
  f2 <- statmod(y ~ te(x, z), distributions7::gaussian1_distrib(), dd)
  nm <- grep("^te", names(f2@spec@terms$mu), value = TRUE)
  expect_error(autoplot(f2, term = nm, plot_type = plot_derivative(), data = dd),
               "plot_surface")
  expect_error(autoplot(fit, term = "s(x)", plot_type = plot_derivative()),
               "give the fitting data")
  expect_error(plot_derivative(order = 3), "'order'")
  expect_error(plot_derivative(n_sim = 5), "'n_sim'")
})

test_that("the surface is the contribution of the term on the grid", {
  f2 <- statmod(y ~ te(x, z), distributions7::gaussian1_distrib(), dd)
  nm <- grep("^te", names(f2@spec@terms$mu), value = TRUE)
  cx <- list(param = "mu", name = nm, vars = c("x", "z"), cols = NULL,
             data = dd)
  s <- plot_data(plot_surface(n_grid = 15), f2@spec@terms$mu[[nm]], f2, cx)
  expect_equal(nrow(s$surface), 15^2)
  # differences of the surface equal differences of predict() on one point
  i <- c(1, 77, 150, 225)
  nd <- data.frame(x = s$surface$x[i], z = s$surface$y[i])
  pr <- predict(f2, what = "mu", newdata = nd, random = "zero")
  expect_equal(diff(s$surface$fit[i]), diff(pr), tolerance = 1e-10)
  # the mask marks the cells whose interval contains zero
  zc <- stats::qnorm(0.975)
  expect_identical(s$surface$insignificant,
                   abs(s$surface$fit) < zc * s$surface$se)
  expect_equal(nrow(s$masked), sum(s$surface$insignificant))
  nm0 <- plot_data(plot_surface(n_grid = 15, mask = FALSE),
                   f2@spec@terms$mu[[nm]], f2, cx)
  expect_equal(nrow(nm0$masked), 0)
  # the slices are the surface at fixed values of the second variable
  sl <- plot_data(plot_surface(n_grid = 15, slices = c(0, 1)),
                  f2@spec@terms$mu[[nm]], f2, cx)
  expect_equal(nrow(sl$curve), 30)
  nd <- data.frame(x = sl$curve$x[c(3, 20)], z = sl$curve$y[c(3, 20)])
  pr <- predict(f2, what = "mu", newdata = nd, random = "zero")
  expect_equal(diff(sl$curve$fit[c(3, 20)]), diff(pr), tolerance = 1e-10)
  expect_equal(sl$curve$upper - sl$curve$fit, zc * sl$curve$se,
               tolerance = 1e-12)
  expect_error(plot_data(plot_surface(slices = 5), f2@spec@terms$mu[[nm]], f2,
                         cx), "range")
  expect_error(autoplot(fit, term = "s(x)", plot_type = plot_surface(), data = dd),
               "draws a function of two")
  pdf(NULL)
  on.exit(dev.off())
  expect_s3_class(print(autoplot(f2, term = nm, plot_type = plot_surface(),
                                 data = dd)), "ggplot")
})

test_that("the criterion at the fitted hyperparameter is the fit's criterion", {
  d <- plot_data(plot_criterion(log_lambda = log(fit@hyper$mu[["s(x)"]]) +
                                  c(-1, 0, 1)), tm, fit, ctx)
  expect_equal(d$fitted$value, fit@criterion, tolerance = 1e-9)
  expect_equal(d$curve$value[2], fit@criterion, tolerance = 1e-9)
  expect_equal(d$fitted$x, unname(log(fit@hyper$mu[["s(x)"]])))
  # the optimum: the neighbors are lower (the criterion is maximized)
  expect_lt(max(d$curve$value[c(1, 3)]), fit@criterion)
  # the edf at the fitted value is the stored one
  e <- fit@edf
  expect_equal(d$fitted$edf, e$edf[e$parameter == "mu" & e$term == "s(x)"],
               tolerance = 1e-6)
  # a larger smoothing parameter gives a smaller edf
  expect_lt(d$curve$edf[3], d$curve$edf[1])
  # a random effect: its scale is the hyperparameter
  dr <- plot_data(plot_criterion(log_lambda = log(fit@hyper$mu[["random(~1 | id)"]])
                                 + c(-0.5, 0, 0.5)),
                  fit@spec@terms$mu[["random(~1 | id)"]], fit,
                  list(param = "mu", name = "random(~1 | id)", vars = "id",
                       cols = NULL, data = dd))
  expect_equal(dr$curve$value[2], fit@criterion, tolerance = 1e-9)
})

test_that("a prediction-error criterion is drawn by refitting", {
  f <- statmod(y ~ z + s(x), distributions7::gaussian1_distrib(), dd,
               outer_criterion = aic())
  cx <- list(param = "mu", name = "s(x)", vars = "x", cols = NULL, data = dd)
  d <- plot_data(plot_criterion(log_lambda = log(f@hyper$mu[["s(x)"]]) +
                                  c(-1, 0, 1)),
                 f@spec@terms$mu[["s(x)"]], f, cx)
  expect_equal(d$curve$value[2], f@criterion, tolerance = 1e-9)
  expect_gt(min(d$curve$value[c(1, 3)]), f@criterion)
  expect_error(plot_data(plot_criterion(method = "quadratic"),
                         f@spec@terms$mu[["s(x)"]], f, cx), "reml\\(\\) or ml\\(\\)")
})

test_that("the quadratic route matches the refit near the optimum", {
  u0 <- log(fit@hyper$mu[["s(x)"]])
  g <- u0 + seq(-0.25, 0.25, length.out = 5)
  r <- plot_data(plot_criterion(log_lambda = g), tm, fit, ctx)$curve
  q <- plot_data(plot_criterion(log_lambda = g, method = "quadratic"), tm, fit,
                 ctx)$curve
  expect_lt(max(abs(r$value - q$value)), 0.01)
  expect_true(all(is.na(q$edf)))
  # the expansion is exact in value and slope at the fit
  expect_equal(q$value[3], fit@criterion, tolerance = 1e-12)
  # the curvature is the exact outer Hessian: compare with a refit difference
  h <- 0.05
  rr <- plot_data(plot_criterion(log_lambda = u0 + c(-h, 0, h)), tm, fit,
                  ctx)$curve$value
  qq <- plot_data(plot_criterion(log_lambda = u0 + c(-h, 0, h),
                                 method = "quadratic"), tm, fit, ctx)$curve$value
  expect_equal(rr[1] - 2 * rr[2] + rr[3], qq[1] - 2 * qq[2] + qq[3],
               tolerance = 5e-3)
})

test_that("the default grid and the errors of the criterion graph", {
  d <- plot_data(plot_criterion(method = "quadratic"), tm, fit, ctx)
  expect_equal(nrow(d$curve), 41)
  expect_equal(range(d$curve$x) - log(fit@hyper$mu[["s(x)"]]), c(-5, 5))
  fh <- statmod(y ~ x + s(z, hyper = 3), distributions7::gaussian1_distrib(), dd)
  nm <- grep("^s", names(fh@spec@terms$mu), value = TRUE)
  expect_error(autoplot(fh, term = nm, plot_type = plot_criterion(), data = dd),
               "no hyperparameter")
  expect_error(autoplot(fit, term = "s(x)", plot_type = plot_criterion(hyper = "zz")),
               "no estimated hyperparameter 'zz'")
  expect_error(plot_criterion(method = "newton"), "should be one of")
  expect_error(plot_criterion(log_lambda = 1), "log_lambda")
  fn <- statmod(y ~ s(x), distributions7::gaussian1_distrib(), dd,
               outer_criterion = NULL)
  expect_error(autoplot(fn, term = "s(x)", plot_type = plot_criterion()),
               "no outer")
  pdf(NULL)
  on.exit(dev.off())
  expect_s3_class(print(autoplot(fit, term = "s(x)", plot_type = plot_criterion(),
                                 data = dd)), "ggplot")
})
