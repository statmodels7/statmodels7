# The repairs of 2026-10-07 found while measuring the limits of the book's
# chapter 14 (misure_cap14/RISULTATI.md).

test_that("an intercept-only nu run to its limit starts at the data-based value", {
  # a t on a response lighter-tailed than a gaussian: the intercept-only fit
  # puts nu past e^16, where a marginal criterion is flat in it
  set.seed(5)
  d <- data.frame(y = runif(400, 0, 10))
  spec <- statmod_spec(y ~ 1, distributions7::student_t1_distrib(), d, NULL, NULL)
  fd <- distributions7::fit_distrib(spec@distrib, d$y)
  skip_if(log(fd@coefficients[["nu"]]) < intercept_chart_limit(),
          "the intercept-only nu stayed inside the chart here")
  e <- statmod_intercepts(spec)
  ds <- distributions7::distrib_start(spec@distrib, d$y, 1L)[[1L]]$nu
  expect_equal(e$nu, log(ds))
  expect_lt(abs(e$nu), intercept_chart_limit())
  # a location on the identity chart is an estimate whatever its size
  d2 <- data.frame(y = rnorm(400, 1e8, 2))
  s2 <- statmod_spec(y ~ 1, distributions7::gaussian1_distrib(), d2, NULL, NULL)
  expect_gt(statmod_intercepts(s2)$mu, 1e7)
})

test_that("a coefficient a lasso holds at zero carries no variance into predict()", {
  set.seed(8)
  n <- 300
  d <- data.frame(t = runif(n), x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  d$y <- sin(2 * pi * d$t) + 0.8 * d$x1 + rnorm(n, sd = 0.3)
  fit <- statmod(y ~ s(t, bspline_smooth(k = 6)) + lasso(~ x1 + x2 + x3),
                 distributions7::gaussian1_distrib(), d,
                 sparse_criterion = bic())
  cf <- coef(fit)$mu
  zero <- names(cf)[grepl("^lasso", names(cf)) & cf == 0]
  skip_if(!length(zero), "no lasso coefficient was set to zero here")
  nd <- data.frame(t = c(0.2, 0.7), x1 = 0.5, x2 = 1, x3 = -1)
  pr <- predict(fit, what = "mu", newdata = nd, se = TRUE)
  expect_true(all(is.finite(pr$se)))
  # the covariates of the zeroed coefficients do not move the standard error
  nd0 <- nd
  for (v in sub("^lasso\\.", "", zero)) nd0[[v]] <- 0
  expect_equal(predict(fit, what = "mu", newdata = nd0, se = TRUE)$se, pr$se)
})

test_that("predict() of one parameter does not need the other equations' covariates", {
  set.seed(9)
  n <- 200
  d <- data.frame(x = runif(n), z = runif(n))
  d$y <- 1 + 2 * d$x + rnorm(n, sd = exp(-1 + d$z))
  fit <- statmod(y ~ x | sigma ~ z, distributions7::gaussian1_distrib(), d)
  a <- predict(fit, what = "mu", newdata = data.frame(x = c(0.2, 0.8)), se = TRUE)
  b <- predict(fit, what = "mu", newdata = data.frame(x = c(0.2, 0.8), z = 0.5),
               se = TRUE)
  expect_equal(a, b)
  # the parameter whose equation reads z still needs it
  expect_error(predict(fit, what = "sigma", newdata = data.frame(x = 0.2)))
  expect_identical(predict_unread(fit, "mu", data.frame(x = 0.2)), "sigma")
  expect_identical(predict_unread(fit, "mu", data.frame(x = 0.2, z = 1)),
                   character(0))
  expect_identical(predict_unread(fit, "parameter", data.frame(x = 0.2)),
                   character(0))
})

test_that("plot() on a fit names what to use instead", {
  set.seed(10)
  d <- data.frame(x = runif(50))
  d$y <- 1 + d$x + rnorm(50, sd = 0.2)
  fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), d)
  expect_error(plot(fit), "no plot() method", fixed = TRUE)
  expect_error(plot(fit), "predict(fit, newdata = grid)", fixed = TRUE)
})

test_that("Pearson residuals are NA, with a warning, where the variance fails", {
  set.seed(11)
  d <- data.frame(x = runif(60))
  d$y <- 1 + d$x + rnorm(60, sd = 0.2)
  fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), d)
  local_mocked_bindings(std_dev = function(x, theta, ...) stop("no variance"),
                        .package = "distributions7")
  expect_warning(r <- residuals(fit, type = "pearson"),
                 "could not be computed")
  expect_true(all(is.na(r)))
})

test_that("a family's own method without ... is named before the fit", {
  Lap <- S7::new_class("LapNoDots", parent = distributions7::continuous_distrib,
                       package = NULL)
  pdf_generic <- distributions7::distrib_pdf
  S7::method(pdf_generic, Lap) <- function(distrib, y, theta, log = FALSE) {
    ld <- -log(2 * theta[[2]]) - abs(y - theta[[1]]) / theta[[2]]
    if (log) ld else exp(ld)
  }
  d <- Lap(distrib_name = "laplace without dots", dimension = "univariate",
           bounds = c(-Inf, Inf), params = c("mu", "b"),
           params_interpretation = c(mu = "location", b = "scale"), n_params = 2,
           params_bounds = list(mu = c(-Inf, Inf), b = c(0, Inf)),
           link_params = list(mu = linkfunctions7::identity_link(),
                              b = linkfunctions7::log_link()))
  set.seed(12)
  dd <- data.frame(x = runif(40))
  dd$y <- dd$x + rnorm(40)
  expect_error(statmod(y ~ x, distrib = d, data = dd),
               "distrib_pdf() method of 'laplace without dots' has no `...`",
               fixed = TRUE)
  expect_error(statmod(y ~ x, distrib = d, data = dd),
               "function(distrib, y, theta, log = FALSE, ...)", fixed = TRUE)
})
