pt <- plot_test_fit()
fit <- pt$fit
dd <- pt$data

dctx <- function(data = NULL) {
  list(param = NULL, name = NULL, vars = NULL, cols = NULL, data = data)
}
dd_data <- function(type, fit, data = NULL) {
  plot_data(type, fit@spec@distrib, fit, ctx = dctx(data))
}

# a plain gaussian model, correct by construction
set.seed(11)
n1 <- 400
d1 <- data.frame(x = runif(n1))
d1$y <- 1 + 2 * d1$x + rnorm(n1, sd = 0.5)
fit1 <- statmod(y ~ x, distributions7::gaussian1_distrib(), d1)
mu1 <- fit1@fitted$mu
sg1 <- fit1@fitted$sigma

test_that("diagnostics() returns the four named graphs and autoplot() composes them", {
  d <- diagnostics(fit)
  expect_named(d, c("qq", "worm", "calibration", "residuals"))
  expect_true(all(vapply(d, ggplot2::is_ggplot, logical(1))))
  p <- autoplot(fit)
  expect_s3_class(p, "patchwork")
  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(plot(fit))
})

test_that("the QQ points are the sorted residuals and the order band is the beta band", {
  d <- dd_data(plot_qq(), fit)
  r <- residuals(fit, residual_type = "quantile")
  n <- length(r)
  expect_equal(d$points$y, sort(r), tolerance = 1e-12)
  expect_equal(d$points$x, qnorm(ppoints(n)), tolerance = 1e-12)
  i <- seq_len(n)
  expect_equal(d$band$lower, qnorm(qbeta(0.025, i, n + 1 - i)), tolerance = 1e-12)
  expect_equal(d$band$upper, qnorm(qbeta(0.975, i, n + 1 - i)), tolerance = 1e-12)
  d9 <- dd_data(plot_qq(level = 0.9), fit)
  expect_true(all(d9$band$upper < d$band$upper))
})

test_that("the simulated band is the envelope of the residuals of simulated responses", {
  d <- dd_data(plot_qq(band = "simulated", n_sim = 30, seed = 5), fit1)
  sims <- simulate(fit1, nsim = 30, seed = 5)
  R <- vapply(sims, function(v) qnorm(pnorm((v - mu1) / sg1)), numeric(n1))
  S <- apply(R, 2, sort)
  expect_equal(d$band$lower, apply(S, 1, quantile, probs = 0.025, names = FALSE),
               tolerance = 1e-9)
  expect_equal(d$band$upper, apply(S, 1, quantile, probs = 0.975, names = FALSE),
               tolerance = 1e-9)
  again <- dd_data(plot_qq(band = "simulated", n_sim = 30, seed = 5), fit1)
  expect_identical(again$band, d$band)
})

test_that("the worm plot is the QQ detrended, also within classes", {
  q <- dd_data(plot_qq(), fit)
  w <- dd_data(plot_worm(), fit)
  expect_equal(w$points$y, q$points$y - q$points$x, tolerance = 1e-12)
  expect_equal(w$band$upper, q$band$upper - q$band$x, tolerance = 1e-12)
  wc <- dd_data(plot_worm(over = "x", n_inter = 4), fit, dd)
  expect_equal(nlevels(wc$points$class), 4)
  expect_equal(nrow(wc$points), nrow(dd))
  cl <- cut(dd$x, quantile(dd$x, seq(0, 1, length.out = 5)), include.lowest = TRUE)
  r <- residuals(fit, residual_type = "quantile")
  k <- levels(wc$points$class)[2]
  nk <- sum(cl == levels(cl)[2])
  xk <- qnorm(ppoints(nk))
  expect_equal(wc$points$y[wc$points$class == k], sort(r[cl == levels(cl)[2]]) - xk,
               tolerance = 1e-12)
  expect_error(dd_data(plot_worm(over = "x"), fit), "give the fitting data")
  fac <- dd_data(plot_worm(over = "g"), fit, dd)
  expect_equal(levels(fac$points$class), levels(dd$g))
})

test_that("the residual graph uses the rank of the predicted median", {
  d <- dd_data(plot_residuals(), fit1)
  expect_equal(d$points$x, rank(mu1), tolerance = 1e-12)
  expect_equal(d$points$y, residuals(fit1), tolerance = 1e-12)
  expect_true(nrow(d$smooth) > 0)
  fi <- dd_data(plot_residuals(against = "fitted", param = "sigma"), fit1)
  expect_equal(fi$points$x, sg1)
  ix <- dd_data(plot_residuals(against = "index", smooth = FALSE), fit1)
  expect_equal(ix$points$x, seq_len(n1))
  expect_equal(nrow(ix$smooth), 0)
  cv <- dd_data(plot_residuals(against = "x"), fit1, d1)
  expect_equal(cv$points$x, d1$x)
  expect_error(plot_residuals(against = "fitted"), "needs 'param'")
  expect_error(dd_data(plot_residuals(against = "x"), fit1), "give the fitting data")
  expect_error(dd_data(plot_residuals(against = "fitted", param = "nu"), fit1),
               "must be one of")
})

test_that("the PIT values of a gaussian are pnorm of the standardized response", {
  d <- dd_data(plot_pit(bins = 10), fit1)
  u <- pnorm((d1$y - mu1) / sg1)
  expect_equal(sort(d$values$u), sort(u), tolerance = 1e-9)
  expect_equal(d$bars$density, tabulate(findInterval(u, seq(0, 1, 0.1),
                                                     rightmost.closed = TRUE),
                                        nbins = 10) / (n1 * 0.1),
               tolerance = 1e-9)
  expect_equal(sum(d$bars$density * 0.1), 1, tolerance = 1e-12)
  expect_equal(d$band$upper, rep(qbinom(0.975, n1, 0.1) / (n1 * 0.1), 10))
})

test_that("a correct model is calibrated within the band", {
  d <- dd_data(plot_calibration(), fit1)$points
  u <- pnorm((d1$y - mu1) / sg1)
  expect_equal(d$observed, vapply(d$p, function(p) mean(u <= p), 0),
               tolerance = 1e-9)
  expect_equal(d$upper, pmin(d$p + qnorm(0.975) * sqrt(d$p * (1 - d$p) / n1), 1),
               tolerance = 1e-12)
  expect_gte(mean(d$observed >= d$lower & d$observed <= d$upper), 0.85)
  dc <- dd_data(plot_calibration(over = "x", n_inter = 3), fit1, d1)$points
  expect_equal(nlevels(dc$class), 3)
})

test_that("the score graph is the contribution to the gradient on the link scale", {
  r <- d1$y - mu1
  s_mu <- dd_data(plot_score(param = "mu"), fit1)
  expect_equal(s_mu$points$y, r / sg1^2, tolerance = 1e-9)
  s_sg <- dd_data(plot_score(param = "sigma", against = "x"), fit1, d1)
  expect_equal(s_sg$points$y, r^2 / sg1^2 - 1, tolerance = 1e-9)
  expect_equal(s_sg$points$x, d1$x)
  expect_equal(nrow(s_sg$top), 5)
  expect_equal(sort(abs(s_sg$top$y), decreasing = TRUE),
               sort(abs(s_sg$points$y), decreasing = TRUE)[1:5], tolerance = 1e-12)
  expect_error(dd_data(plot_score(param = "nu"), fit1), "must be one of")
})

test_that("Cook distance and leverage agree with lm and with deletion refits", {
  set.seed(3)
  m <- 50
  d <- data.frame(x1 = rnorm(m), x2 = runif(m))
  d$y <- 1 + d$x1 - d$x2 + rnorm(m, sd = 0.7)
  d$y[7] <- d$y[7] + 4
  f <- statmod(y ~ x1 + x2, distributions7::gaussian1_distrib(), d)
  l <- lm(y ~ x1 + x2, d)
  inf <- influence_deltas(f, "mu", distributions7::observed_info(), adjust = TRUE)
  # leverage is the hat matrix diagonal
  expect_equal(inf$leverage, unname(hatvalues(l)), tolerance = 1e-5)
  expect_equal(inf$p_eff, 3, tolerance = 1e-5)
  # the adjusted step is the exact change of the coefficients
  expect_equal(unname(inf$delta), unname(lm.influence(l)$coefficients),
               tolerance = 1e-4)
  # Cook: the exact relation is D = cooks.distance * s^2 / sigma_fit^2
  s2 <- sum(residuals(l)^2) / l$df.residual
  expect_equal(inf$cook, unname(cooks.distance(l)) * s2 / f@fitted$sigma[1]^2,
               tolerance = 1e-4)
  # brute force: refit without observations 7 and 20
  for (i in c(7, 20)) {
    fi <- statmod(y ~ x1 + x2, distributions7::gaussian1_distrib(), d[-i, ])
    expect_equal(unname(f@coefficients$mu - fi@coefficients$mu),
                 unname(inf$delta[i, ]), tolerance = 1e-4)
  }
  # without the adjustment the step is V x s, smaller by 1 - h
  raw <- influence_deltas(f, "mu", distributions7::observed_info(), adjust = FALSE)
  expect_equal(raw$delta * (1 / (1 - inf$leverage)), inf$delta, tolerance = 1e-10)
  g <- dd_data(plot_influence(param = "mu", measure = "leverage"), f)
  expect_equal(g$points$y, inf$leverage)
  expect_equal(dd_data(plot_influence(param = "mu"), f)$top$label[1], "7")
})

test_that("the autocorrelation graph is stats::acf of the residuals", {
  set.seed(8)
  ts <- data.frame(y = as.numeric(arima.sim(list(ar = 0.6), n = 120)))
  f <- statmod(y ~ 1, distributions7::gaussian1_distrib(), ts)
  r <- residuals(f)
  d <- dd_data(plot_acf(lag_max = 15), f)
  expect_equal(d$acf$acf, as.numeric(stats::acf(r, lag.max = 15, plot = FALSE)$acf)[-1],
               tolerance = 1e-12)
  expect_equal(d$limit$value, 1.96 / sqrt(120))
  expect_gt(d$acf$acf[1], d$limit$value)
  d2 <- dd_data(plot_acf(lag_max = 15, squared = TRUE), f)
  expect_equal(d2$acf$acf,
               as.numeric(stats::acf(r^2, lag.max = 15, plot = FALSE)$acf)[-1],
               tolerance = 1e-12)
})

test_that("a discrete family is randomized and the seed makes it reproducible", {
  set.seed(4)
  dp <- data.frame(x = runif(150))
  dp$y <- rpois(150, exp(1 + dp$x))
  fp <- statmod(y ~ x, distributions7::poisson_distrib(), dp)
  a <- dd_data(plot_qq(seed = 2), fp)
  b <- dd_data(plot_qq(seed = 2), fp)
  expect_identical(a$points, b$points)
  expect_equal(a$points$y, sort(residuals(fp, seed = 2)), tolerance = 1e-12)
  p1 <- dd_data(plot_pit(seed = 3), fp)
  expect_identical(p1$values, dd_data(plot_pit(seed = 3), fp)$values)
  expect_false(identical(p1$values, dd_data(plot_pit(seed = 4), fp)$values))
  # counts per value: the observed column is the table of the response
  pc <- dd_data(plot_ppc(n_sim = 40, seed = 1), fp)
  tb <- table(dp$y)
  expect_equal(pc$counts$observed[match(as.integer(names(tb)), pc$counts$value)],
               as.numeric(tb))
  expect_true(all(pc$counts$lower <= pc$counts$mean &
                    pc$counts$mean <= pc$counts$upper))
})

test_that("the posterior predictive check is reproducible and reads the statistic", {
  a <- dd_data(plot_ppc(n_sim = 20, seed = 9), fit1)
  b <- dd_data(plot_ppc(n_sim = 20, seed = 9), fit1)
  expect_identical(a, b)
  expect_equal(length(unique(a$sim$sim)), 20)
  expect_equal(nrow(a$observed), 256)
  s <- dd_data(plot_ppc(statistic = max, n_sim = 25, seed = 9), fit1)
  sims <- simulate(fit1, nsim = 25, seed = 9)
  expect_equal(s$sim$value, vapply(sims, max, 0), ignore_attr = TRUE)
  expect_equal(s$observed$value, max(d1$y))
  expect_error(plot_ppc(statistic = 1), "must be NULL or a function")
  expect_error(dd_data(plot_ppc(statistic = function(y) y, n_sim = 3), fit1),
               "single number")
})

test_that("every graph builds and prints", {
  pdf(NULL)
  on.exit(dev.off())
  types <- list(plot_qq(band = "simulated", n_sim = 9, seed = 1),
                plot_worm(over = "x", n_inter = 3),
                plot_residuals(against = "x"), plot_residuals(against = "g"),
                plot_pit(), plot_calibration(over = "g"),
                plot_score(param = "sigma", against = "x"),
                plot_influence(param = "mu"), plot_acf(squared = TRUE),
                plot_ppc(n_sim = 9, seed = 1),
                plot_ppc(statistic = median, n_sim = 9, seed = 1))
  for (ty in types) {
    p <- autoplot(fit, type = ty, data = dd)
    expect_true(ggplot2::is_ggplot(p), info = ty@label)
    expect_no_error(print(p))
  }
})

test_that("arguments are checked and a multivariate family is not offered", {
  expect_error(plot_qq(level = 2), "between 0 and 1")
  expect_error(plot_qq(band = "x"))
  expect_error(plot_worm(n_inter = 1), "at least 2")
  expect_error(plot_pit(bins = 1), "at least 2")
  expect_error(plot_calibration(probs = c(0, 0.5)), "between 0 and 1")
  expect_error(plot_acf(lag_max = 0), "at least 1")
  expect_error(plot_qq(seed = "a"), "single number")
  expect_error(plot_influence(param = "mu", measure = "x"))
  expect_error(autoplot(fit, term = "z", type = plot_qq()), "takes neither")
  for (cl in list(PlotQq, PlotWorm, PlotResiduals, PlotPit, PlotCalibration,
                  PlotScore, PlotInfluence, PlotAcf, PlotPpc)) {
    m <- S7::method(plot_data, list(cl, distributions7::multivariate_distrib))
    expect_error(m(plot_qq(), NULL, NULL, NULL), "multivariate")
  }
})
