pt <- plot_test_fit()
fit <- pt$fit
dd <- pt$data

# the labels of every legend of a graph, whatever the aesthetic
legend_labels <- function(p) {
  unlist(lapply(c("colour", "fill", "linetype", "shape", "alpha"), function(a) {
    g <- ggplot2::get_guide_data(p, a)
    if (is.null(g)) character(0) else as.character(g$.label)
  }))
}

test_that("the effect graph names the estimate with its CI and the residuals", {
  lab <- legend_labels(autoplot(fit, term = "s(x)", data = dd))
  expect_setequal(lab, c("estimated effect\n(95% CI)", "partial residuals"))
  p <- autoplot(fit, term = "s(x)", data = dd)
  expect_identical(p$labels$y, "effect on the link scale of mu")
  lb <- legend_labels(autoplot(fit, term = "x", param = "sigma",
                               plot_type = plot_effect(binned = 5,
                                                       level = 0.9),
                               data = dd))
  expect_true("estimated effect\n(90% CI)" %in% lb)
  expect_true(any(grepl("^mean of the partial residuals in [0-9]+ classes$",
                        lb)))
})

test_that("the reference level of a factor is labeled", {
  p <- autoplot(fit, term = "g", data = dd)
  expect_true("estimated effect\n(95% CI)" %in% legend_labels(p))
  txt <- Filter(function(l) inherits(l$geom, "GeomText"), p$layers)
  expect_length(txt, 1L)
  expect_identical(as.character(txt[[1L]]$data$x), "a")
})

test_that("an effect crossed with a factor keeps the levels", {
  fi <- statmod(y ~ z * g, distributions7::gaussian1_distrib(), dd)
  lab <- legend_labels(autoplot(fi, term = "z", data = dd))
  expect_true(all(c("a", "b", "c", "estimated effect\n(95% CI)",
                    "partial residuals") %in% lab))
})

test_that("legend = FALSE removes the legend", {
  p <- autoplot(fit, term = "s(x)", data = dd,
                plot_type = plot_effect(legend = FALSE))
  expect_identical(p$theme$legend.position, "none")
  expect_error(plot_qq(legend = NA), "TRUE or FALSE")
  expect_false(plot_worm(legend = FALSE)@legend)
})

test_that("the diagnostic graphs name their elements and not the references", {
  q <- legend_labels(autoplot(fit, plot_type = plot_qq(), data = dd))
  expect_setequal(q, c("quantile residuals",
                       "95% pointwise band of the order statistics"))
  r <- legend_labels(autoplot(fit, plot_type = plot_residuals(), data = dd))
  expect_setequal(r, c("quantile residuals", "smoothed mean of the residuals"))
  s <- legend_labels(autoplot(fit, plot_type = plot_score(param = "mu"),
                              data = dd))
  expect_identical(s, "score of an observation")
  c1 <- legend_labels(autoplot(fit, plot_type = plot_ppc(n_sim = 5, seed = 1),
                               data = dd))
  expect_setequal(c1, c("density of the observed response",
                        "densities of the simulated datasets"))
  cal <- legend_labels(autoplot(fit, plot_type = plot_calibration(),
                                data = dd))
  expect_false("exact calibration" %in% cal)
})

test_that("the panels of a graph by classes name the covariate", {
  p <- autoplot(fit, plot_type = plot_worm(over = "x", n_inter = 3), data = dd)
  strips <- ggplot2::ggplot_build(p)$layout$facet$params$labeller(
    data.frame(class = "[0,1]"))
  expect_identical(strips[[1L]], "x: [0,1]")
})

test_that("the model graphs name each quantile and a parameter with its CI", {
  q <- legend_labels(autoplot(fit, plot_type = plot_quantiles(over = "x"),
                              data = dd))
  expect_true(all(c("median", "quantiles 0.25 and 0.75",
                    "quantiles 0.05 and 0.95", "observations") %in% q))
  p <- legend_labels(autoplot(fit, plot_type = plot_param(param = "sigma",
                                                          over = "x"),
                              data = dd))
  expect_identical(p, "fitted sigma\n(95% CI)")
})

test_that("the penalized, random and smoothing graphs carry their legends", {
  s <- legend_labels(autoplot(fit, term = "s(x)",
                              plot_type = plot_shrinkage(), data = dd))
  expect_setequal(s, c("penalized estimate", "non-penalized estimate",
                       "shrinkage"))
  rn <- grep("^random", names(fit@spec@terms$mu), value = TRUE)
  cp <- autoplot(fit, term = rn, data = dd)
  expect_true("predicted effect with 95% interval" %in% legend_labels(cp))
  expect_identical(cp$labels$y, "level of id")
  cr <- legend_labels(autoplot(fit, term = "s(x)",
                               plot_type = plot_criterion(), data = dd))
  expect_true(all(c("fitted lambda", "edf of the term") %in% cr))
  dv <- legend_labels(autoplot(fit, term = "s(x)",
                               plot_type = plot_derivative(), data = dd))
  expect_identical(dv, "estimated derivative\n(95% CI)")
})

test_that("the pooled autocorrelation is the one of acf() for one series", {
  set.seed(3)
  e <- as.numeric(arima.sim(list(ar = 0.6), n = 200))
  r <- pooled_acf(e, list(seq_along(e)), 10)
  expect_equal(r, as.numeric(stats::acf(e, lag.max = 10, plot = FALSE)$acf)[-1],
               tolerance = 1e-12)
  expect_equal(durbin_levinson(r),
               as.numeric(stats::pacf(e, lag.max = 10, plot = FALSE)$acf),
               tolerance = 1e-12)
  # two series stored one after the other: no pair crosses from one to the
  # other
  s <- list(1:100, 101:200)
  manual <- (sum((e[2:100] - mean(e)) * (e[1:99] - mean(e))) +
               sum((e[102:200] - mean(e)) * (e[101:199] - mean(e)))) /
    sum((e - mean(e))^2)
  expect_equal(pooled_acf(e, s, 1), manual, tolerance = 1e-12)
})

test_that("plot_acf() reads the series from by and time", {
  dd2 <- dd
  dd2$id <- rep(1:4, each = 50)
  dd2$t <- rep(50:1, 4)
  ser <- acf_series(plot_acf(by = "id", time = "t"), fit,
                    list(data = dd2))
  expect_length(ser, 4L)
  expect_identical(ser[[1L]], 50:1)
  p <- autoplot(fit, plot_type = plot_pacf(by = "id", time = "t", lag_max = 8),
                data = dd2)
  expect_true("partial autocorrelation" %in% legend_labels(p))
  expect_error(autoplot(fit, plot_type = plot_acf(by = "id")),
               "does not keep its data")
})
