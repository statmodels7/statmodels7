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

test_that("the effect graph names its curve, band, residuals and rug", {
  lab <- legend_labels(autoplot(fit, term = "s(x)", data = dd))
  expect_true(all(c("estimated effect", "partial residuals", "observed values",
                    "95% pointwise confidence band") %in% lab))
  lb <- legend_labels(autoplot(fit, term = "x", param = "sigma",
                               plot_type = plot_effect(binned = 5,
                                                       level = 0.9),
                               data = dd))
  expect_true("90% pointwise confidence band" %in% lb)
  expect_true(any(grepl("^mean of the partial residuals in [0-9]+ classes$",
                        lb)))
  lg <- legend_labels(autoplot(fit, term = "g", data = dd))
  expect_true("estimate with 95% confidence interval" %in% lg)
})

test_that("an effect crossed with a factor keeps the levels and names the band", {
  fi <- statmod(y ~ z * g, distributions7::gaussian1_distrib(), dd)
  lab <- legend_labels(autoplot(fi, term = "z", data = dd))
  expect_true(all(c("a", "b", "c", "estimated effect", "partial residuals",
                    "95% pointwise confidence band") %in% lab))
})

test_that("the diagnostic graphs name their reference elements", {
  q <- legend_labels(autoplot(fit, plot_type = plot_qq(), data = dd))
  expect_true(all(c("quantile residuals", "line y = x",
                    "95% pointwise band of the order statistics") %in% q))
  r <- legend_labels(autoplot(fit, plot_type = plot_residuals(), data = dd))
  expect_true(all(c("quantile residuals", "smoothed mean of the residuals",
                    "-2 and 2 (95% of a standard normal)") %in% r))
  p <- legend_labels(autoplot(fit, plot_type = plot_pit(), data = dd))
  expect_true(all(c("histogram of F(y)", "uniform density") %in% p))
  a <- legend_labels(autoplot(fit, plot_type = plot_acf(), data = dd))
  expect_true(any(grepl("^95% bounds of an uncorrelated series", a)))
  s <- legend_labels(autoplot(fit, plot_type = plot_score(param = "mu"),
                              data = dd))
  expect_true("the five largest in absolute value, with the row number" %in% s)
  c1 <- legend_labels(autoplot(fit, plot_type = plot_ppc(n_sim = 5, seed = 1),
                               data = dd))
  expect_true(all(c("density of the observed response",
                    "density of each of 5 simulated responses") %in% c1))
})

test_that("the panels of a graph by classes name the covariate", {
  p <- autoplot(fit, plot_type = plot_worm(over = "x", n_inter = 3), data = dd)
  strips <- ggplot2::ggplot_build(p)$layout$facet$params$labeller(
    data.frame(class = "[0,1]"))
  expect_identical(strips[[1L]], "x: [0,1]")
})

test_that("the model graphs name each quantile and the band of a parameter", {
  q <- legend_labels(autoplot(fit, plot_type = plot_quantiles(over = "x"),
                              data = dd))
  expect_true(all(c("median", "quantiles 0.25 and 0.75",
                    "quantiles 0.05 and 0.95", "observations") %in% q))
  p <- legend_labels(autoplot(fit, plot_type = plot_param(param = "sigma",
                                                          over = "x"),
                              data = dd))
  expect_true(all(c("fitted sigma", "95% pointwise confidence band") %in% p))
})

test_that("the penalized, random and smoothing graphs carry their legends", {
  s <- legend_labels(autoplot(fit, term = "s(x)",
                              plot_type = plot_shrinkage(), data = dd))
  expect_true(all(c("penalized estimate",
                    "estimate without the penalty of the term",
                    "shrinkage by the penalty") %in% s))
  rn <- grep("^random", names(fit@spec@terms$mu), value = TRUE)
  cp <- autoplot(fit, term = rn, data = dd)
  expect_true("predicted effect with 95% interval" %in% legend_labels(cp))
  expect_identical(cp$labels$y, "level of id")
  cr <- legend_labels(autoplot(fit, term = "s(x)",
                               plot_type = plot_criterion(), data = dd))
  expect_true(all(c("fitted lambda", "edf of the term (right axis)") %in% cr))
  dv <- legend_labels(autoplot(fit, term = "s(x)",
                               plot_type = plot_derivative(), data = dd))
  expect_true(all(c("estimated derivative", "zero derivative",
                    "95% pointwise confidence band") %in% dv))
})
