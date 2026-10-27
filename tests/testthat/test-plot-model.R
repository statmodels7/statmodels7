pt <- plot_test_fit()
fit <- pt$fit
dd <- pt$data

set.seed(3)
dp <- data.frame(x = runif(150), g = factor(sample(c("a", "b"), 150, TRUE)))
dp$y <- rpois(150, exp(1 + dp$x + 0.3 * (dp$g == "b")))
fp <- statmod(y ~ x + g, distributions7::poisson_distrib(), dp)

test_that("plot_param equals predict() on the grid", {
  p <- autoplot(fit, type = plot_param(param = "sigma", over = "x",
                                       n_grid = 25), data = dd)
  expect_s3_class(p, "ggplot")
  d <- plot_data(plot_param(param = "sigma", over = "x", n_grid = 25),
                 fit@spec@distrib, fit,
                 ctx = list(data = dd))
  nd <- data.frame(x = d$curve$x, z = median(dd$z),
                   g = factor("a", levels(dd$g)),
                   id = factor("1", levels(dd$id)))
  pr <- predict(fit, what = "sigma", newdata = nd, se = TRUE,
                random = "zero")
  expect_equal(d$curve$fit, pr$fit, tolerance = 1e-12)
  expect_equal(d$curve$lower, pr$lower, tolerance = 1e-12)
  expect_equal(d$curve$upper, pr$upper, tolerance = 1e-12)
  # a factor gives one estimate per level
  pf <- plot_data(plot_param(param = "mu", over = "g"), fit@spec@distrib, fit,
                  ctx = list(data = dd))
  expect_equal(nrow(pf$curve), 3L)
  # at changes the other covariates
  d2 <- plot_data(plot_param(param = "mu", over = "x", n_grid = 5,
                             at = list(z = 1)), fit@spec@distrib, fit,
                  ctx = list(data = dd))
  d3 <- plot_data(plot_param(param = "mu", over = "x", n_grid = 5,
                             at = list(z = 0)), fit@spec@distrib, fit,
                  ctx = list(data = dd))
  expect_equal(d2$curve$fit - d3$curve$fit,
               rep(coef(fit, readable = FALSE)$mu[["z"]], 5), tolerance = 1e-9)
})

test_that("quantile curves are mu + sigma qnorm(p) for a gaussian", {
  pq <- plot_quantiles(over = "x", n_grid = 20)
  d <- plot_data(pq, fit@spec@distrib, fit, ctx = list(data = dd))
  nd <- data.frame(x = unique(d$curves$x), z = median(dd$z),
                   g = factor("a", levels(dd$g)), id = factor("1", levels(dd$id)))
  mu <- predict(fit, what = "mu", newdata = nd, random = "zero")
  sg <- predict(fit, what = "sigma", newdata = nd, random = "zero")
  for (p in pq@probs) {
    cv <- d$curves[d$curves$prob == p, ]
    expect_equal(cv$q, mu + sg * qnorm(p), tolerance = 1e-8)
  }
  expect_equal(nrow(d$points), nrow(dd))
  pdf(NULL)
  on.exit(dev.off())
  expect_s3_class(print(autoplot(fit, type = pq, data = dd)), "ggplot")
})

test_that("quantile curves of a Poisson fit are qpois and steps", {
  pq <- plot_quantiles(over = "x", n_grid = 15, probs = c(0.1, 0.5, 0.9))
  d <- plot_data(pq, fp@spec@distrib, fp, ctx = list(data = dp))
  nd <- data.frame(x = unique(d$curves$x), g = factor("a", levels(dp$g)))
  mu <- predict(fp, what = "mu", newdata = nd)
  for (p in pq@probs) {
    expect_equal(d$curves$q[d$curves$prob == p], qpois(p, mu))
  }
  expect_true(isTRUE(attr(d, "discrete")))
  p <- autoplot(fp, type = pq, data = dp)
  expect_true(any(vapply(p$layers, function(l) inherits(l$geom,
                                                        "GeomStep"), TRUE)))
})

test_that("densities integrate to one and masses sum to one", {
  at <- data.frame(x = c(0.1, 0.5, 0.9))
  d <- plot_data(plot_density(at = at, over = "x"), fit@spec@distrib, fit,
                 ctx = list(data = dd))
  for (i in 1:3) {
    r <- d$dens[d$dens$index == i, ]
    expect_equal(sum(diff(r$y) * (head(r$dens, -1) + tail(r$dens, -1)) / 2),
                 0.998, tolerance = 2e-3)
    expect_equal(r$x0[1], at$x[i])
    expect_equal(max(r$x - r$x0) >= 0, TRUE)
  }
  # the density is the gaussian at the fitted parameters
  nd <- data.frame(x = 0.5, z = median(dd$z), g = factor("a", levels(dd$g)),
                   id = factor("1", levels(dd$id)))
  mu <- predict(fit, what = "mu", newdata = nd, random = "zero")
  sg <- predict(fit, what = "sigma", newdata = nd, random = "zero")
  r <- d$dens[d$dens$index == 2, ]
  expect_equal(r$dens, dnorm(r$y, mu, sg), tolerance = 1e-8)
  expect_s3_class(autoplot(fit, type = plot_density(at = at), data = dd),
                  "ggplot")
  # discrete
  ap <- data.frame(x = c(0.2, 0.8), g = c("a", "b"))
  dp2 <- plot_data(plot_density(at = ap, over = "x"), fp@spec@distrib, fp,
                   ctx = list(data = dp))
  s <- tapply(dp2$dens$dens, dp2$dens$index, sum)
  expect_equal(as.numeric(s), c(0.998, 0.998), tolerance = 5e-3)
  expect_s3_class(autoplot(fp, type = plot_density(at = ap), data = dp),
                  "ggplot")
  expect_s3_class(autoplot(fp, type = plot_density(at = ap, over = "x"),
                           data = dp), "ggplot")
})

test_that("moments of a gaussian are mu and sigma squared", {
  pm <- plot_moments(over = "x", n_grid = 12)
  d <- plot_data(pm, fit@spec@distrib, fit, ctx = list(data = dd))
  nd <- data.frame(x = unique(d$curve$x), z = median(dd$z),
                   g = factor("a", levels(dd$g)), id = factor("1", levels(dd$id)))
  mu <- predict(fit, what = "mu", newdata = nd, random = "zero")
  sg <- predict(fit, what = "sigma", newdata = nd, random = "zero")
  expect_equal(d$curve$value[d$curve$moment == "mean"], mu, tolerance = 1e-10)
  expect_equal(d$curve$value[d$curve$moment == "variance"], sg^2,
               tolerance = 1e-10)
  p <- autoplot(fit, type = pm, data = dd)
  expect_s3_class(p$facet, "FacetWrap")
  expect_error(plot_moments(over = "x", moments = "median"), "must be among")
})

test_that("arguments and data are checked", {
  expect_error(plot_quantiles(), "'over'")
  expect_error(plot_param(param = "mu"), "'over'")
  expect_error(plot_moments(), "'over'")
  expect_error(plot_density(at = 1), "data frame")
  expect_error(autoplot(fit, type = plot_quantiles(over = "x")), "data")
  expect_error(autoplot(fit, type = plot_param(param = "nu", over = "x"),
                        data = dd), "'param' must be one of")
  expect_error(autoplot(fit, type = plot_quantiles(over = "x",
                                                   at = list(x = 1)),
                        data = dd), "variable on the axis")
  expect_error(autoplot(fit, term = "s(x)", type = plot_moments(over = "x"),
                        data = dd), "neither 'term' nor 'param'")
})
