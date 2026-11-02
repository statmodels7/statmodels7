# modelterms7::scte() through statmod(): the shapes of the fitted surface, the
# surface graph against term_value(), the summary, and an equation for sigma.

test_that("a doubly increasing surface is fitted and drawn", {
  set.seed(16004)
  d <- data.frame(x = runif(300), z = runif(300))
  d$y <- pmin(d$x, 0.5) + 2 * pmin(d$z, 0.4) + rnorm(300, sd = 0.2)
  fit <- statmod(y ~ modelterms7::scte(x, z, shapes = list(modelterms7::shape_increasing(),
                                                           modelterms7::shape_increasing())),
                 distrib = distributions7::gaussian1_distrib(), data = d)
  gr <- expand.grid(x = seq(0.05, 0.95, length.out = 21), z = seq(0.05, 0.95, length.out = 21))
  mu <- matrix(predict(fit, what = "mu", newdata = gr), 21, 21)
  expect_true(all(diff(mu) >= -1e-10))
  expect_true(all(t(diff(t(mu))) >= -1e-10))
  expect_lt(sqrt(mean((mu - (pmin(gr$x, 0.5) + 2 * pmin(gr$z, 0.4)))^2)), 0.08)
  nm <- names(fit@spec@terms$mu)[2L]
  tm <- fit@spec@terms$mu[[nm]]
  expect_true(S7::S7_inherits(default_plot(tm), PlotSurface))
  pd <- plot_data(plot_surface(), tm, fit,
                  list(param = "mu", name = nm, vars = c("x", "z"), cols = NULL, data = d))
  tv <- modelterms7::term_value(tm, coef = coef(fit)$mu[modelterms7::term_coef_names(tm)],
                                newdata = data.frame(x = pd$surface$x, z = pd$surface$y))
  expect_equal(pd$surface$fit, tv, tolerance = 1e-12)
  expect_identical(term_block_kind(tm), "smooth")
})

test_that("a surface monotone in one covariate is fitted in sigma", {
  set.seed(16005)
  d <- data.frame(x = runif(300), z = runif(300))
  d$y <- rnorm(300, sd = exp(-1.5 + 1.2 * d$x + 0.3 * sin(2 * pi * d$z)))
  fit <- statmod(y ~ 1 | sigma ~ modelterms7::scte(x, z),
                 distrib = distributions7::gaussian1_distrib(), data = d)
  gr <- expand.grid(x = seq(0.05, 0.95, length.out = 21), z = seq(0.05, 0.95, length.out = 21))
  s <- matrix(predict(fit, what = "sigma", newdata = gr), 21, 21)
  expect_true(all(diff(s) >= -1e-10))
})

test_that("scte() under a factor or a numeric by is fitted and drawn", {
  set.seed(16006)
  d <- data.frame(x = runif(160), z = runif(160), g = gl(2, 80, labels = c("a", "b")),
                  w = rnorm(160))
  d$y <- c(a = 0, b = 0.5)[as.character(d$g)] + ifelse(d$g == "a", 1, 2) * pmin(d$x, 0.6) +
    sin(2 * pi * d$z) / 2 + d$w * (0.5 + d$x) + rnorm(160, sd = 0.2)
  sm <- basis7::bspline_smooth(k = 5)
  fit <- statmod(y ~ g + modelterms7::scte(x, z, smooths = sm, by = g) +
                   modelterms7::scte(x, z, smooths = sm, by = w),
                 distrib = distributions7::gaussian1_distrib(), data = d)
  nms <- names(fit@spec@terms$mu)
  tw <- fit@spec@terms$mu[[nms[3L]]]
  expect_identical(term_variables(tw, nms[3L]), c("x", "z", "w"))
  expect_identical(numeric_by_variable(tw), "w")
  pd <- plot_data(plot_surface(), tw, fit,
                  list(param = "mu", name = nms[3L], vars = c("x", "z", "w"), cols = NULL, data = d))
  tv <- modelterms7::term_value(tw, coef = coef(fit)$mu[modelterms7::term_coef_names(tw)],
                                newdata = data.frame(x = pd$surface$x, z = pd$surface$y, w = 1))
  expect_equal(pd$surface$fit, tv, tolerance = 1e-12)
  expect_match(attr(pd, "labels")$title, "at w = 1")
  # the coefficient of w is increasing in x
  gr <- data.frame(x = seq(0.05, 0.95, length.out = 31), z = 0.5, w = 1)
  expect_true(all(diff(modelterms7::term_value(tw, coef = coef(fit)$mu[modelterms7::term_coef_names(tw)],
                                               newdata = gr)) >= -1e-10))
})
