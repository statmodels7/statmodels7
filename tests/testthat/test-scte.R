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
