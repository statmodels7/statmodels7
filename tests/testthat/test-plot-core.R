pt <- plot_test_fit()
fit <- pt$fit
dd <- pt$data

test_that("autoplot returns a ggplot and plot returns it invisibly", {
  p <- autoplot(fit, term = "s(x)", data = dd)
  expect_s3_class(p, "ggplot")
  pdf(NULL)
  on.exit(dev.off())
  expect_invisible(q <- plot(fit, term = "s(x)", data = dd))
  expect_s3_class(q, "ggplot")
})

test_that("the working residuals of a gaussian mean are the residuals", {
  wr <- working_residuals(fit, "mu", distributions7::expected_info())
  expect_equal(wr, fit@spec@response - fit@fitted$mu, tolerance = 1e-12)
  # log sigma: score (r^2/sigma^2 - 1), expected information 2
  r2 <- (fit@spec@response - fit@fitted$mu)^2 / fit@fitted$sigma^2
  ws <- working_residuals(fit, "sigma", distributions7::expected_info())
  expect_equal(ws, (r2 - 1) / 2, tolerance = 1e-12)
  wo <- working_residuals(fit, "sigma", distributions7::observed_info())
  expect_equal(wo, (r2 - 1) / (2 * r2), tolerance = 1e-9)
})

test_that("the effect of a smooth matches predict() on the same grid", {
  d <- plot_data(plot_effect(), fit@spec@terms$mu[["s(x)"]], fit,
                 ctx = list(param = "mu", name = "s(x)", vars = "x",
                            cols = NULL, data = dd))
  # the contribution of s(x) is the difference of two predictions that
  # differ only in x
  nd <- data.frame(x = d$curve$x, z = 0, g = factor("a", levels(dd$g)),
                   id = factor("1", levels(dd$id)))
  nd0 <- nd
  nd0$x <- d$curve$x[1]
  pr <- predict(fit, what = "mu", newdata = nd, random = "zero")
  pr0 <- predict(fit, what = "mu", newdata = nd0, random = "zero")
  expect_equal(d$curve$fit - d$curve$fit[1], pr - pr0, tolerance = 1e-10)
  # the partial residuals are the contribution plus the residual
  ofit <- term_contribution(fit, "mu", "s(x)", fit@spec@terms$mu[["s(x)"]],
                            NULL, NULL, bayesian_vcov())
  expect_equal(d$points$y, ofit$fit + fit@spec@response - fit@fitted$mu,
               tolerance = 1e-12)
  # the band is the delta method on the block of the term
  expect_true(all(d$curve$upper > d$curve$fit & d$curve$lower < d$curve$fit))
})

test_that("a variable selects its columns of the linear part", {
  p <- autoplot(fit, term = "z", data = dd)
  d <- p$layers[[2]]$data
  b <- coef(fit, readable = FALSE)$mu[["z"]]
  expect_equal(d$fit, b * d$x, tolerance = 1e-12)
  g <- autoplot(fit, term = "g", data = dd)
  gd <- g$layers[[2]]$data
  expect_equal(gd$fit, c(0, coef(fit, readable = FALSE)$mu[c("gb", "gc")]),
               ignore_attr = TRUE, tolerance = 1e-12)
})

test_that("the response scale is the fitted parameter", {
  p <- autoplot(fit, term = "x", param = "sigma",
                type = plot_effect(scale = "response"), data = dd)
  d <- p$layers[[1]]$data
  b <- coef(fit, readable = FALSE)$sigma
  expect_equal(d$fit, exp(b[[1]] + b[[2]] * d$x), tolerance = 1e-10)
})

test_that("binned means are the means of the partial residuals", {
  d <- plot_data(plot_effect(binned = 5), fit@spec@terms$sigma$linpar, fit,
                 ctx = list(param = "sigma", name = "linpar", vars = "x",
                            cols = 2L, data = dd))
  expect_equal(nrow(d$bins), 5)
  expect_equal(mean(d$bins$y), mean(d$points$y), tolerance = 0.2)
})

test_that("errors name what is missing", {
  expect_error(autoplot(fit, term = "s(x)"), "does not keep its data")
  expect_error(autoplot(fit, term = "nope", data = dd), "neither a term")
  expect_error(autoplot(fit, term = "x", data = dd), "give 'param'")
  expect_error(autoplot(fit, term = "s(x)", data = dd[1:10, ]),
               "10 rows")
  expect_error(autoplot(fit, term = "s(x)", type = 3), "plot_type object")
  expect_error(autoplot(fit, param = "mu"), "give 'term'")
  expect_error(plot_effect(level = 2), "'level'")
  expect_error(plot_effect(scale = "other"))
  expect_error(plot_effect(information = 1), "expected_info")
})

test_that("a term name ignores spaces", {
  p <- autoplot(fit, term = "random(~1|id)", data = dd)
  expect_s3_class(p, "ggplot")
})

test_that("a variable in a linear and a nonlinear term sums their values", {
  set.seed(4)
  dn <- data.frame(x = runif(150, 0, 3))
  dn$y <- 0.5 * dn$x + 5 * exp(-0.8 * dn$x) + rnorm(150, sd = 0.2)
  fn <- statmod(y ~ x + nl(~ a * exp(-r * x), start = list(a = 4, r = 1),
                           links = list(r = linkfunctions7::log_link())),
                distributions7::gaussian1_distrib(), dn)
  p <- autoplot(fn, term = "x", data = dn)
  d <- p$layers[[2]]$data
  b <- coef(fn, readable = FALSE)$mu
  expect_equal(d$fit, b[["x"]] * d$x + b[["nl.a"]] * exp(-exp(b[["nl.r"]]) * d$x),
               tolerance = 1e-10)
})

test_that("a factor is crossed only where it interacts with the variable", {
  p <- autoplot(fit, term = "z", data = dd)
  expect_false("group" %in% names(p$layers[[2]]$data))
  set.seed(3)
  di <- data.frame(z = rnorm(120), g = factor(sample(c("a", "b"), 120, TRUE)))
  di$y <- di$z * c(1, -1)[di$g] + rnorm(120, sd = 0.3)
  fi <- statmod(y ~ z * g, distributions7::gaussian1_distrib(), di)
  q <- autoplot(fi, term = "z", data = di)
  d <- q$layers[[2]]$data
  expect_setequal(levels(d$group), c("a", "b"))
  b <- coef(fi, readable = FALSE)$mu
  db <- d[d$group == "b", ]
  expect_equal(db$fit, (b[["z"]] + b[["z:gb"]]) * db$x + b[["gb"]] * 0,
               tolerance = 1e-10)
})
