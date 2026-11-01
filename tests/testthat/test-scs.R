# A fit with modelterms7::scs(): the shape on a grid in each parameter, and
# the summary block of a smooth (the smoothing parameter, not the free
# coefficients).

scs_data <- function() {
  set.seed(31)
  x <- runif(200)
  f <- ifelse(x < 0.6, 1 - (1 - x / 0.6)^2, 1)
  data.frame(x = x, y = f + exp(-2.3 + 1.2 * f) * rnorm(200))
}

test_that("a fit with scs() in mu and sigma is monotone in both", {
  d <- scs_data()
  fit <- statmod(y ~ modelterms7::scs(x) | sigma ~ modelterms7::scs(x),
                 distrib = distributions7::gaussian1_distrib(), data = d)
  nd <- data.frame(x = seq(0, 1, length.out = 401))
  mu <- predict(fit, what = "mu", newdata = nd)
  sg <- predict(fit, what = "sigma", newdata = nd)
  expect_true(all(diff(mu) >= -1e-10))
  expect_true(all(diff(sg) >= -1e-10))
  expect_true(is.finite(as.numeric(logLik(fit))))
})

test_that("a convex scs() gives a convex fitted mean", {
  set.seed(32)
  d <- data.frame(x = runif(150))
  d$y <- (d$x - 0.4)^2 * 3 + rnorm(150, sd = 0.1)
  fit <- statmod(y ~ modelterms7::scs(x, modelterms7::shape_convex()),
                 distrib = distributions7::gaussian1_distrib(), data = d)
  mu <- predict(fit, what = "mu", newdata = data.frame(x = seq(0, 1, length.out = 401)))
  expect_true(all(diff(diff(mu)) >= -1e-10))
})

test_that("scs() is summarized as a smooth", {
  d <- scs_data()
  fit <- statmod(y ~ modelterms7::scs(x),
                 distrib = distributions7::gaussian1_distrib(), data = d)
  tm <- fit@spec@terms$mu[["modelterms7::scs(x)"]]
  if (is.null(tm)) tm <- Filter(function(t) S7::S7_inherits(t, modelterms7::ScsTerm),
                                fit@spec@terms$mu)[[1L]]
  expect_identical(term_block_kind(tm), "smooth")
  out <- paste(capture.output(print(summary(fit))), collapse = "\n")
  expect_match(out, "lambda [reml]", fixed = TRUE)
  expect_false(grepl("\\.g1", out))
})

test_that("a value at a point holds in a fit without an intercept and is inert with one", {
  set.seed(35)
  d <- data.frame(x = runif(150, 0.5, 3))
  d$y <- 2 * (1 - exp(-1.2 * d$x)) + rnorm(150, sd = 0.15)
  sm <- basis7::bspline_smooth(k = 8, lower = 0, upper = 3)
  g <- distributions7::gaussian1_distrib()
  f0 <- statmod(y ~ 0 + modelterms7::scs(
    x, modelterms7::shape_increasing_concave(value = c(at = 0, value = 0)), sm),
    distrib = g, data = d)
  expect_equal(predict(f0, what = "mu", newdata = data.frame(x = 0)), 0)
  mu <- predict(f0, what = "mu", newdata = data.frame(x = seq(0, 3, length.out = 301)))
  expect_true(all(diff(mu) >= -1e-10))
  expect_true(all(diff(diff(mu)) <= 1e-10))
  fits <- lapply(c(0, 5), function(v) statmod(
    y ~ modelterms7::scs(x, modelterms7::shape_increasing(value = c(at = 1, value = v)), sm),
    distrib = g, data = d))
  fc <- statmod(y ~ modelterms7::scs(x, modelterms7::shape_increasing(), sm),
                distrib = g, data = d)
  nd <- data.frame(x = c(0.7, 1.5, 2.8))
  for (f in fits) {
    expect_equal(as.numeric(logLik(f)), as.numeric(logLik(fc)), tolerance = 1e-8)
    expect_equal(predict(f, what = "mu", newdata = nd),
                 predict(fc, what = "mu", newdata = nd), tolerance = 1e-8)
  }
  expect_equal(coef(fits[[1L]])$mu[["(Intercept)"]] - coef(fits[[2L]])$mu[["(Intercept)"]],
               5, tolerance = 1e-8)
})

test_that("a slope at a point holds in a fitted convex mean", {
  set.seed(36)
  d <- data.frame(x = runif(150))
  d$y <- exp(2 * d$x) - 2 * d$x + rnorm(150, sd = 0.2)
  fit <- statmod(y ~ modelterms7::scs(
    x, modelterms7::shape_convex(slope = c(at = 0.5, value = 0.718))),
    distrib = distributions7::gaussian1_distrib(), data = d)
  h <- 1e-4
  mu <- function(x) predict(fit, what = "mu", newdata = data.frame(x = x))
  expect_equal((mu(0.5 + h) - mu(0.5 - h)) / (2 * h), 0.718, tolerance = 1e-6)
  grid <- seq(min(d$x), max(d$x), length.out = 401)
  expect_true(all(diff(diff(mu(grid))) >= -1e-10))
})

test_that("a factor by fits one monotone curve per level", {
  set.seed(37)
  d <- data.frame(x = runif(180), g = gl(3, 60, labels = c("a", "b", "c")))
  d$y <- c(a = 0, b = 0.5, c = 1)[as.character(d$g)] +
    ifelse(d$g == "a", 2 * d$x, ifelse(d$g == "b", d$x^2, plogis((d$x - 0.5) / 0.1))) +
    rnorm(180, sd = 0.12)
  sm <- basis7::bspline_smooth(k = 7)
  fit <- statmod(y ~ g + modelterms7::scs(x, modelterms7::shape_increasing(), sm, by = g),
                 distrib = distributions7::gaussian1_distrib(), data = d)
  expect_identical(statmod_certificate(fit)$state, "converged")
  grid <- expand.grid(x = seq(min(d$x), max(d$x), length.out = 201), g = levels(d$g))
  mu <- predict(fit, what = "mu", newdata = grid)
  expect_true(all(tapply(mu, grid$g, function(v) min(diff(v))) >= -1e-10))
  expect_equal(predict(fit, what = "mu", newdata = d), predict(fit, what = "mu"),
               tolerance = 1e-10)
})

test_that("a numeric by equal to 2 fits as the smooth without by", {
  set.seed(38)
  d <- data.frame(x = runif(120), z = 2)
  d$y <- exp(2 * d$x) / 3 + rnorm(120, sd = 0.2)
  sm <- basis7::bspline_smooth(k = 7)
  g <- distributions7::gaussian1_distrib()
  fz <- statmod(y ~ 0 + modelterms7::scs(x, modelterms7::shape_increasing(), sm, by = z),
                distrib = g, data = d)
  f1 <- statmod(y ~ modelterms7::scs(x, modelterms7::shape_increasing(), sm),
                distrib = g, data = d)
  expect_equal(as.numeric(logLik(fz)), as.numeric(logLik(f1)), tolerance = 1e-8)
  expect_equal(predict(fz, what = "mu"), predict(f1, what = "mu"), tolerance = 1e-8)
  # the free level enters z f(x) scaled by 2: the REML constant is -log 2
  expect_equal(fz@criterion - f1@criterion, -log(2), tolerance = 1e-6)
})

test_that("the modes of a factor by are searched one level at a time", {
  set.seed(39)
  d <- data.frame(x = runif(160), g = gl(2, 80, labels = c("a", "b")))
  d$y <- exp(-((d$x - c(a = 0.3, b = 0.7)[as.character(d$g)]) / 0.2)^2) +
    rnorm(160, sd = 0.1)
  sm <- basis7::bspline_smooth(k = 6)
  fit <- statmod(y ~ g + modelterms7::scs(x, modelterms7::shape_unimodal(), sm, by = g),
                 distrib = distributions7::gaussian1_distrib(), data = d)
  va <- fit@history$variants
  lev <- sub(":.*", "", va$version)
  expect_identical(unique(lev), c("a", "b"))
  expect_identical(sum(va$chosen), 1L)
  expect_identical(va$score[va$chosen], min(va$score))
  tm <- Filter(function(t) S7::S7_inherits(t, modelterms7::ScsTerm), fit@spec@terms$mu)[[1L]]
  idx <- tm@spec$index_by
  expect_identical(names(idx), c("a", "b"))
  # the chosen version carries the indices of the term that was fitted
  best <- va$version[va$chosen]
  expect_identical(idx[[sub(":.*", "", best)]], as.integer(sub(".*M = ", "", best)))
  mu <- predict(fit, what = "mu", newdata = data.frame(x = seq(0, 1, length.out = 201),
                                                       g = "b"))
  expect_gt(seq(0, 1, length.out = 201)[which.max(mu)], 0.5)
})

test_that("the graphs of a numeric by draw the coefficient at z = 1", {
  set.seed(40)
  d <- data.frame(x = runif(150), z = rnorm(150))
  d$y <- 1 + d$z * (0.5 + 1.5 * d$x^2) + rnorm(150, sd = 0.2)
  sm <- basis7::bspline_smooth(k = 7)
  g <- distributions7::gaussian1_distrib()
  for (fo in list(y ~ modelterms7::scs(x, modelterms7::shape_increasing(), sm, by = z),
                  y ~ s(x, smoother = sm, by = z))) {
    fit <- statmod(fo, distrib = g, data = d)
    nm <- names(fit@spec@terms$mu)[2L]
    tm <- fit@spec@terms$mu[[nm]]
    expect_identical(term_variables(tm, nm), c("x", "z"))
    expect_identical(numeric_by_variable(tm), "z")
    pd <- plot_data(plot_effect(), tm, fit,
                    list(param = "mu", name = nm, vars = term_variables(tm, nm),
                         cols = NULL, data = d))
    cf <- coef(fit)$mu[modelterms7::term_coef_names(tm)]
    ref <- modelterms7::term_value(tm, coef = cf,
                                   newdata = data.frame(x = pd$curve$x, z = 1))
    expect_equal(pd$curve$fit, ref, tolerance = 1e-12)
    expect_identical(nrow(pd$points), 0L)
    expect_match(attr(pd, "labels")$title, "at z = 1")
    expect_s3_class(plot(fit, plot_type = plot_derivative(), term = nm, data = d), "ggplot")
  }
})

test_that("the warnings of a discarded version stay in the history", {
  set.seed(39)
  d <- data.frame(x = runif(80))
  d$y <- exp(-((d$x - 0.3) / 0.2)^2) + rnorm(80, sd = 0.1)
  fit <- expect_no_warning(statmod(y ~ modelterms7::scs(x, modelterms7::shape_unimodal(),
                                                        basis7::bspline_smooth(k = 6)),
                                   distrib = distributions7::gaussian1_distrib(), data = d))
  va <- fit@history$variants
  expect_true("warning" %in% names(va))
  expect_identical(va$warning[va$chosen], "")
})

test_that("the settings of a shape may be shared or differ across the levels", {
  set.seed(41)
  d <- data.frame(x = runif(200), g = gl(2, 100, labels = c("a", "b")))
  d$y <- exp(-((d$x - c(a = 0.35, b = 0.6)[as.character(d$g)]) / 0.2)^2) + rnorm(200, sd = 0.1)
  sm <- basis7::bspline_smooth(k = 7)
  g0 <- distributions7::gaussian1_distrib()
  fs <- statmod(y ~ g + modelterms7::scs(x, modelterms7::shape_unimodal(mode_by = "shared"), sm,
                                         by = g), distrib = g0, data = d)
  va <- fs@history$variants
  expect_true(all(grepl("^M = ", va$version)))
  tm <- Filter(function(t) S7::S7_inherits(t, modelterms7::ScsTerm), fs@spec@terms$mu)[[1L]]
  expect_identical(length(tm@build_info$pieces), 1L)
  # a convex mean with one left slope for both levels
  d$y2 <- exp(1.5 * d$x) + c(a = 0, b = 1)[as.character(d$g)] * d$x^2 + rnorm(200, sd = 0.1)
  ft <- statmod(y2 ~ g + modelterms7::scs(x, modelterms7::shape_convex(left_slope_by = "shared"),
                                          sm, by = g), distrib = g0, data = d)
  expect_identical(statmod_certificate(ft)$state, "converged")
  grid <- expand.grid(x = seq(min(d$x), max(d$x), length.out = 101), g = c("a", "b"))
  mu <- predict(ft, what = "mu", newdata = grid)
  expect_true(all(tapply(mu, grid$g, function(v) min(diff(diff(v)))) >= -1e-10))
  # minimum slopes by level
  d$y3 <- c(a = 0.5, b = 2)[as.character(d$g)] * d$x + rnorm(200, sd = 0.1)
  fm <- statmod(y3 ~ g + modelterms7::scs(x, modelterms7::shape_increasing(min_slope = c(a = 0.2, b = 1)),
                                          sm, by = g), distrib = g0, data = d)
  sl <- tapply(predict(fm, what = "mu", newdata = grid), grid$g, function(v) min(diff(v)) / diff(grid$x[1:2]))
  expect_true(sl[["a"]] >= 0.2 - 1e-8 && sl[["b"]] >= 1 - 1e-8)
})

test_that("an increasing shape in sigma on a constant truth fits", {
  set.seed(1)
  d <- data.frame(x = runif(240))
  d$y <- sin(3 * d$x) + rnorm(240, sd = 0.15)
  sm <- basis7::bspline_smooth(k = 8)
  fit <- suppressWarnings(statmod(y ~ s(x, smoother = sm) |
                                    sigma ~ modelterms7::scs(x, modelterms7::shape_increasing(), sm),
                                  distrib = distributions7::gaussian1_distrib(), data = d))
  s <- predict(fit, what = "sigma", newdata = data.frame(x = c(0.1, 0.5, 0.9)))
  expect_true(all(diff(s) >= -1e-10))
  expect_lt(abs(log(s[3] / 0.15)), 0.2)
})
