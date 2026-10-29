# The object that statmod_profile() returns, its intervals and its graphs.
# The reference of the intervals is the inversion of the likelihood ratio
# test in confint.StatmodFit(), an independent route to the same limits.

lin_fit <- function() {
  set.seed(1)
  dd <- data.frame(x = runif(80))
  dd$y <- 1 + 2 * dd$x + rnorm(80, sd = 0.4)
  statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
}

smooth_fit <- function() {
  set.seed(3)
  d <- data.frame(z = runif(200), x = rnorm(200))
  d$y <- rpois(200, exp(1 + sin(2 * pi * d$z) + 0.4 * d$x))
  statmod(y ~ s(z) + x, distributions7::poisson_distrib(), d)
}

test_that("the object records the labels, their kinds and the profile column", {
  fit <- smooth_fit()
  p <- statmod_profile(fit, c("mu/s(z)/lambda", "mu:s(z).z3", "mu:x"),
                       values = profile_range(width = 1, n = 3), joint = FALSE)
  expect_true(S7::S7_inherits(p, StatmodProfile))
  expect_identical(p@profiled_labels, c("mu/s(z)/lambda", "mu:s(z).z3", "mu:x"))
  q <- p@profiled_quantities
  expect_true(S7::S7_inherits(q[["mu/s(z)/lambda"]], ProfiledHyperparameter))
  expect_true(S7::S7_inherits(q[["mu:s(z).z3"]], ProfiledPenalizedCoefficient))
  expect_true(S7::S7_inherits(q[["mu:x"]], ProfiledCoefficient))
  expect_identical(unname(p@profile_column),
                   c("outer_criterion", "profile_pen_loglik", "profile_pen_loglik"))
  expect_true(all(p@is_likelihood_ratio))
  expect_false(p@joint)
  expect_identical(names(p@points)[1:2], c("profiled_label", "value"))
  expect_identical(nrow(p@points), 9L)
  expect_identical(p@unrestricted_point$profiled_label, p@profiled_labels)
  expect_identical(as.data.frame(p), p@points)
  expect_output(print(p), "penalized coefficient")
  expect_identical(p@hyper_handling, "held")
  expect_identical(p@outer_method@kind, "reml")
})

test_that("a separate profile equals the profile of each label by itself", {
  fit <- lin_fit()
  v <- profile_range(width = 2, n = 5)
  p <- statmod_profile(fit, c("mu:x", "mu:(Intercept)"), values = v, joint = FALSE)
  for (lab in p@profiled_labels) {
    one <- statmod_profile(fit, lab, values = v)
    rows <- p@points$profiled_label == lab
    expect_equal(p@points$value[rows], one@points[[lab]])
    expect_equal(p@points$profile_pen_loglik[rows], one@points$profile_pen_loglik,
                 tolerance = 1e-10)
  }
  expect_error(statmod_profile(fit, c("mu:x", "mu:(Intercept)"),
                               values = list("mu:z" = 1:3), joint = FALSE),
               "not in 'which'")
})

test_that("confint() inverts the likelihood ratio like confint() of the fit", {
  fit <- lin_fit()
  p <- statmod_profile(fit, "mu:x", values = profile_range(width = 3, n = 13))
  ci <- confint(p)
  ref <- confint(fit, parm = "mu:x", test = lr_test(), readable = FALSE)
  expect_equal(ci$lower, ref$lower, tolerance = 1e-6)
  expect_equal(ci$upper, ref$upper, tolerance = 1e-6)
  expect_identical(rownames(ci), "mu:x")
  ci9 <- confint(p, level = 0.9)
  ref9 <- confint(fit, parm = "mu:x", level = 0.9, test = lr_test(), readable = FALSE)
  expect_equal(c(ci9$lower, ci9$upper), c(ref9$lower, ref9$upper), tolerance = 1e-6)
})

test_that("simultaneous intervals use the chi-squared quantile with k degrees of freedom", {
  fit <- lin_fit()
  p <- statmod_profile(fit, c("mu:x", "mu:(Intercept)"),
                       values = profile_range(width = 4, n = 17), joint = FALSE)
  s <- confint(p, type = "simultaneous")
  lev <- stats::pchisq(stats::qchisq(0.95, 2), 1)
  ref <- confint(fit, parm = c("mu:x", "mu:(Intercept)"), level = lev,
                 test = lr_test(), readable = FALSE)
  expect_equal(s[rownames(ref), "lower"], ref$lower, tolerance = 1e-5)
  expect_equal(s[rownames(ref), "upper"], ref$upper, tolerance = 1e-5)
  m <- confint(p)
  expect_true(all(s$lower < m$lower & s$upper > m$upper))
})

test_that("a joint profile gives intervals from the maximum over the other labels", {
  fit <- lin_fit()
  j <- statmod_profile(fit, c("mu:x", "mu:(Intercept)"),
                       values = profile_range(width = 3, n = 9))
  ci <- confint(j)
  ref <- confint(fit, parm = c("mu:x", "mu:(Intercept)"), test = lr_test(),
                 readable = FALSE)
  # the grid of the other label resolves the maximum to a few per cent of the
  # width of the interval
  w <- ref$upper - ref$lower
  expect_lt(max(abs(ci[rownames(ref), "lower"] - ref$lower) / w), 0.05)
  expect_lt(max(abs(ci[rownames(ref), "upper"] - ref$upper) / w), 0.05)
})

test_that("a limit outside the values is missing, with a warning", {
  fit <- lin_fit()
  p <- statmod_profile(fit, "mu:x", values = profile_range(width = 0.5, n = 5))
  expect_warning(ci <- confint(p), "Widen the values")
  expect_true(is.na(ci$lower) && is.na(ci$upper))
  expect_error(confint(p, parm = "mu:z"), "does not hold")
  expect_error(confint(p, level = 1.2), "between 0 and 1")
})

test_that("confint() rejects a profile whose drop is not a likelihood ratio", {
  fit <- smooth_fit()
  p <- statmod_profile(fit, c("mu/s(z)/lambda", "mu:x"),
                       values = list("mu/s(z)/lambda" = c(1, 1.3),
                                     "mu:x" = c(0.35, 0.4)))
  expect_false(any(p@is_likelihood_ratio))
  expect_identical(unname(p@profile_column), rep("outer_criterion", 2))
  expect_error(confint(p), "not a likelihood ratio")
  expect_output(print(p), "not a likelihood ratio")
})

test_that("each kind of label is drawn by its own graph", {
  fit <- smooth_fit()
  p <- statmod_profile(fit, c("mu/s(z)/lambda", "mu:x"),
                       values = profile_range(width = 1, n = 3), joint = FALSE)
  g <- autoplot(p)
  expect_s3_class(g, "patchwork")
  one <- statmod_profile(fit, "mu:x", values = profile_range(width = 1, n = 3))
  g1 <- autoplot(one)
  expect_s3_class(g1, "ggplot")
  expect_identical(g1$labels$subtitle, "coefficient")
  gh <- autoplot(statmod_profile(fit, "mu/s(z)/lambda",
                                 values = profile_range(width = 1, n = 3)))
  expect_identical(gh$labels$subtitle, "hyperparameter")
  set.seed(5)
  w <- runif(120, 0, 10)
  db <- data.frame(w = w, y = 1 + 0.8 * pmax(w - 4, 0) + rnorm(120))
  fb <- statmod(y ~ seg(w), distributions7::gaussian1_distrib(), db)
  pb <- statmod_profile(fb, "mu:seg.psi1", values = profile_range(width = 1, n = 3))
  expect_true(S7::S7_inherits(pb@profiled_quantities[[1L]], ProfiledBreakpoint))
  expect_length(pb@profiled_quantities[[1L]]@covariate_values, 120L)
  gb <- autoplot(pb)
  expect_identical(gb$labels$subtitle, "break-point position")
  expect_true(any(vapply(gb$layers, function(l) inherits(l$geom, "GeomRug"), TRUE)))
})

test_that("a joint profile is drawn as contours, one panel for each pair", {
  fit <- lin_fit()
  v <- profile_range(width = 2, n = 3)
  j2 <- statmod_profile(fit, c("mu:x", "mu:(Intercept)"), values = v)
  g2 <- autoplot(j2)
  expect_s3_class(g2, "ggplot")
  expect_true(any(vapply(g2$layers, function(l) inherits(l$geom, "GeomContour"), TRUE)))
  j3 <- statmod_profile(fit, c("mu:x", "mu:(Intercept)", "sigma:(Intercept)"),
                        values = v)
  g3 <- autoplot(j3)
  expect_s3_class(g3, "patchwork")
  expect_length(g3$patches$plots, 2L)
})
