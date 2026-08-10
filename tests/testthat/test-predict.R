# Prediction, and simulating a response from a written model.

set.seed(2)
n <- 200
dd <- data.frame(x = runif(n), z = runif(n),
                 g = factor(rep(c("a", "b", "c"), length.out = n)))

test_that("a written model draws a response and a fit recovers it", {
  sim <- rstatmod(y ~ x + g, distributions7::gaussian1_distrib(), dd,
                  par = list(mu = c(1, 2, -0.5, 0.8), sigma = log(0.4)))
  expect_true("y" %in% names(sim))
  expect_identical(nrow(sim), 200L)
  fit <- statmod(y ~ x + g, distributions7::gaussian1_distrib(), sim)
  # 200 observations and a scale of 0.4: the estimate is within a few standard
  # errors of the truth, and the tolerance says so rather than pretending to
  # be exact
  expect_equal(fit@coefficients$mu, c(1, 2, -0.5, 0.8), tolerance = 0.2)
  expect_equal(exp(fit@coefficients$sigma), 0.4, tolerance = 0.1)
})

test_that("a factor needs no declaring", {
  # the design comes from the same interpreter a fit uses, so contrasts
  # appear by themselves
  sim <- rstatmod(y ~ g, distributions7::gaussian1_distrib(), dd,
                  par = list(mu = c(0, 1, 2), sigma = log(0.3)))
  expect_identical(names(attr(sim, "par")$mu),
                   c("(Intercept)", "gb", "gc"))
  fit <- statmod(y ~ g, distributions7::gaussian1_distrib(), sim)
  expect_equal(fit@coefficients$mu, c(0, 1, 2), tolerance = 0.15)
})

test_that("a count model is drawn and recovered", {
  sim <- rstatmod(counts ~ x, distributions7::poisson_distrib(), dd,
                  par = list(mu = c(0.5, 1.5)))
  expect_true(all(sim$counts >= 0))
  expect_true(all(sim$counts == round(sim$counts)))
  fit <- statmod(counts ~ x, distributions7::poisson_distrib(), sim)
  expect_equal(fit@coefficients$mu, c(0.5, 1.5), tolerance = 0.25)
})

test_that("every parameter can be simulated and recovered", {
  sim <- rstatmod(y ~ x | sigma ~ z, distributions7::gaussian1_distrib(), dd,
                  par = list(mu = c(1, 2), sigma = c(-1, 1.5)))
  fit <- statmod(y ~ x | sigma ~ z, distributions7::gaussian1_distrib(), sim)
  expect_equal(fit@coefficients$mu, c(1, 2), tolerance = 0.2)
  expect_equal(fit@coefficients$sigma, c(-1, 1.5), tolerance = 0.5)
})

test_that("the coefficients are drawn when none are given, and reported", {
  set.seed(11)
  sim <- rstatmod(y ~ x + g, distributions7::gaussian1_distrib(), dd)
  p <- attr(sim, "par")
  expect_named(p, c("mu", "sigma"))
  expect_length(p$mu, 4L)
  expect_true(all(is.finite(unlist(p))))
  # what was used is reported, so the simulation is reproducible from its own
  # output
  again <- rstatmod(y ~ x + g, distributions7::gaussian1_distrib(), dd,
                    par = p)
  expect_equal(attr(again, "theta")$mu, attr(sim, "theta")$mu)
})

test_that("what rstatmod refuses", {
  expect_error(rstatmod(y ~ x, distributions7::gaussian1_distrib(),
                        list(x = 1)), "must be a data frame")
  expect_error(rstatmod(y ~ x, "gaussian", dd), "distributions7")
  expect_error(rstatmod(y ~ x, distributions7::gaussian1_distrib(), dd,
                        par = list(wrong = 1)), "not a parameter")
  expect_error(rstatmod(y ~ x, distributions7::gaussian1_distrib(), dd,
                        par = list(mu = 1)), "length 1")
})

test_that("prediction does not need the response", {
  sim <- rstatmod(y ~ x + g, distributions7::gaussian1_distrib(), dd,
                  par = list(mu = c(1, 2, -0.5, 0.8), sigma = log(0.4)))
  fit <- statmod(y ~ x + g, distributions7::gaussian1_distrib(), sim)
  # new data routinely has no response column, and a prediction that demanded
  # one would be useless for the thing prediction is for
  nd <- dd[1:5, ]
  expect_false("y" %in% names(nd))
  p <- predict(fit, nd, "link")
  expect_named(p, c("mu", "sigma"))
  expect_length(p$mu, 5L)
  expect_true(all(is.finite(unlist(p))))
})

test_that("the three types of prediction agree with each other", {
  sim <- rstatmod(y ~ x, distributions7::gaussian1_distrib(), dd,
                  par = list(mu = c(1, 2), sigma = log(0.4)))
  fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), sim)
  nd <- dd[1:8, ]
  eta <- predict(fit, nd, "link")
  th <- predict(fit, nd, "parameter")
  # the parameter is the inverse link of the predictor
  expect_equal(th$sigma, exp(eta$sigma), tolerance = 1e-12)
  expect_equal(th$mu, eta$mu, tolerance = 1e-12)
  # and for a gaussian the mean IS the location
  expect_equal(predict(fit, nd, "response"), th$mu, tolerance = 1e-12)
  # at the fitting data, prediction and fitted() are the same thing
  expect_equal(predict(fit), fitted(fit), tolerance = 1e-12)
})

test_that("a level absent from the new data keeps its contrasts", {
  # the blueprint is what makes this work: rebuilding from the new frame would
  # give a different design and silently different predictions
  sim <- rstatmod(y ~ x + g, distributions7::gaussian1_distrib(), dd,
                  par = list(mu = c(1, 2, -0.5, 0.8), sigma = log(0.4)))
  fit <- statmod(y ~ x + g, distributions7::gaussian1_distrib(), sim)
  one_level <- dd[dd$g == "a", ][1:4, ]
  expect_length(levels(droplevels(one_level$g)), 1L)
  p <- predict(fit, one_level, "link")
  expect_length(p$mu, 4L)
  expect_true(all(is.finite(p$mu)))
})

test_that("a mean that does not exist is reported as NaN, not as a number", {
  # the Cauchy has no mean, and distributions7 answers NaN rather than
  # returning the location, which would be a plausible wrong number. The
  # prediction carries that through instead of hiding it.
  sim <- rstatmod(y ~ x, distributions7::cauchy_distrib(), dd,
                  par = list(mu = c(0, 1), sigma = log(0.5)))
  fit <- statmod(y ~ x, distributions7::cauchy_distrib(), sim)
  expect_true(all(is.nan(predict(fit, dd[1:3, ], "response"))))
  # while the parameters are perfectly available, which is what to ask for
  expect_true(all(is.finite(predict(fit, dd[1:3, ], "parameter")$mu)))
})
