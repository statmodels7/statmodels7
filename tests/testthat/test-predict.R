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
  p <- predict(fit, "link", nd)
  expect_named(p, c("mu", "sigma"))
  expect_length(p$mu, 5L)
  expect_true(all(is.finite(unlist(p))))
})

test_that("a parameter can be asked for by its own name", {
  sim <- rstatmod(y ~ x | sigma ~ z, distributions7::gaussian1_distrib(), dd,
                  par = list(mu = c(1, 2), sigma = c(-1, 1.5)))
  fit <- statmod(y ~ x | sigma ~ z, distributions7::gaussian1_distrib(), sim)
  nd <- dd[1:8, ]
  th <- predict(fit, "parameter", nd)
  eta <- predict(fit, "link", nd)

  expect_equal(predict(fit, "mu", nd), th$mu, tolerance = 1e-14)
  expect_equal(predict(fit, "sigma", nd), th$sigma, tolerance = 1e-14)
  # the prefix asks for the predictor of that one parameter
  expect_equal(predict(fit, "link:sigma", nd), eta$sigma, tolerance = 1e-14)
  # the parameter is the inverse link of the predictor
  expect_equal(th$sigma, exp(eta$sigma), tolerance = 1e-12)

  # a parameter that does not vary still comes back one value per observation,
  # so a prediction has the length the data has whatever the model
  fit2 <- statmod(y ~ x, distributions7::gaussian1_distrib(), sim)
  expect_length(predict(fit2, "sigma", nd), 8L)
})

test_that("a moment can be asked for, and agrees with the parameters", {
  sim <- rstatmod(y ~ x | sigma ~ z, distributions7::gaussian1_distrib(), dd,
                  par = list(mu = c(1, 2), sigma = c(-1, 1.5)))
  fit <- statmod(y ~ x | sigma ~ z, distributions7::gaussian1_distrib(), sim)
  nd <- dd[1:8, ]
  th <- predict(fit, "parameter", nd)

  expect_equal(predict(fit, "mean", nd), th$mu, tolerance = 1e-12)
  expect_equal(predict(fit, "variance", nd), th$sigma^2, tolerance = 1e-10)
  expect_equal(predict(fit, "std_dev", nd), th$sigma, tolerance = 1e-10)
  # a gaussian's shape does not depend on the fit at all
  expect_equal(predict(fit, "skewness", nd), rep(0, 8), tolerance = 1e-12)
  expect_equal(predict(fit, "kurtosis", nd), rep(0, 8), tolerance = 1e-12)
})

test_that("a moment of a count model varies with the fitted mean", {
  sim <- rstatmod(counts ~ x, distributions7::poisson_distrib(), dd,
                  par = list(mu = c(0.5, 1.5)))
  fit <- statmod(counts ~ x, distributions7::poisson_distrib(), sim)
  nd <- dd[1:6, ]
  mu <- predict(fit, "mu", nd)
  # the Poisson's variance is its mean, which is a check of the routing
  expect_equal(predict(fit, "variance", nd), mu, tolerance = 1e-10)
  expect_equal(predict(fit, "skewness", nd), 1 / sqrt(mu), tolerance = 1e-8)
})

test_that("a level absent from the new data keeps its contrasts", {
  # the blueprint is what makes this work: rebuilding from the new frame would
  # give a different design and silently different predictions
  sim <- rstatmod(y ~ x + g, distributions7::gaussian1_distrib(), dd,
                  par = list(mu = c(1, 2, -0.5, 0.8), sigma = log(0.4)))
  fit <- statmod(y ~ x + g, distributions7::gaussian1_distrib(), sim)
  one_level <- dd[dd$g == "a", ][1:4, ]
  expect_length(levels(droplevels(one_level$g)), 1L)
  p <- predict(fit, "link", one_level)
  expect_length(p$mu, 4L)
  expect_true(all(is.finite(p$mu)))
})

test_that("a moment that does not exist is reported as NaN, not as a number", {
  # the Cauchy has no mean, and distributions7 answers NaN rather than
  # returning the location, which would be a plausible wrong number. The
  # prediction carries that through instead of hiding it.
  sim <- rstatmod(y ~ x, distributions7::cauchy_distrib(), dd,
                  par = list(mu = c(0, 1), sigma = log(0.5)))
  fit <- statmod(y ~ x, distributions7::cauchy_distrib(), sim)
  expect_true(all(is.nan(suppressWarnings(predict(fit, "mean", dd[1:3, ])))))
  expect_true(all(is.nan(suppressWarnings(
    predict(fit, "variance", dd[1:3, ])))))
  # while the parameters are perfectly available, which is the whole reason a
  # parameter can be asked for by name: a family with no moments still has a
  # location, and a fit still estimates it
  expect_true(all(is.finite(predict(fit, "mu", dd[1:3, ]))))
  expect_true(all(is.finite(predict(fit, "sigma", dd[1:3, ]))))
})

test_that("at the fitting data, prediction and fitted() are the same thing", {
  sim <- rstatmod(y ~ x, distributions7::gaussian1_distrib(), dd,
                  par = list(mu = c(1, 2), sigma = log(0.4)))
  fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), sim)
  expect_equal(predict(fit), fitted(fit), tolerance = 1e-12)
})

test_that("an unrecognized target names what is available", {
  sim <- rstatmod(y ~ x, distributions7::gaussian1_distrib(), dd,
                  par = list(mu = c(1, 2), sigma = log(0.4)))
  fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), sim)
  expect_error(predict(fit, "median"), "neither a parameter")
  expect_error(predict(fit, "median"), "mu, sigma")
  expect_error(predict(fit, "link:nu"), "neither a parameter")
  expect_error(predict(fit, c("mu", "sigma")), "a single string")
})

test_that("a data frame passed second is named rather than failing inside", {
  # the argument order departs from predict.lm on purpose, so the departure has
  # to report itself: without this a data frame lands in 'what' and the failure
  # surfaces several frames away from the mistake
  sim <- rstatmod(y ~ x, distributions7::gaussian1_distrib(), dd,
                  par = list(mu = c(1, 2), sigma = log(0.4)))
  fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), sim)
  expect_error(predict(fit, dd[1:3, ]), "is 'what', not")
  # and the spelling the message suggests works
  expect_length(predict(fit, newdata = dd[1:3, ])$mu, 3L)
})


test_that("new data reapplies each block instead of rebuilding it", {
  # A term records how its block was made -- a factor's levels and contrasts,
  # a spline's knots, a basis reparametrization -- and term_predict() replays
  # that record. Rebuilding gives a block of the same shape multiplying the
  # same coefficients and meaning something else, so the identity that has to
  # hold is that rows the model was fitted to predict to their fitted values.
  # Measured before the fix: 0.237 on 40 arbitrary rows and 1.19 on the 51
  # with |x| < 0.5, where rebuilt knots move furthest. The whole data handed
  # back agreed exactly, which is why nothing noticed.
  set.seed(5)
  n <- 200L
  dn <- data.frame(x = stats::runif(n, -2, 2),
                   g = factor(sample(letters[1:3], n, TRUE)))
  dn$y <- sin(1.4 * dn$x) + stats::rnorm(n, sd = 0.3)
  fit <- statmod(y ~ s(x, k = 10) + g, distributions7::gaussian1_distrib(), dn)
  full <- predict(fit, "mu")

  for (rows in list(1:40, which(abs(dn$x) < 0.5), c(3L, 7L, 100L))) {
    got <- predict(fit, "mu", dn[rows, , drop = FALSE])
    expect_equal(got, full[rows], tolerance = 1e-12)
  }
  # and the whole thing, which agreed even when it was wrong
  expect_equal(predict(fit, "mu", dn), full, tolerance = 1e-12)

  # the log-likelihood on a subset is the sum of that subset's contributions
  rows <- 1:60
  part <- as.numeric(loglik(fit, data = dn[rows, , drop = FALSE]))
  th <- predict(fit, "parameter", dn[rows, , drop = FALSE])
  by_hand <- sum(log(distributions7::distrib_pdf(
    fit@spec@distrib, dn$y[rows], th)))
  expect_equal(part, by_hand, tolerance = 1e-10)

  # a factor level absent from the new rows keeps its column, the levels
  # being the fitted term's and not the new data's
  drop_c <- droplevels(dn[dn$g != "c", , drop = FALSE])
  expect_equal(predict(fit, "mu", drop_c),
               full[dn$g != "c"], tolerance = 1e-12)
})
