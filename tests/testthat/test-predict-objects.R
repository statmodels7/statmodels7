# The interval and predictive objects of predict().

pobj_fit <- function(seed = 2) {
  set.seed(seed)
  dd <- data.frame(x = runif(30))
  dd$y <- 1 + dd$x + rnorm(30, sd = 0.5)
  statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
}

test_that("the constructors build objects of the right classes", {
  expect_true(S7::S7_inherits(confidence_interval(), interval_type))
  expect_true(S7::S7_inherits(group_interval(), interval_type))
  pi <- prediction_interval()
  expect_true(S7::S7_inherits(pi, interval_type))
  expect_true(S7::S7_inherits(pi@predictive, PredAveraged))
  expect_true(S7::S7_inherits(prediction_interval(pred_plugin())@predictive,
                              PredPlugin))
  b <- pred_bootstrap()
  expect_identical(b@n, 200L)
  expect_identical(b@refit, "coefficients")
  b2 <- pred_bootstrap(n = 50, refit = "full")
  expect_identical(b2@n, 50L)
  expect_identical(b2@refit, "full")
  expect_true(S7::S7_inherits(b2, predictive_type))
})

test_that("the abstract parents cannot be built", {
  expect_error(interval_type("x"))
  expect_error(predictive_type("x"))
})

test_that("invalid arguments are rejected", {
  expect_error(pred_bootstrap(n = 0), "positive whole number")
  expect_error(pred_bootstrap(n = 2.5), "positive whole number")
  expect_error(pred_bootstrap(n = c(1, 2)), "positive whole number")
  expect_error(pred_bootstrap(refit = "partial"))
  expect_error(prediction_interval(predictive = "plugin"),
               "no longer accepted")
  expect_error(prediction_interval(predictive = 1),
               "must be a predictive object")
})

test_that("a string for interval signals the error naming the object", {
  fit <- pobj_fit()
  for (s in c("confidence", "group", "prediction")) {
    expect_error(predict(fit, what = "mu", interval = s),
                 "must be an interval object such as confidence_interval[(][)]")
    expect_error(predict(fit, what = "mu", interval = s), "no longer accepted")
  }
  expect_error(predict(fit, what = "mu", interval = 1),
               "must be an interval object")
})

test_that("the print methods describe the object", {
  expect_output(print(confidence_interval()), "confidence")
  expect_output(print(prediction_interval(pred_bootstrap(n = 7))),
                "bootstrap, 7 replicas")
  expect_output(print(pred_plugin()), "plugin")
})

test_that("the default interval is the confidence interval", {
  fit <- pobj_fit()
  nd <- data.frame(x = c(0.2, 0.7))
  expect_identical(
    predict(fit, what = "mu", newdata = nd, se = TRUE),
    predict(fit, what = "mu", newdata = nd, se = TRUE,
            interval = confidence_interval()))
})

test_that("the default prediction interval is the averaged one", {
  fit <- pobj_fit()
  nd <- data.frame(x = 0.5)
  expect_identical(
    predict(fit, what = "response", newdata = nd,
            interval = prediction_interval()),
    predict(fit, what = "response", newdata = nd,
            interval = prediction_interval(predictive = pred_averaged())))
})

test_that("the bootstrap object reproduces the seeded replicas", {
  fit <- pobj_fit()
  nd <- data.frame(x = 0.5)
  set.seed(11)
  a <- predict(fit, what = "response", newdata = nd,
               interval = prediction_interval(pred_bootstrap(n = 12)))
  set.seed(11)
  b <- predict(fit, what = "response", newdata = nd,
               interval = prediction_interval(
                 pred_bootstrap(n = 12, refit = "coefficients")))
  expect_identical(a, b)
})
