## The test statistic and the variance convention as objects.
##
## `test` and `vcov_type` used to be character strings. Both are now S7 objects, a
## string is rejected with an error that names the constructor, and the
## default object has to give the numbers that the default string gave.

small_fit <- function() {
  set.seed(1)
  d <- data.frame(x = stats::runif(60))
  d$y <- 1 + 2 * d$x + stats::rnorm(60, sd = 0.4)
  statmod(y ~ x, distributions7::gaussian1_distrib(), d)
}

test_that("the four test constructors build objects of the right classes", {
  tests <- list(wald = wald_test(), lr = lr_test(), score = score_test(),
                gradient = gradient_test())
  for (nm in names(tests)) {
    tt <- tests[[nm]]
    expect_true(S7::S7_inherits(tt, test_statistic))
    expect_identical(tt@name, nm)
    expect_true(nzchar(tt@label))
  }
  expect_identical(vapply(tests, function(tt) tt@label, character(1L),
                          USE.NAMES = FALSE),
                   c("Wald", "Likelihood-ratio", "Rao score",
                     "Terrell gradient"))
  expect_error(test_statistic(name = "x", label = "X"), "abstract")
})

test_that("test_restricted() is TRUE for exactly the three refitting tests", {
  expect_false(test_restricted(wald_test()))
  expect_true(test_restricted(lr_test()))
  expect_true(test_restricted(score_test()))
  expect_true(test_restricted(gradient_test()))
})

test_that("the three variance constructors build objects of the right classes", {
  types <- list(bayesian = bayesian_vcov(), frequentist = frequentist_vcov(),
                unconditional = unconditional_vcov())
  for (nm in names(types)) {
    expect_true(S7::S7_inherits(types[[nm]], vcov_type))
    expect_identical(types[[nm]]@name, nm)
  }
  expect_error(vcov_type(name = "x"), "abstract")
})

test_that("the objects print one line", {
  expect_identical(utils::capture.output(print(lr_test())), "<test> Likelihood-ratio")
  expect_identical(utils::capture.output(print(frequentist_vcov())),
                   "<vcov> frequentist")
})

test_that("a string for test is rejected with the constructor to use", {
  fit <- small_fit()
  expect_error(statmod_test(fit, param = "mu", coefname = "x", test = "lr"),
               "'test' must be a test_statistic object.*\"lr\" is no longer accepted.*lr_test\\(\\)")
  expect_error(summary(fit, test = "wald"), "wald_test\\(\\)")
  expect_error(confint(fit, test = "score", readable = FALSE), "score_test\\(\\)")
  expect_error(statmod_invert(fit, "mu", "x", test = "gradient"),
               "gradient_test\\(\\)")
  ## a value that is neither a string with a known name nor an object
  expect_error(statmod_test(fit, param = "mu", coefname = "x", test = 1),
               "wald_test\\(\\), lr_test\\(\\), score_test\\(\\) or gradient_test\\(\\)")
  expect_error(summary(fit, test = bayesian_vcov()), "must be a test_statistic")
})

test_that("a string for vcov_type is rejected with the constructor to use", {
  fit <- small_fit()
  expect_error(vcov(fit, vcov_type = "frequentist"),
               "'vcov_type' must be a vcov_type object.*\"frequentist\" is no longer accepted.*frequentist_vcov\\(\\)")
  expect_error(confint(fit, vcov_type = "bayesian"), "bayesian_vcov\\(\\)")
  expect_error(summary(fit, vcov_type = "unconditional"), "unconditional_vcov\\(\\)")
  expect_error(statmod_test(fit, param = "mu", coefname = "x",
                            vcov_type = "bayesian"), "bayesian_vcov\\(\\)")
  expect_error(predict(fit, what = "mu", se = TRUE, vcov_type = "frequentist"),
               "frequentist_vcov\\(\\)")
  expect_error(vcov(fit, vcov_type = lr_test()), "must be a vcov_type")
  expect_error(vcov(fit, vcov_type = 2), "bayesian_vcov\\(\\), frequentist_vcov\\(\\)")
})

test_that("the default objects give what the defaults gave", {
  fit <- small_fit()
  expect_identical(vcov(fit), vcov(fit, vcov_type = bayesian_vcov()))
  expect_identical(confint(fit), confint(fit, vcov_type = bayesian_vcov(),
                                         test = wald_test()))
  expect_identical(summary(fit)@tables,
                   summary(fit, vcov_type = bayesian_vcov(), test = wald_test())@tables)
  expect_identical(summary(fit)@vcov_type@name, "bayesian")
  expect_identical(summary(fit)@test@name, "wald")
  a <- statmod_test(fit, param = "mu", coefname = "x", value = 1)
  b <- statmod_test(fit, param = "mu", coefname = "x", value = 1,
                    test = wald_test(), vcov_type = bayesian_vcov())
  expect_identical(a@statistic, b@statistic)
  expect_identical(a@test@name, "wald")
  ## the Wald statistic is the squared estimate over the variance of the
  ## matrix that the convention selects
  for (ty in list(bayesian_vcov(), frequentist_vcov())) {
    v <- vcov(fit, vcov_type = ty)[["mu:x", "mu:x"]]
    s <- statmod_test(fit, param = "mu", coefname = "x", value = 1,
                      vcov_type = ty)@statistic
    expect_equal(s, (fit@coefficients$mu[[2L]] - 1)^2 / v)
  }
})

test_that("confint() reads whether the test is restricted from the object", {
  fit <- small_fit()
  expect_error(confint(fit, test = lr_test()), "test = lr_test\\(\\)")
  expect_error(confint(fit, test = score_test()), "readable = FALSE")
  expect_error(confint(fit, test = gradient_test()), "test = gradient_test\\(\\)")
  expect_no_error(confint(fit, test = wald_test()))
  expect_no_error(confint(fit, parm = "mu:x", test = lr_test(), readable = FALSE))
})

test_that("the result of a test keeps the object, and prints its label", {
  fit <- small_fit()
  tt <- statmod_test(fit, param = "mu", coefname = "x", value = 1,
                     test = score_test())
  expect_true(S7::S7_inherits(tt@test, test_statistic))
  expect_identical(tt@test@name, "score")
  out <- utils::capture.output(print(tt))
  expect_true(any(grepl("Rao score test on a statmod fit", out, fixed = TRUE)))
  s <- summary(fit, test = gradient_test(), vcov_type = frequentist_vcov())
  out <- utils::capture.output(print(s))
  expect_true(any(grepl("frequentist variance;  gradient test, signed root",
                        out, fixed = TRUE)))
  expect_error(StatmodTest(test = "lr"), "test")
})

test_that("the old argument name type is rejected by name in every method", {
  fit <- small_fit()
  expect_error(vcov(fit, type = bayesian_vcov()), "'type' is now 'vcov_type'")
  expect_error(confint(fit, type = bayesian_vcov()), "'type' is now 'vcov_type'")
  expect_error(summary(fit, type = bayesian_vcov()), "'type' is now 'vcov_type'")
  expect_error(predict(fit, what = "mu", se = TRUE, type = bayesian_vcov()),
               "'type' is now 'vcov_type'")
  expect_error(residuals(fit, type = "pearson"), "'type' is now 'residual_type'")
  expect_error(logLik(fit, type = "marginal"), "'type' is now 'loglik_type'")
  ## any other unknown argument of predict() is still named
  expect_error(predict(fit, what = "mu", bogus = 3), "takes no argument 'bogus'")
})

test_that("residual_type and loglik_type select as type did", {
  fit <- small_fit()
  expect_identical(residuals(fit, residual_type = "response"),
                   residuals(fit, "response"))
  expect_gt(abs(sum(residuals(fit, residual_type = "pearson"))), 0)
  expect_error(residuals(fit, residual_type = "deviance"), "should be one of")
  expect_s3_class(logLik(fit, loglik_type = "conditional"), "logLik")
  expect_true(is.finite(as.numeric(logLik(fit, loglik_type = "marginal"))))
  expect_identical(predict(fit, what = "mu", se = TRUE,
                           vcov_type = frequentist_vcov())$se,
                   predict(fit, what = "mu", se = TRUE,
                           vcov_type = frequentist_vcov())$se)
  expect_error(predict(fit, what = "mu", vcov_type = frequentist_vcov()),
               "read only with")
})
