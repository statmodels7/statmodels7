# A method whose `...` nothing reads rejects an argument by name, so that a
# misspelled or nonexistent argument is not accepted and ignored.

test_that("the methods of a fit reject an argument they do not take", {
  set.seed(2)
  d <- data.frame(x = stats::runif(40))
  d$y <- 1 + d$x + stats::rnorm(40, sd = 0.4)
  fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), d)
  bad <- list(
    coef = function() coef(fit, bogus = 1),
    vcov = function() vcov(fit, bogus = 1),
    logLik = function() logLik(fit, bogus = 1),
    residuals = function() residuals(fit, bogus = 1),
    predict = function() predict(fit, what = "mu", bogus = 1),
    nobs = function() nobs(fit, bogus = 1),
    formula = function() formula(fit, bogus = 1),
    family = function() family(fit, bogus = 1),
    weights = function() weights(fit, bogus = 1),
    df.residual = function() df.residual(fit, bogus = 1),
    sigma = function() sigma(fit, bogus = 1),
    model.matrix = function() model.matrix(fit, bogus = 1),
    simulate = function() simulate(fit, bogus = 1),
    loglik = function() loglik(fit, bogus = 1),
    gradient = function() gradient(fit, bogus = 1),
    hessian = function() hessian(fit, bogus = 1)
  )
  for (nm in names(bad)) {
    expect_error(bad[[nm]](), sprintf("%s\\(\\) takes no argument 'bogus'", nm),
                 info = nm)
  }
  # an unnamed argument is reported as such
  expect_error(coef(fit, TRUE, 3), "coef\\(\\) takes no argument '<unnamed>'")
  # confint() and summary() pass their dots to vcov(), which names the argument
  expect_error(confint(fit, bogus = 1), "takes no argument 'bogus'")
  expect_error(summary(fit, bogus = 1), "takes no argument 'bogus'")
  # and the arguments that are taken still work
  expect_true(is.matrix(vcov(fit, readable = FALSE, param = "mu")))
  expect_s3_class(logLik(fit, loglik_type = "conditional"), "logLik")
})
