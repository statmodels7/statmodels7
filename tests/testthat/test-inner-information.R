# The matrix of the inner fit is the argument inner_information of statmod().
# iwls() is a Newton method on whatever curvature the objective offers, and an
# optimizers7 optimizer receives the Hessian that the argument selects.

gam_data <- function(n = 120L, seed = 21L) {
  set.seed(seed)
  x <- stats::runif(n)
  data.frame(x = x, y = distributions7::distrib_rng(
    distributions7::gamma1_distrib(), n,
    list(mu = exp(0.5 + x), phi = exp(-0.5 + 0.6 * x))))
}

test_that("inner_information defaults to auto_info() and the fit records both objects", {
  dd <- gam_data()
  f <- statmod(y ~ x | phi ~ x, distributions7::gamma1_distrib(), dd)
  expect_identical(info_label(f@methods$inner_information), "auto")
  # gamma1 writes its expected information out at an acceptable cost
  expect_identical(info_label(f@methods$smooth@hessian), "expected")
  expect_true(f@methods$smooth@switching)
  o <- statmod(y ~ x | phi ~ x, distributions7::gamma1_distrib(), dd,
               inner_information = observed_info())
  expect_identical(info_label(o@methods$inner_information), "observed")
  expect_identical(info_label(o@methods$smooth@hessian), "observed")
  expect_null(o@methods$smooth@backup)
  expect_false(o@methods$smooth@switching)
  expect_equal(unlist(coef(f)), unlist(coef(o)), tolerance = 1e-5)
})

test_that("a named observed information never falls back to the expected one", {
  set.seed(3)
  n <- 100
  d <- data.frame(x = stats::runif(n))
  d$y <- 1 + d$x + 0.4 * stats::rt(n, 3)
  f <- statmod(y ~ x | sigma ~ x, distributions7::student_t1_distrib(), d,
               outer_criterion = NULL, inner_information = observed_info())
  expect_null(f@methods$smooth@backup)
  expect_identical(info_label(f@methods$smooth@hessian), "observed")
  # the choice of auto_info() on the same family does arm the fallback
  a <- statmod(y ~ x | sigma ~ x, distributions7::student_t1_distrib(), d,
               outer_criterion = NULL)
  expect_false(is.null(a@methods$smooth@backup) &&
                 info_label(a@methods$smooth@hessian) == "observed")
})

test_that("an optimizers7 optimizer receives the information that inner_information names", {
  dd <- gam_data()
  g <- distributions7::gamma1_distrib()
  fo <- statmod(y ~ x | phi ~ x, g, dd, outer_criterion = NULL,
                inner_optimizer = optimizers7::newton(),
                inner_information = observed_info())
  fe <- statmod(y ~ x | phi ~ x, g, dd, outer_criterion = NULL,
                inner_optimizer = optimizers7::newton(),
                inner_information = expected_info())
  fa <- statmod(y ~ x | phi ~ x, g, dd, outer_criterion = NULL,
                inner_optimizer = optimizers7::newton())
  expect_false(inner_settings(fo@methods$smooth)$expected)
  expect_true(inner_settings(fe@methods$smooth)$expected)
  # auto_info() is settled against the family as it is for iwls(): gamma1
  # writes its expected information out, so that is the matrix
  expect_true(inner_settings(fa@methods$smooth)$expected)
  expect_identical(info_label(attr(fa@methods$smooth, "inner_information")),
                   "expected")
  # the matrix changes the path and not the mode
  expect_equal(unlist(coef(fo)), unlist(coef(fe)), tolerance = 1e-4)
  expect_equal(unlist(coef(fa)), unlist(coef(fe)), tolerance = 1e-10)
  # a family without a closed form settles auto_info() on the observed matrix
  fp <- statmod(y ~ x, pig_bare_distrib(),
                data.frame(x = dd$x, y = stats::rpois(nrow(dd), 3)),
                outer_criterion = NULL,
                inner_optimizer = optimizers7::newton())
  expect_false(inner_settings(fp@methods$smooth)$expected)
})

test_that("the first step of an optimizer depends on the matrix it receives", {
  dd <- gam_data(200L, 5L)
  g <- distributions7::gamma1_distrib()
  one <- function(info) {
    f <- suppressWarnings(statmod(y ~ x | phi ~ x, g, dd, outer_criterion = NULL,
                                  inner_optimizer = optimizers7::newton(maxit = 1),
                                  inner_information = info))
    unlist(coef(f))
  }
  expect_false(isTRUE(all.equal(one(observed_info()), one(expected_info()),
                                tolerance = 1e-8)))
  expect_identical(one(expected_info()), one(expected_info()))
})

test_that("the refits of a fit read the information it recorded", {
  dd <- gam_data()
  g <- distributions7::gamma1_distrib()
  f <- statmod(y ~ x | phi ~ x, g, dd, outer_criterion = NULL,
               inner_optimizer = optimizers7::newton(),
               inner_information = expected_info())
  s <- statmod_test(f, param = "mu", coefname = "x", test = score_test())
  expect_true(is.finite(s@statistic))
  l <- statmod_test(f, param = "mu", coefname = "x", test = lr_test())
  expect_true(is.finite(l@statistic))
  p <- predict(f, what = "response", newdata = data.frame(x = c(0.2, 0.7)),
               interval = prediction_interval(predictive = pred_bootstrap(
                 n = 5, refit = "coefficients")))
  expect_equal(nrow(p), 2L)
  # a clone of the call reproduces the fit
  f2 <- eval(f@call)
  expect_identical(info_label(f2@methods$inner_information), "expected")
})

test_that("a structural term is fitted on the observed information whatever inner_information says", {
  skip_on_cran()
  set.seed(12)
  n <- 120
  y <- as.numeric(stats::arima.sim(list(ar = 0.5), n)) + 2
  d <- data.frame(y = y)
  g <- distributions7::gaussian1_distrib()
  a <- statmod(y ~ gas(1, 1), g, d, outer_criterion = NULL)
  b <- statmod(y ~ gas(1, 1), g, d, outer_criterion = NULL,
               inner_information = expected_info())
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)), tolerance = 1e-8)
})

test_that("inner_information must be an information object", {
  dd <- gam_data()
  g <- distributions7::gamma1_distrib()
  expect_error(statmod(y ~ x, g, dd, inner_information = "expected"),
               "'inner_information' must be an information object")
  expect_error(statmod(y ~ x, g, dd, inner_information = iwls()),
               "'inner_information' must be an information object")
})
