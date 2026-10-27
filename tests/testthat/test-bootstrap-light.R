g1l <- distributions7::gaussian1_distrib()

light_fit <- function(n = 300, seed = 3) {
  set.seed(seed)
  d <- data.frame(x = runif(n))
  d$y <- sin(2 * pi * d$x) + rnorm(n, sd = 0.3)
  # the formula carries the global environment, so that the environment of the
  # term is the shim and its parent only and not the frame of this function
  fml <- y ~ s(x, basis7::bspline_smooth(k = 10))
  environment(fml) <- globalenv()
  statmod(fml, g1l, d, outer_criterion = reml())
}

ser_mb <- function(x) length(serialize(x, NULL)) / 2^20

# A package installed with its source kept (R_KEEP_PKG_SOURCE=yes, as on CI)
# serializes the source of the whole package with any closure that carries a
# source reference. The size of one function of the package is that baseline.
src_mb <- function() ser_mb(statmodels7:::boot_light)

test_that("the evaluation environment of the formula serializes in a few kilobytes", {
  env <- statmodels7:::terms_first(globalenv())
  expect_lt(ser_mb(env), 0.5)
  # every function export of modelterms7 is bound, to the object it names
  ex <- getNamespaceExports("modelterms7")
  is_fn <- vapply(ex, function(nm)
    is.function(get(nm, envir = asNamespace("modelterms7"))), logical(1))
  expect_setequal(ls(env), ex[is_fn])
  for (nm in c("s", "lasso", "random", "term_build", "LinparTerm")) {
    expect_identical(get(nm, envir = env), get(nm, envir = asNamespace("modelterms7")))
    expect_true(exists(nm, envir = env, mode = "function"))
  }
  expect_true(bindingIsActive("s", env))
  # the terms of a fit keep the environment, so the specification stays light
  f <- light_fit()
  fenv <- environment(f@spec@terms$mu[["linpar"]]@formula)
  expect_true(bindingIsActive("s", fenv))
  expect_lt(ser_mb(fenv), 0.5)
})

test_that("a formula with several constructors is read through the bindings", {
  set.seed(2)
  d <- data.frame(g = factor(rep(1:8, each = 6)), x1 = rnorm(48), x2 = rnorm(48))
  d$y <- 1 + d$x1 + rnorm(8)[d$g] + rnorm(48, sd = 0.5)
  f <- suppressWarnings(statmod(y ~ lasso(~ x1 + x2) + random(~ 1 | g), g1l, d))
  tm <- f@spec@terms$mu
  expect_true(any(vapply(tm, S7::S7_inherits, logical(1),
                         modelterms7::PenalizedTerm)))
  expect_true(any(vapply(tm, S7::S7_inherits, logical(1),
                         modelterms7::RandomTerm)))
})

test_that("the fit the workers receive keeps what a refit reads", {
  f <- light_fit()
  lf <- statmodels7:::boot_light(f)
  expect_null(lf@methods$search)
  keep <- setdiff(names(f@methods), "search")
  expect_identical(lf@methods[keep], f@methods[keep])
  expect_identical(lf@coefficients, f@coefficients)
  expect_identical(lf@hyper, f@hyper)
  expect_lt(ser_mb(lf), ser_mb(f))
  expect_lt(ser_mb(lf), 4 + src_mb())
})

test_that("the replicates do not depend on the light fit", {
  f <- light_fit(n = 150)
  for (resample in list(case_resample(), parametric_resample(),
                        bayesian_resample())) {
    for (refit in c("coefficients", "full")) {
      set.seed(9)
      a <- statmodels7:::bootstrap_engine(f, 6, resample, refit)
      set.seed(9)
      b <- testthat::with_mocked_bindings(
        statmodels7:::bootstrap_engine(f, 6, resample, refit),
        boot_light = function(object) object, .package = "statmodels7")
      expect_identical(a$coefficients, b$coefficients)
      expect_identical(a$hyper, b$hyper)
      expect_identical(a$status, b$status)
    }
  }
})
