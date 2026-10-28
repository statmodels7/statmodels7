# The weight of a stored fit. outer_fit() stores the optimizer of the outer
# search in fit@methods$search, and a function in it that holds the frame of
# the search takes the design, the inner fit and the derivative caches into
# every serialization of the fit (31 MB for a smooth of 2000 observations).

# Size of the serialization with every source file replaced by a name, so that
# the figure does not depend on whether the package was installed with its
# source kept (R_KEEP_PKG_SOURCE=yes, as on CI), which adds the source of the
# package to any closure that carries a source reference.
ser_mb <- function(x) {
  hook <- function(e) if (inherits(e, "srcfile")) "srcfile" else NULL
  length(serialize(x, NULL, refhook = hook)) / 2^20
}

weight_data <- function(n = 1000, seed = 3) {
  set.seed(seed)
  d <- data.frame(x = runif(n), g = factor(rep(1:20, length.out = n)))
  d$y <- sin(2 * pi * d$x) + rnorm(20)[as.integer(d$g)] + rnorm(n, sd = 0.3)
  d
}

# the formula carries the global environment, so that the environment of the
# term is the shim and not the frame of the caller
weight_fit <- function(d, rhs, ...) {
  fml <- stats::as.formula(paste("y ~", rhs), env = globalenv())
  statmod(fml, distributions7::gaussian1_distrib(), d, outer_criterion = reml(),
          ...)
}

test_that("the optimizer a fit stores holds no frame of the search", {
  d <- weight_data()
  for (rhs in c("s(x, basis7::bspline_smooth(k = 20))",
                "s(x, basis7::bspline_smooth(k = 10)) + random(~ 1 | g)")) {
    f <- weight_fit(d, rhs)
    ls <- f@methods$search@line_search
    expect_false(is.function(ls@resolution))
    expect_true(is.numeric(ls@resolution) && length(ls@resolution) == 1L &&
                  is.finite(ls@resolution) && ls@resolution >= 0)
    expect_lt(ser_mb(f@methods$search), 0.5)
    expect_lt(ser_mb(f), 3)
  }
})

test_that("a fit read back from a file predicts as the fit did", {
  d <- weight_data()
  f <- weight_fit(d, "s(x, basis7::bspline_smooth(k = 10))")
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))
  saveRDS(f, path)
  # the size of the file, compressed as saveRDS() does, without the sources
  gz <- memCompress(serialize(f, NULL, refhook = function(e)
    if (inherits(e, "srcfile")) "srcfile" else NULL), type = "gzip")
  expect_lt(length(gz) / 2^20, 1)
  f2 <- readRDS(path)
  expect_identical(coef(f2), coef(f))
  new <- d[1:5, ]
  expect_equal(predict(f2, what = "mu", newdata = new),
               predict(f, what = "mu", newdata = new))
})

test_that("freeze_resolution stores the number of a function the package set", {
  freeze <- statmodels7:::freeze_resolution
  nw <- optimizers7::newton()
  nw <- S7::set_props(nw, line_search = S7::set_props(
    nw@line_search, resolution = function() 0.5))
  expect_identical(freeze(nw, TRUE, 0.04)@line_search@resolution, 0.04)
  # a value that is not finite is no resolution, which a line search holds as 0
  expect_identical(freeze(nw, TRUE, NA_real_)@line_search@resolution, 0)
  expect_identical(freeze(nw, TRUE, Inf)@line_search@resolution, 0)
  # the rest of the optimizer is the same object
  expect_identical(freeze(nw, TRUE, 0.04)@criterion, nw@criterion)
  expect_identical(freeze(nw, TRUE, 0.04)@line_search@max_step,
                   nw@line_search@max_step)
  # an optimizer the caller supplied, one with a number, and one with no line
  # search are returned as they are
  expect_identical(freeze(nw, FALSE, 0.04), nw)
  num <- S7::set_props(nw, line_search = S7::set_props(
    nw@line_search, resolution = 0.3))
  expect_identical(freeze(num, TRUE, 0.04), num)
  nm <- optimizers7::nelder_mead()
  expect_identical(freeze(nm, TRUE, 0.04), nm)
})

test_that("an optimizer the caller names keeps the function it carries", {
  d <- weight_data(n = 300)
  fn <- function() 1e-10
  o <- optimizers7::newton()
  o <- S7::set_props(o, line_search = S7::set_props(
    o@line_search, resolution = fn))
  f <- weight_fit(d, "s(x, basis7::bspline_smooth(k = 8))",
                  outer_optimizer = o)
  expect_identical(f@methods$search@line_search@resolution, fn)
})
