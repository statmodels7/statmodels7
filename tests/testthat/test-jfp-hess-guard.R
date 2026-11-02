# The observed outer Hessian reads a Jeffreys block's fourth derivative only
# where jfp_outer_ok() admits order 2; a term whose block moves without its
# own term_block_deriv3() falls back to differencing the exact gradient

test_that("the observed outer Hessian beside jfp() checks order 2", {
  set.seed(6)
  n <- 200
  d <- data.frame(x = stats::runif(n, -2, 2), h = stats::rbinom(n, 1, 0.3),
                  g = factor(rep(1:20, 10)))
  d$y <- stats::rbinom(n, 1, stats::plogis(2 * (d$x - 0.3) +
                                             stats::rnorm(20, sd = 0.5)[d$g]))
  B <- distributions7::bernoulli_distrib()
  fm <- y ~ 0 + nl(~ b * (x - e), start = list(b = 1, e = 0)) + jfp(~ h) +
    random(~1 | g)
  method <- ml(information = observed_info())
  fit0 <- statmod(fm, distrib = B, data = d, outer_criterion = NULL)
  spec <- fit0@spec
  design <- statmod_design(spec)
  idx <- outer_hyper_index(spec, statmod_blocks(spec, design))
  basis <- integrated_basis(spec, design, method@kind)
  cf <- fit0@coefficients
  hy <- fit0@hyper
  calls <- 0L
  stencil <- statmod_hess_stencil
  spy <- function(...) {
    calls <<- calls + 1L
    stencil(...)
  }
  # nl() has its own third derivative of the block: the assembly answers
  expect_true(jfp_outer_ok(spec, design, 2L))
  H <- testthat::with_mocked_bindings(
    statmod_marginal_hess(spec, design, cf, hy, method, idx, basis),
    statmod_hess_stencil = spy)
  expect_identical(calls, 0L)
  # a term without it (here jfp_refresh_ok() answers FALSE at order 2, as
  # it does for a class whose term_block_deriv3() is the base method)
  Hs <- testthat::with_mocked_bindings(
    statmod_marginal_hess(spec, design, cf, hy, method, idx, basis),
    statmod_hess_stencil = spy,
    jfp_refresh_ok = function(spec, design, order = 1L) order < 2L)
  expect_identical(calls, 1L)
  expect_equal(Hs, stencil(spec, design, cf, hy, method, idx, basis),
               tolerance = 1e-12)
  # the two routes read the same second derivative
  expect_lt(max(abs(H - Hs)) / max(abs(H)), 1e-5)
})

test_that("beside scs() without its own third derivative the Hessian falls back", {
  set.seed(6)
  n <- 200
  d <- data.frame(z = stats::rbinom(n, 1, 0.3), x = stats::runif(n))
  d$y <- stats::rpois(n, exp(0.3 + 1.2 * d$x^2 + 0.5 * d$z))
  P <- distributions7::poisson_distrib()
  fm <- y ~ modelterms7::scs(x, modelterms7::shape_increasing(),
                             basis7::bspline_smooth(k = 5)) +
    modelterms7::jfp(~ z)
  method <- ml(information = observed_info())
  fit0 <- statmod(fm, distrib = P, data = d, outer_criterion = NULL)
  spec <- fit0@spec
  design <- statmod_design(spec)
  idx <- outer_hyper_index(spec, statmod_blocks(spec, design))
  basis <- integrated_basis(spec, design, method@kind)
  cf <- fit0@coefficients
  hy <- fit0@hyper
  calls <- 0L
  stencil <- statmod_hess_stencil
  spy <- function(...) {
    calls <<- calls + 1L
    stencil(...)
  }
  hess <- function() testthat::with_mocked_bindings(
    statmod_marginal_hess(spec, design, cf, hy, method, idx, basis),
    statmod_hess_stencil = spy)
  # scs() has its own term_block_deriv3(): the assembly answers
  expect_true(jfp_outer_ok(spec, design, 2L))
  H <- hess()
  expect_identical(calls, 0L)
  # its method removed, so that the class inherits the base one, which
  # returns zeros
  rf <- attr(design, "refresh")[[1L]]
  cls <- S7::S7_class(attr(design, "state")$terms[[rf$param]][[rf$term]])
  gen <- modelterms7::term_block_deriv3
  own <- S7::method(gen, cls)
  tab <- gen@methods
  key <- Filter(function(k) identical(get(k, envir = tab), own), ls(tab))
  expect_length(key, 1L)
  rm(list = key, envir = tab)
  on.exit(assign(key, own, envir = tab), add = TRUE)
  expect_false(jfp_outer_ok(spec, design, 2L))
  Hs <- hess()
  expect_identical(calls, 1L)
  expect_identical(Hs, stencil(spec, design, cf, hy, method, idx, basis))
  # what the guard keeps out: the assembly on the base method's zeros
  Hz <- testthat::with_mocked_bindings(
    statmod_marginal_hess(spec, design, cf, hy, method, idx, basis),
    jfp_outer_ok = function(spec, design, order = 1L) TRUE)
  expect_gt(max(abs(Hz - H)) / max(abs(H)), 1e-4)
  assign(key, own, envir = tab)
  expect_true(jfp_outer_ok(spec, design, 2L))
})
