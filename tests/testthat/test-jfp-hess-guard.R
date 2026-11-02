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
