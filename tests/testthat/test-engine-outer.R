# The outer criterion on the compiled route against the R route: the value,
# the gradient and the Hessian at the fitted point, on each class of
# penalized unit that the compiled context covers.

outer_twin <- function(fit, method = fit@methods$outer) {
  sp <- fit@spec
  dz <- statmod_design(sp)
  attr(dz, "pinned") <- fit@methods$pinned
  gam <- marginal_coords(sp, dz, method)
  if (length(gam$where)) {
    sp <- statmod_hold(sp, gam, unlist(fit@coefficients[sp@distrib@params],
                                       use.names = FALSE)[gam$where])
  }
  idx <- outer_hyper_index(sp, statmod_blocks(sp, dz))
  basis <- integrated_basis(sp, dz, method@kind, gamma = length(gam$where) > 0)
  cf <- fit@coefficients
  hy <- fit@hyper
  sr <- sp
  sr@engine <- "r"
  ctx <- outer_context(sp, dz, cf, hy, NULL)
  covered <- !is.null(engine_outer_at(ctx, sp, dz, cf, hy, method, basis))
  list(
    covered = covered,
    v = c(statmod_marginal(sp, dz, cf, hy, method, basis = basis, ctx = ctx)$value,
          statmod_marginal(sr, dz, cf, hy, method, basis = basis)$value),
    g = list(statmod_marginal_grad(sp, dz, cf, hy, method, idx, basis, ctx = ctx,
                                   gam = gam),
             statmod_marginal_grad(sr, dz, cf, hy, method, idx, basis, gam = gam)),
    h = list(statmod_marginal_hess(sp, dz, cf, hy, method, idx, basis, ctx = ctx,
                                   gam = gam),
             statmod_marginal_hess(sr, dz, cf, hy, method, idx, basis, gam = gam)))
}

# the gate's measure: the largest difference against the size of the
# reference, with a floor of one (a gradient at the optimum is near zero)
twin_rel <- function(a, b) max(abs(a - b)) / max(1, abs(b))

expect_outer_twin <- function(fit, ...) {
  tw <- outer_twin(fit, ...)
  expect_true(tw$covered)
  expect_lt(twin_rel(tw$v[1], tw$v[2]), 1e-12)
  expect_lt(twin_rel(tw$g[[1]], tw$g[[2]]), 1e-9)
  expect_lt(twin_rel(tw$h[[1]], tw$h[[2]]), 1e-9)
}

outer_data <- function(n = 900, groups = 30, seed = 11) {
  set.seed(seed)
  g <- factor(rep_len(seq_len(groups), n))
  x <- runif(n)
  z <- runif(n, -1, 1)
  x2 <- runif(n, -1, 1)
  # three correlated random effects, intercept and two slopes
  L <- t(chol(0.64 * (diag(3) * 0.5 + 0.5)))
  b <- t(L %*% matrix(rnorm(3 * groups), 3))
  eta <- sin(5 * x) + b[g, 1] + b[g, 2] * z + b[g, 3] * x2
  data.frame(g = g, x = x, z = z, x2 = x2,
             y = eta + rnorm(n, sd = 0.3),
             yt = sin(4 * x) * cos(3 * z) + rnorm(n, sd = 0.3),
             yg = rgamma(n, shape = 3, rate = 3 / exp(0.3 * eta)))
}

test_that("a smooth in two equations: reml and ml", {
  d <- outer_data()
  f <- statmod(yg ~ s(x) | phi ~ s(z), distrib = distributions7::gamma1_distrib(),
               data = d)
  expect_outer_twin(f)
  f2 <- statmod(yg ~ s(x) | phi ~ s(z), distrib = distributions7::gamma1_distrib(),
                data = d, outer_criterion = ml())
  expect_outer_twin(f2)
})

test_that("a random intercept and a tensor product", {
  d <- outer_data()
  f <- statmod(y ~ s(x) + random(~ 1 | g), distrib = distributions7::gaussian1_distrib(),
               data = d)
  expect_outer_twin(f)
  f2 <- statmod(yt ~ te(x, z), distrib = distributions7::gaussian1_distrib(), data = d)
  expect_outer_twin(f2)
})

test_that("correlated and structured random effects", {
  d <- outer_data()
  f <- statmod(y ~ z + random(~ z | g), distrib = distributions7::gaussian1_distrib(),
               data = d)
  expect_outer_twin(f)
  pr <- distributions7::fixed(distributions7::mvgaussian2_distrib(
    n_dim = 2, omega = parameters7::log_cholesky(2)), mu = 0)
  cs <- distributions7::fixed(distributions7::mvgaussian1_distrib(
    n_dim = 3, sigma = parameters7::compound_symmetry(3)), mu = 0)
  ar <- distributions7::fixed(distributions7::mvgaussian1_distrib(
    n_dim = 3, sigma = parameters7::ar1(3)), mu = 0)
  expect_outer_twin(statmod(y ~ z + random(~ z | g, distrib = pr),
                            distrib = distributions7::gaussian1_distrib(), data = d))
  expect_outer_twin(statmod(y ~ z + random(~ z + x2 | g, distrib = cs),
                            distrib = distributions7::gaussian1_distrib(), data = d))
  expect_outer_twin(statmod(y ~ z + random(~ z + x2 | g, distrib = ar),
                            distrib = distributions7::gaussian1_distrib(), data = d))
})

test_that("engine = 'r' leaves the outer criterion to the R route", {
  d <- outer_data()
  f <- statmod(y ~ s(x), distrib = distributions7::gaussian1_distrib(), data = d,
               engine = "r")
  sp <- f@spec
  dz <- statmod_design(sp)
  expect_null(engine_outer_at(NULL, sp, dz, f@coefficients, f@hyper,
                              f@methods$outer))
})
