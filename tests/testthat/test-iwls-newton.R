# The Newton step on the full Hessian, which iwls_fit() tries where the
# scoring step is rejected or shrunk below a tenth: the scoring pieces of a
# block moving with its coefficients leave out the score times the block's
# own derivative.

test_that("a scoring curvature far too small is replaced by the full one", {
  Q <- matrix(c(4, 1, 1, 3), 2)
  m <- c(1, -2)
  toy <- list(fn = function(b) 0.5 * sum((b - m) * as.numeric(Q %*% (b - m))),
              gr = function(b) as.numeric(Q %*% (b - m)))
  # the scoring pieces carry a thousandth of the curvature, so their step is
  # a thousand times too long and the line search shrinks it to about that
  small <- function(b) list(R = NULL, C = NULL, A = Q / 1000)
  meth <- iwls(hessian = "observed", maxit = 20)
  r0 <- iwls_fit(toy, c(0, 0), meth, 1, small)
  r1 <- iwls_fit(toy, c(0, 0), meth, 1, small, newton_at = function(b) Q)
  expect_true(r1$converged)
  expect_equal(r1$par, m, tolerance = 1e-8)
  expect_gt(r1$fallback[["newton"]], 0L)
  # without it the run is still far from the minimum after the same budget
  expect_gt(max(abs(r0$par - m)), 1e-3)
  expect_identical(r0$fallback[["newton"]], 0L)
  # a run whose scoring steps are accepted never takes it
  right <- function(b) list(R = NULL, C = NULL, A = Q)
  r2 <- iwls_fit(toy, c(0, 0), meth, 1, right, newton_at = function(b) Q)
  expect_identical(r2$fallback[["newton"]], 0L)
})

test_that("an indefinite Hessian is repaired into a descent direction", {
  toy <- list(fn = function(b) sum(b^4) - sum(b^2),
              gr = function(b) 4 * b^3 - 2 * b)
  st <- iwls_newton_step(toy, c(0.3, -0.2), toy$fn(c(0.3, -0.2)),
                         toy$gr(c(0.3, -0.2)), diag(12 * c(0.3, -0.2)^2 - 2),
                         integer(0), 30L)
  expect_true(st$ok)
  expect_lt(st$vnew, toy$fn(c(0.3, -0.2)))
  # a held coordinate does not move
  st2 <- iwls_newton_step(toy, c(0.3, -0.2), toy$fn(c(0.3, -0.2)),
                          toy$gr(c(0.3, -0.2)), diag(12 * c(0.3, -0.2)^2 - 2),
                          2L, 30L)
  expect_identical(st2$cand[2], -0.2)
})

test_that("the full Hessian of a smoothed break-point is the objective's", {
  # Measured on a jump() under the quintic: in the break-point the scoring
  # curvature is 6.1 and the true one 7853, and the sum the Newton step
  # reads agrees with a difference of the analytic gradient to 1e-7.
  set.seed(3)
  d <- data.frame(x = stats::runif(200))
  d$y <- 1 + (d$x > 0.5) + stats::rnorm(200, sd = 0.4)
  fit <- suppressWarnings(statmod(
    y ~ jump(x, smoothed = numericals7::smooth_quintic()),
    distributions7::gaussian1_distrib(), d, outer_criterion = NULL))
  spec <- fit@spec
  design <- statmod_design(spec)
  obj <- statmod_objective(spec, fit@hyper, design, FALSE, "opg")
  par <- obj$stack(fit@coefficients)
  par[3] <- par[3] + 0.002
  npar <- obj$npar
  D <- mode_curvature(spec, design, obj$split(par), spec@distrib@params, npar,
                      cumsum(npar) - npar, sum(npar))
  H <- as.matrix(obj$he(par)) + D
  Hd <- sapply(seq_along(par), function(j) {
    e <- replace(numeric(length(par)), j, 1e-6)
    (obj$gr(par + e) - obj$gr(par - e)) / 2e-6
  })
  Hd <- (Hd + t(Hd)) / 2
  expect_lt(max(abs(H - Hd)) / max(abs(Hd)), 1e-6)
  # and the scoring curvature alone is not it
  expect_gt(max(abs(as.matrix(obj$he(par)) - Hd)) / max(abs(Hd)), 1e-2)
})
