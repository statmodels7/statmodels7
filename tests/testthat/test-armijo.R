## The sufficient-decrease condition of the scoring loop.
##
## The increment solves (H + S) delta = -g, so g'delta is NEGATIVE and
## Armijo's bound f(x) + c1 s g'delta sits BELOW f(x): the test asks for a
## decrease. Written with the term SUBTRACTED it changes sign and permits an
## increase of c1 s |g'delta| instead, which is what the loop did while its
## own comment said otherwise.
##
## None of the fits in this suite relied on that licence -- measured, 0 of 18
## accepted steps of the model that motivated it increased the objective --
## so a test over a fit could not have caught the sign. The condition is
## named so that it can be asked directly, and every case below carries the
## form it replaced as its negative control.

test_that("armijo_ok asks for a decrease and never permits an increase", {
  gd <- -1          # a descent direction: g'delta < 0
  old <- function(vnew, value, step) vnew <= value - 1e-4 * step * gd

  ## an increase, however small, is refused
  expect_false(armijo_ok(10 + 1e-6, 10, 1, gd))
  expect_true(old(10 + 1e-6, 10, 1))          # the form it replaced took it

  ## and so is standing still
  expect_false(armijo_ok(10, 10, 1, gd))
  expect_true(old(10, 10, 1))

  ## a decrease of the share asked for is accepted, one just short is not
  expect_true(armijo_ok(10 - 1e-4, 10, 1, gd))
  expect_false(armijo_ok(10 - 0.99e-4, 10, 1, gd))

  ## an ordinary Newton step, whose decrease is about |g'delta|/2, passes
  expect_true(armijo_ok(10 - 0.5, 10, 1, gd))

  ## the bound scales with the step, so a short step asks for little
  expect_true(armijo_ok(10 - 1e-8, 10, 1e-4, gd))
  expect_false(armijo_ok(10 - 1e-9, 10, 1e-4, gd))

  ## a non-finite trial point is never acceptable
  expect_false(armijo_ok(NaN, 10, 1, gd))
  expect_false(armijo_ok(Inf, 10, 1, gd))
  expect_false(armijo_ok(-Inf, 10, 1, gd))    # -Inf is not finite either
})

test_that("g'delta really is negative in a fit", {
  ## the premise the sign rests on, measured rather than assumed: delta
  ## solves (H + S) delta = -g with H + S positive definite, so
  ## g'delta = -g'(H + S)^-1 g < 0.
  seen <- new.env()
  seen$gd <- numeric(0)
  orig <- iwls_solve
  local_mocked_bindings(
    iwls_solve = function(pieces, u, how, damp = 0) {
      out <- orig(pieces, u, how, damp)
      seen$gd <- c(seen$gd, -sum(u * out$delta))   # u = -g
      out
    })
  set.seed(5)
  n <- 200
  dd <- data.frame(x = stats::runif(n, -1, 1))
  dd$y <- stats::rpois(n, exp(0.4 + 0.9 * dd$x))
  statmod(y ~ x, distributions7::poisson_distrib(), dd)
  expect_gt(length(seen$gd), 0L)
  expect_true(all(seen$gd < 0))
})

test_that("no accepted step increases the objective", {
  ## the property the condition delivers, read off the loop's own history
  set.seed(6)
  n <- 300
  dd <- data.frame(x = stats::runif(n, -2, 2))
  dd$y <- stats::rgamma(n, shape = 3, rate = 3 / exp(0.5 + 0.6 * dd$x))
  fit <- statmod(y ~ x | phi ~ 1, distributions7::gamma1_distrib(), dd)
  hh <- fit@history$inner
  expect_s3_class(hh, "data.frame")
  expect_true(nrow(hh) > 1L)
  ## within each block of each pass, the objective never rises
  for (k in split(hh, list(hh$pass, hh$block), drop = TRUE)) {
    if (nrow(k) > 1L) expect_true(all(diff(k$objective) <= 0))
  }
})
