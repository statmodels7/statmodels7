# jfp() beside or inside a filter (gas()): the derivatives of the penalty on
# the joint vector of the coefficients and the filter's free parameters, and
# the subformulas in which a block is rejected. Split from test-jfp-model.R.

test_that("a Jeffreys block inside a filter's subformula is exact", {
  skip_if_not_installed("numDeriv")
  set.seed(4)
  n <- 80
  d <- data.frame(t = 1:n, g = factor(rep(1:2, n / 2)), x = rnorm(n))
  e <- numeric(n)
  for (i in 2:n) e[i] <- 0.7 * e[i - 1] + rnorm(1, sd = 0.6)
  d$y <- 0.3 * d$x + e + rnorm(n, sd = 0.5)
  g1 <- distributions7::gaussian1_distrib()
  # the point: the fit of the same model without jfp(), whose coordinates
  # are in the same order
  f0 <- statmod(y ~ x + gas(p = 1, q = 1, kappa1 ~ 1 + g, time = t), g1, d)
  f <- statmod(y ~ x + gas(p = 1, q = 1, kappa1 ~ 1 + jfp(~ g), time = t),
               g1, d)
  sp <- f@spec
  des <- statmod_design(sp)
  b <- jfp_blocks(sp, des)
  expect_identical(b[[1L]]$kind, "filter")
  expect_true(b[[1L]]$const)
  sst <- statmod_structural_state(des)
  key <- names(sst$zeta)[1]
  free <- setdiff(names(sst$zeta[[key]]), sst$held[[key]])
  nb <- sum(vapply(des, function(z) z$npar, integer(1)))
  # the fit's own point, read before the evaluations below move the state
  u1 <- c(unlist(f@coefficients[sp@distrib@params], use.names = FALSE),
          as.numeric(sst$zeta[[key]][free]))
  s0 <- statmod_structural_state(statmod_design(f0@spec))
  u0 <- c(unlist(f0@coefficients[sp@distrib@params], use.names = FALSE),
          as.numeric(s0$zeta[[1L]][setdiff(names(s0$zeta[[1L]]),
                                           s0$held[[1L]])]))
  expect_length(u0, nb + length(free))
  obj <- statmod_objective(sp, f@hyper, des)
  jp <- statmod_joint_pieces(sp, des, obj, f@hyper)
  cf_at <- function(u) {
    jp$setz(u[nb + seq_along(free)])
    obj$split(u[seq_len(nb)])
  }
  val <- function(u) jfp_penalty(sp, des, cf_at(u), "value")
  gr <- function(u) jfp_penalty(sp, des, cf_at(u), "gradient", joint = TRUE)
  g <- gr(u0)
  gn <- numDeriv::grad(val, u0)
  expect_lt(max(abs(g - gn)) / max(abs(gn)), 1e-7)
  hn <- numDeriv::jacobian(gr, u0)
  hn <- (hn + t(hn)) / 2
  h <- jfp_penalty(sp, des, cf_at(u0), "hessian", joint = TRUE)
  expect_lt(max(abs(h - hn)) / max(abs(hn)), 1e-7)
  # the negative control: without the movement of the block's rows
  st <- jfp_state_at(sp, des, cf_at(u0))
  hp <- jfp_hessian(st)
  expect_gt(max(abs(-(hp$hessian - hp$refresh) - hn)) / max(abs(hn)), 1e-1)
  # the outer orders: D3 against a central difference of the exact Hessian,
  # D4 against one of the exact D3, along random directions
  set.seed(11)
  k <- length(u0)
  v <- rnorm(k)
  w <- rnorm(k)
  hh <- function(u) jfp_penalty(sp, des, cf_at(u), "hessian", joint = TRUE)
  d3 <- function(u) -jfp_deriv3_dir(jfp_state_at(sp, des, cf_at(u)), v)
  cdif <- function(f, dir, h = 1e-4) {
    a <- (f(u0 + h * dir) - f(u0 - h * dir)) / (2 * h)
    b <- (f(u0 + 2 * h * dir) - f(u0 - 2 * h * dir)) / (4 * h)
    (4 * a - b) / 3
  }
  r3 <- cdif(hh, v)
  r4 <- cdif(d3, w)
  D3 <- d3(u0)
  D4 <- -jfp_deriv4_dir(jfp_state_at(sp, des, cf_at(u0)), v, w)
  expect_lt(max(abs(D3 - r3)) / max(abs(r3)), 1e-8)
  expect_lt(max(abs(D4 - r4)) / max(abs(r4)), 1e-8)
  M <- crossprod(matrix(rnorm(k * k), k)) / k
  expect_equal(-jfp_deriv4(jfp_state_at(sp, des, cf_at(u0)), M, v, w),
               sum(M * D4), tolerance = 1e-10)
  # the negative control: without the filter's derivatives of the rows
  ns <- environment(jfp_mv_rows)
  real <- jfp_mv_rows
  unlockBinding("jfp_mv_rows", ns)
  assign("jfp_mv_rows", function(cx, i, S)
    lapply(real(cx, i, S), function(x) 0 * x), envir = ns)
  C3 <- d3(u0)
  assign("jfp_mv_rows", real, envir = ns)
  lockBinding("jfp_mv_rows", ns)
  expect_gt(max(abs(C3 - r3)) / max(abs(r3)), 1e-2)
  # the fit is stationary over the joint vector, the Jeffreys part included
  expect_lt(max(abs(jp$gr(u1))), 1e-5)
})

test_that("a block in a joint filter with a scaled score is exact at order 4", {
  # D4 of a block in the subformula of a joint filter reads the fifth order
  # of the filter, whose driving vector with scaling 1 needs the fifth
  # derivative of the expected information
  set.seed(5)
  n <- 100
  d <- data.frame(t = 1:n, z = rnorm(n))
  e <- numeric(n)
  for (i in 2:n) e[i] <- 0.7 * e[i - 1] + rnorm(1, sd = 0.5)
  d$y <- e + 0.3 * d$z + rnorm(n, sd = exp(0.2 * sin(d$t / 9)))
  g1 <- distributions7::gaussian1_distrib()
  # the derivatives are read away from the mode, so the hyperparameters
  # stay at their start and the joint Newton stops after 20 iterations
  f <- suppressWarnings(statmod(
    y ~ gas(p = 1, q = 1, kappa1 ~ 1 + jfp(~ z), time = t, scaling = 1) |
      sigma ~ gas(p = 1, q = 1, time = t, scaling = 1),
    distrib = g1, data = d, outer_criterion = NULL,
    inner_optimizer = optimizers7::newton(maxit = 20)))
  sp <- f@spec
  des <- statmod_design(sp)
  expect_true(jfp_filter_outer_ok(sp, des, 2L))
  sst <- statmod_structural_state(des)
  key <- names(sst$zeta)[1]
  free <- setdiff(names(sst$zeta[[key]]), sst$held[[key]])
  nb <- sum(vapply(des, function(z) z$npar, integer(1)))
  u0 <- c(unlist(f@coefficients[sp@distrib@params], use.names = FALSE),
          as.numeric(sst$zeta[[key]][free]))
  expect_length(u0, nb + length(free))
  obj <- statmod_objective(sp, f@hyper, des)
  jp <- statmod_joint_pieces(sp, des, obj, f@hyper)
  cf_at <- function(u) {
    jp$setz(u[nb + seq_along(free)])
    obj$split(u[seq_len(nb)])
  }
  set.seed(11)
  u0 <- u0 + 0.05 * rnorm(length(u0))
  k <- length(u0)
  v <- rnorm(k)
  w <- rnorm(k)
  d3 <- function(u) -jfp_deriv3_dir(jfp_state_at(sp, des, cf_at(u)), v)
  cdif <- function(f, dir, h = 1e-4) {
    a <- (f(u0 + h * dir) - f(u0 - h * dir)) / (2 * h)
    b <- (f(u0 + 2 * h * dir) - f(u0 - 2 * h * dir)) / (4 * h)
    (4 * a - b) / 3
  }
  r4 <- cdif(d3, w)
  D4 <- -jfp_deriv4_dir(jfp_state_at(sp, des, cf_at(u0)), v, w)
  expect_lt(max(abs(D4 - r4)) / max(abs(r4)), 1e-8)
  # the negative control: without the filter's derivatives of the rows
  ns <- environment(jfp_mv_rows)
  real <- jfp_mv_rows
  unlockBinding("jfp_mv_rows", ns)
  assign("jfp_mv_rows", function(cx, i, S)
    lapply(real(cx, i, S), function(x) 0 * x), envir = ns)
  C4 <- -jfp_deriv4_dir(jfp_state_at(sp, des, cf_at(u0)), v, w)
  assign("jfp_mv_rows", real, envir = ns)
  lockBinding("jfp_mv_rows", ns)
  expect_gt(max(abs(C4 - r4)) / max(abs(r4)), 1e-2)
})

test_that("jfp() inside a structural term other than gas() is rejected", {
  set.seed(1)
  d <- data.frame(t = 1:60, z = rnorm(60))
  d$y <- rnorm(60)
  expect_error(statmod(y ~ regime(k = 2, time = t, mu ~ 1 + jfp(~ z)),
                       distrib = distributions7::gaussian1_distrib(),
                       data = d))
})

test_that("jfp() in the subformula of a filter's level is rejected", {
  set.seed(1)
  d <- data.frame(t = 1:60, z = rnorm(60))
  d$y <- rpois(60, 2)
  pois <- distributions7::poisson_distrib()
  expect_error(statmod(y ~ gas(p = 1, q = 1, omega ~ 1 + jfp(~ z), time = t),
                       distrib = pois, data = d),
               "cannot be integrated at b = 1")
  expect_error(statmod(y ~ gas(p = 1, q = 1, omega ~ 0 + jfp(~ z), time = t),
                       distrib = pois, data = d),
               "subformula of omega")
  # in a joint filter, the level of one of its equations
  g1 <- distributions7::gaussian1_distrib()
  d$y2 <- rnorm(60)
  expect_error(statmod(y2 ~ gas(p = 1, q = 1, omega ~ 1 + jfp(~ z),
                                time = t) |
                         sigma ~ gas(p = 1, q = 1, time = t),
                       distrib = g1, data = d),
               "subformula of omega")
})

test_that("a Jeffreys block beside a filter in another equation is exact", {
  skip_if_not_installed("numDeriv")
  set.seed(8)
  n <- 200
  d <- data.frame(t = 1:n, x = rnorm(n), z = rnorm(n))
  e <- numeric(n)
  for (i in 2:n) e[i] <- 0.6 * e[i - 1] + rnorm(1, sd = 0.3)
  d$y <- 1 + 0.5 * d$x + rnorm(n, sd = exp(0.2 * d$z + e))
  f <- statmod(y ~ jfp(~ x) | sigma ~ z + gas(p = 1, q = 1, time = t),
               distrib = distributions7::gaussian1_distrib(), data = d)
  sp <- f@spec
  des <- statmod_design(sp)
  sst <- statmod_structural_state(des)
  key <- names(sst$zeta)[1]
  free <- setdiff(names(sst$zeta[[key]]), sst$held[[key]])
  nb <- sum(vapply(des, function(z) z$npar, integer(1)))
  u0 <- c(unlist(f@coefficients[sp@distrib@params], use.names = FALSE),
          as.numeric(sst$zeta[[key]][free]))
  obj <- statmod_objective(sp, f@hyper, des)
  jp <- statmod_joint_pieces(sp, des, obj, f@hyper)
  cf_at <- function(u) {
    jp$setz(u[nb + seq_along(free)])
    obj$split(u[seq_len(nb)])
  }
  val <- function(u) jfp_penalty(sp, des, cf_at(u), "value")
  gr <- function(u) jfp_penalty(sp, des, cf_at(u), "gradient", joint = TRUE)
  g <- gr(u0)
  gn <- numDeriv::grad(val, u0)
  expect_lt(max(abs(g - gn)) / max(abs(gn)), 1e-7)
  hn <- numDeriv::jacobian(gr, u0)
  hn <- (hn + t(hn)) / 2
  h <- jfp_penalty(sp, des, cf_at(u0), "hessian", joint = TRUE)
  expect_lt(max(abs(h - hn)) / max(abs(hn)), 1e-7)
  # the negative control: without the filter's second derivative
  st <- jfp_state_at(sp, des, cf_at(u0))
  hp <- jfp_hessian(st)
  expect_gt(max(abs(-(hp$hessian - hp$refresh) - hn)) / max(abs(hn)), 1e-2)
  # the fit is stationary over the joint vector, the Jeffreys part included
  expect_lt(max(abs(jp$gr(u0))), 1e-6)
  expect_gt(max(abs(g[nb + seq_along(free)])), 1e-3)
  expect_identical(statmod_certificate(f)$state, "converged")
})
