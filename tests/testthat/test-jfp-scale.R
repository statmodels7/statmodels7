# The factor that only scales a parameter's information is divided out of the
# Jeffreys penalty (distributions7::information_scale())

test_that("a gaussian with a jfp block keeps the unpenalized REML sigma", {
  set.seed(1)
  d <- data.frame(x = stats::rnorm(30), z = stats::rnorm(30))
  d$y <- 1 + d$x + stats::rnorm(30)
  G <- distributions7::gaussian1_distrib()
  f1 <- statmod(y ~ x + jfp(~ z), distrib = G, data = d)
  expect_equal(unname(exp(coef(f1)$sigma)), summary(stats::lm(y ~ x + z, d))$sigma,
               tolerance = 1e-6)
})

jfp_scale_case <- function(fm, distrib, d) {
  fit <- suppressWarnings(statmod(fm, distrib = distrib, data = d,
                                  outer_criterion = NULL))
  sp <- fit@spec
  des <- statmod_design(sp)
  params <- sp@distrib@params
  cf <- fit@coefficients[params]
  npar <- lengths(cf)
  m <- sum(npar)
  offs <- cumsum(npar) - npar
  split <- function(u) stats::setNames(lapply(seq_along(params), function(a)
    u[offs[a] + seq_len(npar[a])]), params)
  set.seed(1)
  u0 <- unlist(cf, use.names = FALSE) + stats::rnorm(m, sd = 0.05)
  st_at <- function(u) jfp_state(sp, des, split(u), jfp_blocks(sp, des),
                                 emethod = jfp_emethod(sp@distrib))
  st <- st_at(u0)
  st0 <- st_at(u0)
  st0$scale <- list()
  st0$E <- st0$raw0
  st0$blocks <- lapply(st0$blocks, function(b) jfp_block_state(st0, b$block))
  list(st = st, st0 = st0, u0 = u0, st_at = st_at, m = m,
       disp = offs[2L] + seq_len(npar[2L]))
}

test_that("the scaled penalty's derivatives are exact and free of the dispersion", {
  skip_if_not_installed("numDeriv")
  set.seed(2)
  n <- 80
  d <- data.frame(x = stats::runif(n, 0.2, 3), z = stats::rbinom(n, 1, 0.3),
                  w = stats::rnorm(n))
  d$y <- stats::rgamma(n, shape = 3, rate = 3 / exp(0.5 + 0.3 * d$x))
  r <- jfp_scale_case(y ~ x + jfp(~ z) | phi ~ w,
                      distributions7::invgauss1_distrib(), d)
  f <- function(u) jfp_value(r$st_at(u))
  gn <- numDeriv::grad(f, r$u0)
  g <- jfp_gradient(r$st)
  expect_lt(max(abs(g - gn)) / max(abs(gn)), 1e-8)
  hn <- numDeriv::hessian(f, r$u0)
  expect_lt(max(abs(jfp_hessian(r$st)$hessian - hn)) / max(abs(hn)), 1e-6)
  v <- stats::rnorm(r$m)
  n3 <- matrix(numDeriv::jacobian(function(t)
    as.vector(jfp_hessian(r$st_at(r$u0 + t * v))$hessian), 0), r$m)
  expect_lt(max(abs(jfp_deriv3_dir(r$st, v) - n3)) / max(abs(n3)), 1e-7)
  expect_lt(max(abs(g[r$disp])), 1e-12)
  # the negative control: the full information moves with the dispersion
  g0 <- jfp_gradient(r$st0)
  expect_gt(max(abs(g0[r$disp])), 0.1)
  expect_gt(max(abs(g0 - gn)) / max(abs(gn)), 0.1)
})

test_that("a family that declares no factor keeps the full information", {
  set.seed(3)
  n <- 80
  d <- data.frame(x = stats::runif(n, 0.2, 3), z = stats::rbinom(n, 1, 0.3))
  d$y <- stats::rnbinom(n, size = 3, mu = exp(1 - 0.5 * d$x))
  r <- jfp_scale_case(y ~ x + jfp(~ z), distributions7::negbin2_distrib(), d)
  expect_length(r$st$scale, 0L)
  expect_equal(jfp_gradient(r$st), jfp_gradient(r$st0), tolerance = 1e-14)
})
