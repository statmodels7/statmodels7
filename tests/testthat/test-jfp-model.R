# jfp() in statmod(): the penalized fit, its consumers, and the solve with a
# matrix plus a product of low rank.

jfp_sep_data <- function() {
  set.seed(11)
  n <- 60
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  # x1 separates the outcome quasi-completely
  d$y <- as.numeric(d$x1 > 0)
  d$y[1:6] <- 1 - d$y[1:6]
  d$y[d$x1 > 1] <- 1
  d
}

test_that("statmod() with jfp() over every coefficient is Firth's estimate", {
  skip_if_not_installed("numDeriv")
  d <- jfp_sep_data()
  f <- statmod(y ~ jfp(~ x1 + x2), distrib = distributions7::bernoulli_distrib(),
               data = d)
  b <- unlist(coef(f), use.names = FALSE)
  X <- cbind(1, d$x1, d$x2)
  closed <- function(b) {
    mu <- plogis(drop(X %*% b))
    sum(dbinom(d$y, 1, mu, log = TRUE)) +
      0.5 * as.numeric(determinant(crossprod(X * (mu * (1 - mu)), X))$modulus)
  }
  expect_lt(max(abs(numDeriv::grad(closed, b))), 1e-6)
  expect_true(all(is.finite(b)))
  expect_identical(statmod_certificate(f)$state, "converged")
  # the same estimate with the block written as a matrix
  d$X <- cbind(x1 = d$x1, x2 = d$x2)
  f2 <- statmod(y ~ jfp(X), distrib = distributions7::bernoulli_distrib(),
                data = d)
  expect_equal(unlist(coef(f2), use.names = FALSE), b, tolerance = 1e-8)
})

test_that("a factor in jfp() is coded with its contrasts", {
  set.seed(3)
  d <- data.frame(g = factor(rep(c("A", "B"), each = 80)), x = rnorm(160))
  lam <- exp(1.2 + 0.3 * d$x)
  d$y <- ifelse(runif(160) < 0.3, 0L, rpois(160, lam))
  # group B has no zeros: the zero-inflation of group B goes to zero
  d$y[d$g == "B" & d$y == 0] <- 1L
  d$gB <- as.numeric(d$g == "B")
  Z <- distributions7::zero_inflated(distributions7::poisson_distrib())
  f1 <- statmod(y ~ x | zi ~ jfp(~ g), distrib = Z, data = d)
  f2 <- statmod(y ~ x | zi ~ jfp(~ gB), distrib = Z, data = d)
  expect_equal(unlist(coef(f1), use.names = FALSE),
               unlist(coef(f2), use.names = FALSE), tolerance = 1e-8)
  expect_length(f1@aliased, 0L)
  expect_true(all(is.finite(unlist(coef(f1)))))
  expect_identical(statmod_certificate(f1)$state, "converged")
  # without the penalty the same coefficient runs away and is a boundary
  f0 <- statmod(y ~ x | zi ~ g, distrib = Z, data = d)
  expect_identical(statmod_certificate(f0)$state, "boundary")
})

test_that("jfp() beside a random intercept is fitted with reml()", {
  set.seed(4)
  d <- data.frame(id = factor(rep(1:20, each = 8)), x = rnorm(160),
                  z = rbinom(160, 1, 0.2))
  u <- rnorm(20, 0, 0.8)
  d$y <- rbinom(160, 1, plogis(-0.3 + 0.7 * d$x + u[d$id]))
  d$y[d$z == 1] <- 1L
  f <- statmod(y ~ x + jfp(~ z) + random(~ 1 | id),
               distrib = distributions7::bernoulli_distrib(), data = d)
  expect_true(is.finite(coef(f)$mu[["jfp.z"]]))
  expect_lt(coef(f)$mu[["jfp.z"]], 10)
  expect_identical(statmod_certificate(f)$state, "converged")
})

test_that("the consumers of a fit read the Jeffreys penalty", {
  d <- jfp_sep_data()
  f <- statmod(y ~ x2 + jfp(~ x1), distrib = distributions7::bernoulli_distrib(),
               data = d)
  for (vt in list(bayesian_vcov(), frequentist_vcov(), unconditional_vcov(),
                  frequentist_vcov(sandwich = hc0_sandwich()))) {
    V <- vcov(f, vcov_type = vt)
    expect_true(all(is.finite(diag(V))))
  }
  # the bayesian variance is the inverse of the penalized curvature
  sp <- f@spec
  des <- statmod_design(sp)
  K <- as_dense(statmod_information_at(sp, f@coefficients, des, TRUE)) +
    as_dense(statmod_penalty_at(sp, f@coefficients, f@hyper, des, "hessian"))
  expect_equal(unname(vcov(f, readable = FALSE)), unname(solve(K)),
               tolerance = 1e-6)
  s <- summary(f)
  expect_true(S7::S7_inherits(statmod_test(f, param = "mu",
                                           coefname = "jfp.x1",
                                           test = lr_test()), StatmodTest))
  pr <- statmod_profile(f, which = "mu:jfp.x1")
  expect_true(all(is.finite(pr@points$profile_penalty)))
  sim <- simulate(f, nsim = 1, seed = 1)
  expect_equal(nrow(sim), nrow(d))
  e <- f@edf
  expect_true(all(is.finite(e$edf)))
  expect_true(is.finite(AIC(f)))
})

test_that("the effective degrees of freedom can exceed p under separation", {
  d <- jfp_sep_data()
  f <- statmod(y ~ jfp(~ x1 + x2), distrib = distributions7::bernoulli_distrib(),
               data = d)
  expect_true(is.finite(sum(f@edf$edf)))
})

test_that("a block whose information is singular has no usable value", {
  d <- jfp_sep_data()
  sp <- statmod_spec(y ~ x1, distributions7::bernoulli_distrib(), d)
  des <- statmod_design(sp)
  # a column repeated in the block: the determinant is zero
  st <- jfp_state(sp, des, list(mu = c(0.1, 0.2)),
                  list(jfp_block("mu", c(2L, 2L))))
  expect_identical(jfp_value(st), -Inf)
  expect_error(jfp_gradient(st), "not positive definite")
})

test_that("lowrank_solve() equals the dense solve", {
  set.seed(2)
  p <- 30
  M <- crossprod(matrix(rnorm(60 * p), 60, p))
  V <- matrix(rnorm(p * 3), p, 3)
  u <- rnorm(p)
  ref <- solve(M + tcrossprod(V) + 0.5 * diag(p), u)
  expect_equal(lowrank_solve(M, V, u, damp = 0.5), ref, tolerance = 1e-10)
  Ms <- Matrix::Matrix(M, sparse = TRUE)
  expect_equal(lowrank_solve(Ms, V, u, damp = 0.5), ref, tolerance = 1e-10)
  # an indefinite A has no factor and the caller assembles the sum
  expect_null(lowrank_solve(M - 1e4 * diag(p), V, u))
  pc <- list(A = M, V = V)
  expect_equal(iwls_solve(pc, u, "chol")$delta, solve(M + tcrossprod(V), u),
               tolerance = 1e-10)
})

test_that("jfp() inside a subformula is rejected in this phase", {
  set.seed(1)
  d <- data.frame(x = runif(40), g = rnorm(40))
  d$y <- 2 * exp(-d$x) + rnorm(40, sd = 0.1)
  expect_error(statmod(y ~ 0 + nl(~ a * exp(-r * x), a ~ 1 + jfp(~ g)),
                       distrib = distributions7::gaussian1_distrib(),
                       data = d),
               "not implemented yet")
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
