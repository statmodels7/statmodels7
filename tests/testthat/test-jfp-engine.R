# The Jeffreys engine: value, gradient and Hessian of a log|I_tilde| over the
# joint vector, against numDeriv applied to an assembly of the block's
# expected information written here by hand.

jfp_split <- function(design, u) {
  npar <- vapply(design, function(z) z$npar, integer(1))
  offs <- cumsum(npar) - npar
  stats::setNames(lapply(seq_along(npar), function(a)
    u[offs[a] + seq_len(npar[a])]), names(design))
}

# the reference: predictors from statmod_eta(), the block's matrix assembled
# from distrib_expected_hessian() and its log-determinant by determinant()
jfp_reference <- function(spec, design, blocks) function(u) {
  cf <- jfp_split(design, u)
  th <- statmod_eta(spec, design, cf)$theta
  dz <- statmod_design_at(spec, cf, design)
  E <- distributions7::distrib_expected_hessian(spec@distrib, spec@response,
                                                th, scale = "link")
  tot <- 0
  for (b in blocks) {
    X <- as.matrix(dz[[b$param]]$X)[, b$cols, drop = FALSE]
    if (b$const) X <- cbind(1, X)
    e <- -rep_len(E[[paste(b$param, b$param, sep = "_")]], spec@n_obs)
    I <- crossprod(X * (e * spec@weights), X)
    tot <- tot + b$power * as.numeric(determinant(I, logarithm = TRUE)$modulus)
  }
  tot
}

jfp_check_case <- function(formula, family, data, blocks, u0) {
  sp <- statmod_spec(formula, family, data)
  des <- statmod_design(sp)
  ref <- jfp_reference(sp, des, blocks)
  st <- jfp_state(sp, des, jfp_split(des, u0), blocks)
  h <- jfp_hessian(st)
  hn <- numDeriv::hessian(ref, u0)
  list(value = c(jfp_value(st), ref(u0)),
       grad = c(max(abs(jfp_gradient(st) - numDeriv::grad(ref, u0))) /
                  max(abs(numDeriv::grad(ref, u0)))),
       hess = max(abs(h$hessian - hn)) / max(abs(hn)),
       parts = max(abs(as.matrix(h$sparse) - tcrossprod(h$lowrank) -
                         h$hessian)),
       sparse_only = max(abs(as.matrix(h$sparse) - hn)) / max(abs(hn)))
}

jfp_data <- function() {
  set.seed(7)
  n <- 120
  d <- data.frame(x = rnorm(n), z = rnorm(n), w = runif(n),
                  g = factor(rep(1:10, 12)),
                  h = factor(sample(c("a", "b", "c"), n, TRUE)))
  d$yb <- rbinom(n, 1, plogis(0.3 + d$x - 0.5 * d$z))
  d$ynb <- rnbinom(n, mu = exp(1 + 0.3 * d$x), size = exp(0.5 + 0.5 * d$x))
  d$yzip <- ifelse(runif(n) < plogis(-1 + (d$h == "b")), 0,
                   rpois(n, exp(1 + 0.3 * d$x)))
  d$yt <- 1 + d$x + rt(n, df = 5)
  d$ynl <- rpois(n, exp(0.2 * d$z) * 3 * exp(-0.5 * d$w))
  d
}

test_that("the engine's derivatives match numDeriv on several families", {
  skip_if_not_installed("numDeriv")
  d <- jfp_data()
  set.seed(3)
  cases <- list(
    list(yb ~ x + z + random(~1 | g), distributions7::bernoulli_distrib(),
         list(jfp_block("mu", 2:3, const = TRUE))),
    list(yb ~ x + z + random(~1 | g), distributions7::bernoulli_distrib(
      link_mu = linkfunctions7::probit_link()),
      list(jfp_block("mu", 2:3, const = FALSE))),
    list(ynb ~ x + random(~1 | g) | theta ~ x,
         distributions7::negbin2_distrib(),
         list(jfp_block("mu", 2L, const = TRUE),
              jfp_block("theta", 1:2, power = 1))),
    list(yzip ~ x + random(~1 | g) | zi ~ h,
         distributions7::zero_inflated(distributions7::poisson_distrib()),
         list(jfp_block("zi", 2:3, const = TRUE))),
    list(yt ~ x + random(~1 | g) | nu ~ 1,
         distributions7::student_t1_distrib(),
         list(jfp_block("mu", 2L, const = TRUE), jfp_block("nu", 1L)))
  )
  for (cs in cases) {
    sp <- statmod_spec(cs[[1]], cs[[2]], d)
    m <- sum(vapply(statmod_design(sp), function(z) z$npar, integer(1)))
    r <- jfp_check_case(cs[[1]], cs[[2]], d, cs[[3]], rnorm(m, sd = 0.2))
    expect_equal(r$value[1], r$value[2], tolerance = 1e-12)
    expect_lt(r$grad, 1e-6)
    expect_lt(r$hess, 1e-5)
    expect_lt(r$parts, 1e-10)
    # the negative control: the Hessian without its low-rank part is wrong
    expect_gt(r$sparse_only, 1e-3)
  }
})

test_that("the engine's Hessian reads a refreshed block elsewhere in the model", {
  skip_if_not_installed("numDeriv")
  d <- jfp_data()
  r <- jfp_check_case(ynl ~ 0 + z + nl(~ a * exp(-r * w)),
                      distributions7::poisson_distrib(), d,
                      list(jfp_block("mu", 1L, const = TRUE)),
                      c(0.1, 1, 0.5))
  expect_lt(r$grad, 1e-7)
  expect_lt(r$hess, 1e-7)
  expect_gt(r$sparse_only, 1e-2)
})

test_that("Newton on the engine reaches Firth's estimate and the Poisson identity", {
  set.seed(11)
  n <- 60
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  # x1 separates the outcome quasi-completely
  d$y <- as.numeric(d$x1 > 0)
  d$y[1:6] <- 1 - d$y[1:6]
  d$y[d$x1 > 1] <- 1
  d$cnt <- rpois(n, exp(0.2 + 0.5 * d$x1))
  fit <- function(sp, bl) {
    des <- statmod_design(sp)
    u <- numeric(sum(vapply(des, function(z) z$npar, integer(1))))
    for (it in 1:100) {
      cf <- jfp_split(des, u)
      st <- jfp_state(sp, des, cf, bl)
      g <- unlist(statmod_score_at(sp, cf, des)) + jfp_gradient(st)
      if (max(abs(g)) < 1e-11) break
      H <- -as.matrix(statmod_information_at(sp, cf, des, expected = FALSE)) +
        jfp_hessian(st)$hessian
      u <- u - solve(H, g)
    }
    list(u = u, mu = statmod_eta(sp, des, jfp_split(des, u))$theta$mu,
         grad = max(abs(g)))
  }
  # Firth on every coefficient: the closed form's gradient vanishes there
  sp <- statmod_spec(y ~ x1 + x2, distributions7::bernoulli_distrib(), d)
  X <- statmod_design(sp)$mu$X
  expect_true(jfp_constant(NULL, X[, 1, drop = FALSE], X[, 2:3]))
  f <- fit(sp, list(jfp_block("mu", 2:3, const = TRUE)))
  closed <- function(b) {
    mu <- plogis(drop(X %*% b))
    sum(dbinom(d$y, 1, mu, log = TRUE)) +
      0.5 * as.numeric(determinant(crossprod(X * (mu * (1 - mu)), X))$modulus)
  }
  expect_lt(f$grad, 1e-11)
  skip_if_not_installed("numDeriv")
  expect_lt(max(abs(numDeriv::grad(closed, f$u))), 1e-8)
  # the Poisson identity: sum(mu) - sum(y) = a p_tilde
  sp <- statmod_spec(cnt ~ x1 + x2, distributions7::poisson_distrib(), d)
  for (a in c(0.5, 2)) {
    f <- fit(sp, list(jfp_block("mu", 2:3, power = a, const = TRUE)))
    expect_equal(sum(f$mu) - sum(d$cnt), 3 * a, tolerance = 1e-10)
  }
  # the negative control: without an intercept there is no such identity
  sp0 <- statmod_spec(cnt ~ 0 + x1 + x2, distributions7::poisson_distrib(), d)
  f0 <- fit(sp0, list(jfp_block("mu", 1:2, const = TRUE)))
  expect_gt(abs(sum(f0$mu) - sum(d$cnt) - 1.5), 0.1)
})

test_that("jfp_constant() applies the rule of the intercept argument", {
  one <- matrix(1, 10, 1)
  x <- matrix(rnorm(10), 10, 1)
  f <- cbind(rep(c(1, 0), 5), rep(c(0, 1), 5))
  expect_true(jfp_constant(NULL, one, x))
  expect_false(jfp_constant(NULL, NULL, x))
  expect_true(jfp_constant(TRUE, NULL, x))
  expect_false(jfp_constant(FALSE, one, x))
  # a factor with all its levels spans the constant
  expect_true(jfp_constant(NULL, f, x))
  expect_false(jfp_constant(TRUE, one, cbind(f, x)))
  expect_error(jfp_constant(NA, one, x), "NULL, TRUE or FALSE")
})

test_that("the engine rejects what it does not compute", {
  expect_error(jfp_block("mu", 1L, power = 0), "positive")
  expect_error(jfp_block("mu", 1L, power = c(1, 2)), "positive")
  d <- jfp_data()
  sp <- statmod_spec(yb ~ x, distributions7::bernoulli_distrib(), d)
  des <- statmod_design(sp)
  cf <- list(mu = c(0, 0))
  expect_error(jfp_state(sp, des, cf, jfp_block("mu", 2L, kind = "refresh")),
               "not\n  implemented yet")
  expect_error(jfp_state(sp, des, cf, jfp_block("mu", 3L)), "positions")
  expect_error(jfp_state(sp, des, cf, jfp_block("sigma", 1L)), "no parameter")
})
