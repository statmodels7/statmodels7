# Where a fit begins, as a strategy rather than a vector.
#
# The first question is the one that matters most: the default must not have
# moved. Everything else is new machinery, and new machinery that changes an
# existing answer is a defect however good it looks.

set.seed(12)
n <- 120
dd <- data.frame(x = runif(n), z = runif(n))
dd$y <- 1 + 2 * dd$x + stats::rnorm(n, sd = 0.4)

test_that("the default is what it was, to the bit", {
  spec <- statmod_spec(y ~ x + z, distributions7::gaussian1_distrib(), dd)
  design <- statmod_design(spec)
  obj <- statmod_objective(spec, statmod_hyper_start(spec), design)
  a <- statmodels7:::statmod_start(spec, design, obj, NULL)
  b <- statmodels7:::statmod_start(spec, design, obj, start_intercepts())
  expect_identical(a, b)
  # and a fit begun either way is the same fit
  f1 <- statmod(y ~ x + z, distributions7::gaussian1_distrib(), dd)
  f2 <- statmod(y ~ x + z, distributions7::gaussian1_distrib(), dd,
                start = start_intercepts())
  expect_identical(f1@coefficients, f2@coefficients)
  expect_identical(f1@objective, f2@objective)
})

test_that("a list of values still means what it meant", {
  spec <- statmod_spec(y ~ x, distributions7::gaussian1_distrib(), dd)
  design <- statmod_design(spec)
  obj <- statmod_objective(spec, statmod_hyper_start(spec), design)
  v <- statmodels7:::statmod_start(spec, design, obj,
                                   list(mu = c(0.5, 1.5)))
  expect_equal(obj$split(v)$mu, c(0.5, 1.5))
  expect_error(statmodels7:::statmod_start(spec, design, obj,
                                           list(nope = 1)), "not a parameter")
  expect_error(statmodels7:::statmod_start(spec, design, obj,
                                           list(mu = 1)), "has length")
})

test_that("zeros are zeros and are not the default", {
  spec <- statmod_spec(y ~ x, distributions7::gaussian1_distrib(), dd)
  design <- statmod_design(spec)
  obj <- statmod_objective(spec, statmod_hyper_start(spec), design)
  z <- start_at(start_origin(), spec, design, obj)
  expect_true(all(unlist(z) == 0))
  d <- start_at(start_intercepts(), spec, design, obj)
  # the intercept-only fit puts the location near the response's own level,
  # which is the whole reason zero is not the default
  expect_gt(abs(d$mu[1L]), 1)
})

test_that("a random start is centred on the intercept fit unless told not to", {
  spec <- statmod_spec(y ~ x + z, distributions7::gaussian1_distrib(), dd)
  design <- statmod_design(spec)
  obj <- statmod_objective(spec, statmod_hyper_start(spec), design)
  base <- start_at(start_intercepts(), spec, design, obj)

  set.seed(3)
  a <- start_at(start_random(stats::rnorm, sd = 0.1), spec, design, obj)
  # near the intercept fit, and not equal to it
  expect_lt(abs(a$mu[1L] - base$mu[1L]), 1)
  expect_false(isTRUE(all.equal(a$mu, base$mu)))

  set.seed(3)
  b <- start_at(start_random(stats::rnorm, sd = 0.1, center = FALSE),
                spec, design, obj)
  expect_lt(abs(b$mu[1L]), 1)
  # the two differ by exactly the intercept fit, which is what centring means
  expect_equal(a$mu - b$mu, base$mu)

  # the stream is the caller's
  set.seed(5); p <- start_at(start_random(), spec, design, obj)
  set.seed(5); q <- start_at(start_random(), spec, design, obj)
  expect_equal(p, q)
})

test_that("a strategy drives a fit and the fit is still right", {
  # a random start must reach the same answer as the default on a convex
  # problem: if it did not, the fit had not converged
  ref <- statmod(y ~ x + z, distributions7::gaussian1_distrib(), dd)
  set.seed(7)
  f <- statmod(y ~ x + z, distributions7::gaussian1_distrib(), dd,
               start = start_random(stats::rnorm, sd = 0.5))
  expect_true(f@converged)
  expect_equal(f@coefficients$mu, ref@coefficients$mu, tolerance = 1e-6)
  g <- statmod(y ~ x + z, distributions7::gaussian1_distrib(), dd,
               start = start_origin())
  expect_true(g@converged)
  expect_equal(g@coefficients$mu, ref@coefficients$mu, tolerance = 1e-6)
})

test_that("a strategy that answers wrongly is refused where it is asked", {
  Bad <- S7::new_class("Bad", parent = start_strategy)
  start_at_g <- start_at
  S7::method(start_at_g, Bad) <- function(strategy, spec, design, obj, ...) {
    list(mu = c(1, 2, 3))
  }
  spec <- statmod_spec(y ~ x, distributions7::gaussian1_distrib(), dd)
  design <- statmod_design(spec)
  obj <- statmod_objective(spec, statmod_hyper_start(spec), design)
  expect_error(statmodels7:::statmod_start(spec, design, obj,
                                           Bad(label = "bad")),
               "one vector per distribution parameter")
})

test_that("the constructors validate and print", {
  expect_error(start_random(fn = "rnorm"), "must be a function")
  expect_error(start_random(center = NA), "TRUE or FALSE")
  expect_output(print(start_origin()), "zero")
  expect_output(print(start_random()), "random")
  expect_true(S7::S7_inherits(start_intercepts(), start_strategy))
})


# ---------------------------------------------------------------------------
# The intercept-only fit belongs to a PARAMETRIC intercept, and a term with
# parameters of its own is handed the response on the predictor's scale.
# ---------------------------------------------------------------------------

nl_start_data <- function(seed = 456, ni = 25, nt = 20) {
  set.seed(seed)
  u1 <- stats::rnorm(ni, 0, 5); u2 <- stats::rnorm(ni, 0, 1.2)
  u3 <- stats::rnorm(ni, 0, 0.3)
  d <- data.frame(id = factor(rep(seq_len(ni), each = nt)),
                  time = rep(0:(nt - 1), ni))
  p1 <- abs(50 + u1[d$id]); p2 <- abs(10 + u2[d$id]); p3 <- abs(2 + u3[d$id])
  d$y <- p1 / (1 + exp(-(d$time - p2) / p3)) + stats::rnorm(nrow(d), 0, 1.4)
  d
}
nl_start_formula <- function() {
  y ~ 0 + modelterms7::nl(
    ~ phi / (1 + exp(-(time - theta) / sigma)),
    phi ~ 1, theta ~ 1, sigma ~ 1,
    links = list(phi = linkfunctions7::log_link(),
                 theta = linkfunctions7::identity_link(),
                 sigma = linkfunctions7::log_link()))
}

test_that("an equation's intercept is a column of its parametric block", {
  spec <- statmod_spec(y ~ x + z, distributions7::gaussian1_distrib(), dd)
  design <- statmod_design(spec)
  expect_identical(statmodels7:::parametric_intercept(spec, design, "mu"), 1L)
  # a nonlinear term names the intercept of each of its own parameters the
  # same way, and those are not the equation's
  d2 <- nl_start_data()
  s2 <- statmod_spec(nl_start_formula(), distributions7::gaussian1_distrib(), d2)
  g2 <- statmod_design(s2)
  expect_true(endsWith(g2$mu$coef_names[1L], "(Intercept)"))
  expect_true(is.na(statmodels7:::parametric_intercept(s2, g2, "mu")))
})

test_that("the start of a nonlinear term is on the data's scale, not the response's", {
  d <- nl_start_data()
  spec <- statmod_spec(nl_start_formula(), distributions7::gaussian1_distrib(), d)
  design <- statmod_design(spec)
  obj <- statmod_objective(spec, statmod_hyper_start(spec, design), design)
  beta <- unlist(start_at(start_intercepts(), spec, design, obj),
                 use.names = FALSE)
  cf <- obj$split(beta)
  # phi rides a log link: writing mean(y) into its coefficient started it at
  # exp(23.9) = 2.5e10, with an objective of 7e20 and a gradient of 1.4e21
  expect_lt(exp(cf$mu[[1L]]), 10 * max(d$y))
  expect_equal(exp(cf$mu[[1L]]), 50, tolerance = 0.2)
  expect_true(is.finite(obj$fn(beta)))
  expect_lt(obj$fn(beta), 1e6)
  expect_lt(max(abs(obj$gr(beta))), 1e6)
})

test_that("the target exists for a mean and not for a scale", {
  d <- nl_start_data()
  spec <- statmod_spec(y ~ time, distributions7::gaussian1_distrib(), d)
  tg <- statmodels7:::predictor_target(spec, "mu")
  expect_equal(tg, d$y)                       # identity link on the mean
  expect_null(statmodels7:::predictor_target(spec, "sigma"))
})

test_that("the target is carried onto the predictor's scale", {
  set.seed(5)
  d <- data.frame(x = seq(0, 5, length.out = 200))
  d$y <- stats::rpois(200, exp(1 + 0.3 * d$x))
  spec <- statmod_spec(y ~ x, distributions7::poisson_distrib(), d)
  tg <- statmodels7:::predictor_target(spec, "mu")
  expect_true(all(is.finite(tg)))
  # a count of zero sits on the bound of the mean's domain and is moved half
  # way to the smallest admissible value, never differenced past it
  expect_true(all(exp(tg) > 0))
  expect_equal(tg[d$y > 0], log(d$y[d$y > 0]))
})

test_that("an ordinary model's start is untouched", {
  spec <- statmod_spec(y ~ x + z, distributions7::gaussian1_distrib(), dd)
  design <- statmod_design(spec)
  obj <- statmod_objective(spec, statmod_hyper_start(spec, design), design)
  beta <- unlist(start_at(start_intercepts(), spec, design, obj),
                 use.names = FALSE)
  cf <- obj$split(beta)
  eta0 <- statmod_intercepts(spec)
  expect_equal(cf$mu[[1L]], eta0[["mu"]])
  expect_true(all(cf$mu[-1L] == 0))
})
