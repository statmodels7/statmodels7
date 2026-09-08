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


# --- start_from(): another fit's estimates ------------------------------

sf <- local({
  set.seed(31)
  m <- 400
  d <- data.frame(x = stats::runif(m, -1, 1), z = stats::runif(m, -1, 1),
                  w = stats::runif(m, -1, 1))
  d$y <- stats::rpois(m, exp(0.4 + 0.6 * d$x - 0.3 * d$z))
  d$ys <- sin(2 * pi * d$x) + stats::rnorm(m, 0, 0.3)
  list(d = d,
       full = statmod(y ~ x + z + w, distributions7::poisson_distrib(), d),
       sm6 = statmod(ys ~ s(x, k = 6), distributions7::gaussian1_distrib(),
                     d))
})

test_that("the parametric block is matched column by column", {
  # adding or dropping a covariate is what this exists for: the columns the
  # two models share carry their estimates across and the rest falls back
  spec <- statmod_spec(y ~ z + w, distributions7::poisson_distrib(), sf$d)
  design <- statmod_design(spec)
  s <- start_at(start_from(sf$full), spec, design, NULL)
  took <- attr(s, "taken")

  expect_setequal(took$coefficient, c("(Intercept)", "z", "w"))
  expect_identical(unique(took$parameter), "mu")
  # by NAME and not by position: 'z' is the third column of the reference and
  # the second here, so a positional copy would put 'x' where 'z' belongs
  ref <- sf$full@coefficients$mu
  expect_equal(s$mu, c(ref[[1L]], ref[[3L]], ref[[4L]]))
})

test_that("a start does not move where the fit lands", {
  # the property that matters most: on a convex problem the answer is the
  # data's, not the starting point's, and a strategy that changed it would be
  # a defect however fast it was
  a <- statmod(y ~ z + w, distributions7::poisson_distrib(), sf$d)
  b <- statmod(y ~ z + w, distributions7::poisson_distrib(), sf$d,
               start = start_from(sf$full))
  expect_equal(unlist(a@coefficients), unlist(b@coefficients),
               tolerance = 1e-6)
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)),
               tolerance = 1e-8)
})

test_that("a parameter the reference does not have falls back", {
  # a negative binomial started from a Poisson: the mean's coefficients are
  # taken and 'theta', which the reference has no counterpart for, is left to
  # the fallback strategy
  spec <- statmod_spec(y ~ x + z, distributions7::negbin2_distrib(), sf$d)
  design <- statmod_design(spec)
  s <- start_at(start_from(sf$full), spec, design, NULL)
  took <- attr(s, "taken")

  expect_identical(unique(took$parameter), "mu")
  expect_false("theta" %in% took$parameter)
  base <- start_at(start_intercepts(), spec, design, NULL)
  expect_equal(s$theta, base$theta)
  expect_length(s$theta, design$theta$npar)
})

test_that("the same basis is carried across and not re-estimated", {
  # where the coefficients mean the same coordinates they are copied, which is
  # exact and costs nothing; re-estimating them would be work for no answer
  sp6 <- statmod_spec(ys ~ s(x, k = 6),
                      distributions7::gaussian1_distrib(), sf$d)
  s6 <- start_at(start_from(sf$sm6), sp6, statmod_design(sp6), NULL)
  tk <- attr(s6, "taken")
  expect_true(all(c("s(x).lin", "s(x).z1", "s(x).z4") %in% tk$coefficient))
  expect_identical(unique(tk$how), "matched")
  expect_equal(s6$mu, sf$sm6@coefficients$mu)
})

test_that("a basis of another dimension is projected, not copied", {
  # s(x, k = 6) and s(x, k = 10) give the same names to different coordinates
  # -- measured, z1 is 0.0267 against -0.0277, OPPOSITE IN SIGN -- so what
  # carries across is the fitted FUNCTION and not the numbers. The two blocks
  # do not even meet by name: a block is keyed by the term's own call, so
  # 's(x, k = 6)' and 's(x, k = 10)' are two keys, and what pairs them is the
  # stem of their coefficient names.
  sp10 <- statmod_spec(ys ~ s(x, k = 10),
                       distributions7::gaussian1_distrib(), sf$d)
  de10 <- statmod_design(sp10)
  s10 <- start_at(start_from(sf$sm6), sp10, de10, NULL)
  tk <- attr(s10, "taken")

  sm <- tk[startsWith(tk$coefficient, "s(x)"), , drop = FALSE]
  expect_gt(nrow(sm), 0)
  expect_identical(unique(sm$how), "projected")

  # NOT the reference's numbers: a copy by name would put z1 in with the
  # wrong sign, which is why the coefficients are not comparable
  ref_z1 <- sf$sm6@coefficients$mu[[
    match("s(x).z1", statmod_design(sf$sm6@spec)$mu$coef_names)]]
  new_z1 <- s10$mu[[match("s(x).z1", de10$mu$coef_names)]]
  expect_false(isTRUE(all.equal(new_z1, ref_z1)))

  # what IS reproduced is the function: the starting predictor sits on the
  # reference's, where the fallback is a whole standard deviation away
  eta_ref <- stats::predict(sf$sm6, "link")$mu
  base <- start_at(start_intercepts(), sp10, de10, NULL)
  rmse <- function(s) sqrt(mean((as.numeric(de10$mu$X %*% s$mu) - eta_ref)^2))
  expect_lt(rmse(s10), 0.05)
  expect_gt(rmse(base), 0.3)
})

test_that("the projection is exact where the coarse basis is nested", {
  # equally spaced interior knots make s(x, k = 6) a SUBSPACE of s(x, k = 12)
  # -- three intervals refined into nine -- and there the projection
  # reproduces the function to machine precision. k = 10 is not a refinement
  # of k = 6, so there it is a genuine approximation, and the contrast is what
  # says the exactness is the geometry's and not the harness's.
  gauss <- distributions7::gaussian1_distrib()
  eta_ref <- stats::predict(sf$sm6, "link")$mu
  gap <- function(k) {
    sp <- statmod_spec(stats::as.formula(sprintf("ys ~ s(x, k = %d)", k)),
                       gauss, sf$d)
    de <- statmod_design(sp)
    s <- start_at(start_from(sf$sm6), sp, de, NULL)
    sqrt(mean((as.numeric(de$mu$X %*% s$mu) - eta_ref)^2))
  }
  expect_lt(gap(12), 1e-8)
  expect_gt(gap(10), 100 * gap(12))
})

test_that("a reference read at other rows is not projected", {
  # the projection compares two predictors observation by observation, so an
  # identical response is what says the two designs are read at the same rows;
  # without that the comparison is between functions evaluated at different
  # points, which means nothing
  d2 <- sf$d
  d2$ys <- d2$ys + 1
  other <- statmod(ys ~ s(x, k = 6), distributions7::gaussian1_distrib(), d2)
  sp10 <- statmod_spec(ys ~ s(x, k = 10),
                       distributions7::gaussian1_distrib(), sf$d)
  s <- start_at(start_from(other), sp10, statmod_design(sp10), NULL)
  expect_false(any(attr(s, "taken")$how == "projected"))
})

test_that("a projected start does not move where the fit lands", {
  # the same property the matched blocks are held to: a starting strategy that
  # changed the answer would be a defect however close it started
  gauss <- distributions7::gaussian1_distrib()
  a <- statmod(ys ~ s(x, k = 10), gauss, sf$d)
  b <- statmod(ys ~ s(x, k = 10), gauss, sf$d, start = start_from(sf$sm6))
  expect_equal(unlist(a@coefficients), unlist(b@coefficients),
               tolerance = 1e-5)
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)),
               tolerance = 1e-7)
})

test_that("the fallback strategy is the caller's", {
  # I(x^2) is a column the reference does not carry, so it is the fallback
  # that answers for it: start_origin() leaves it at zero where the default
  # would too, and the intercept shows the two strategies apart
  spec <- statmod_spec(y ~ x + I(x^2), distributions7::poisson_distrib(),
                       sf$d)
  design <- statmod_design(spec)
  s <- start_at(start_from(sf$full, rest = start_origin()), spec, design,
                NULL)
  expect_equal(s$mu[[3L]], 0)
  expect_equal(s$mu[[2L]], sf$full@coefficients$mu[[2L]])
  expect_equal(s$mu[[1L]], sf$full@coefficients$mu[[1L]])

  # and with the default fallback the intercept is still the reference's,
  # since a taken coefficient overwrites whatever the fallback put there
  s2 <- start_at(start_from(sf$full), spec, design, NULL)
  expect_equal(s2$mu[[1L]], sf$full@coefficients$mu[[1L]])
  expect_equal(s2$mu[[3L]], 0)
})

test_that("start_from returns a strategy that satisfies the contract", {
  spec <- statmod_spec(y ~ z + w, distributions7::poisson_distrib(), sf$d)
  design <- statmod_design(spec)
  s <- start_at(start_from(sf$full), spec, design, NULL)
  expect_true(S7::S7_inherits(start_from(sf$full), start_strategy))
  expect_setequal(names(s), spec@distrib@params)
  for (p in spec@distrib@params) {
    expect_length(s[[p]], design[[p]]$npar)
    expect_true(is.numeric(s[[p]]) && !anyNA(s[[p]]))
  }
})

test_that("start_from refuses what it cannot use", {
  expect_error(start_from(42), "must be a statmod fit")
  expect_error(start_from(sf$full, rest = "intercepts"),
               "must be a start strategy")
  # a chain of references would be resolved in some order and the order would
  # decide the answer
  expect_error(start_from(sf$full, rest = start_from(sf$full)),
               "cannot itself be")
})


# Where the HYPERPARAMETERS begin ------------------------------------------

test_that("a hyperparameter with a finite bound begins where it did", {
  # the rule is one unit inside whichever ends are finite, and this is the
  # half of it that must not move: a ridge, a smooth and an ordinary random
  # effect all have a hyperparameter bounded below
  expect_equal(unname(penalty_theta_start(penalties7::ridge_penalty(5L))), 1)
  expect_equal(unname(penalty_theta_start(penalties7::lasso_penalty(5L))), 1)
  expect_equal(unname(penalty_theta_start(penalties7::quadratic_penalty(diag(4)))), 1)
  # both ends finite: the midpoint. The elastic net's alpha lives on [0, 1]
  en <- penalty_theta_start(penalties7::elasticnet_penalty(5L))
  expect_equal(unname(en[["alpha"]]), 0.5)
  expect_equal(unname(en[["lambda"]]), 1)
  # bounded below at something other than zero
  expect_equal(unname(penalty_theta_start(penalties7::scad_penalty(5L))[["a"]]), 3)
})

test_that("a free coordinate of a chart begins at the chart's neutral point", {
  # A coordinate unbounded on BOTH sides is not a scale a hyperparameter lives
  # on, it is a free coordinate of a parameters7 chart, and one is not one
  # there: on dr_prod(2) it reads as standard deviations of e and a
  # correlation of -0.664, where zero reads as unit standard deviations and no
  # correlation -- which is what an ordinary random effect starts from.
  pen <- penalties7::structured_penalty(parameters7::dr_prod(2L))
  b <- pen@params_bounds
  expect_true(all(vapply(b, function(z) !is.finite(z[1L]) && !is.finite(z[2L]),
                         logical(1))))
  s <- penalty_theta_start(pen)
  expect_equal(unname(s), rep(0, length(b)))

  # WHAT IT MEANS, read off the chart rather than asserted about the number
  S0 <- parameters7::param_value(parameters7::dr_prod(2L), unlist(s))
  expect_equal(unname(sqrt(diag(S0))), c(1, 1), tolerance = 1e-12)
  expect_equal(S0[1, 2] / sqrt(S0[1, 1] * S0[2, 2]), 0, tolerance = 1e-12)

  # THE NEGATIVE CONTROL, without which the test above passes for any rule
  # returning zeros: the value it replaces means something else entirely
  S1 <- parameters7::param_value(parameters7::dr_prod(2L),
                                 rep(1, length(b)))
  expect_equal(sqrt(S1[1, 1]), exp(1), tolerance = 1e-6)
  expect_lt(S1[1, 2] / sqrt(S1[1, 1] * S1[2, 2]), -0.6)
})
