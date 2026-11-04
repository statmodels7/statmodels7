# The shorthand gas(by = ~f) and the level: beside an intercept the level is
# left out of the development and held, without one it is developed; the
# filter is grouped by the random terms' variable; and a level is held only
# when the UNPENALIZED columns of its equation span the constant.

sim_panel <- function(S = 5, Tn = 30, seed = 3) {
  set.seed(seed)
  id <- factor(rep(seq_len(S), each = Tn))
  t <- rep(seq_len(Tn), S)
  x <- stats::rnorm(S * Tn)
  lev <- stats::rnorm(S)[as.integer(id)]
  y <- numeric(S * Tn)
  for (s in seq_len(S)) {
    ff <- 0
    for (i in which(id == s)) {
      mu <- 1 + lev[i] + 0.5 * x[i] + ff
      y[i] <- mu + stats::rnorm(1, sd = 0.7)
      ff <- 0.8 * ff + 0.2 * (y[i] - mu)
    }
  }
  data.frame(y, x, id, t)
}
dp <- sim_panel()
fit_g <- function(fm) statmod(fm, distributions7::gaussian1_distrib(), dp)
omega_names <- function(fit) grep("omega", names(coef(fit)$mu), value = TRUE)

test_that("beside an intercept the shorthand leaves the level held", {
  a <- fit_g(y ~ x + gas(p = 1, q = 1, time = t,
                         by = ~ random(~ 1 | id)))
  expect_identical(omega_names(a), "gas.omega")
  expect_identical(a@structural[[1L]]$held, "omega")
  expect_identical(coef(a)$mu[["gas.omega"]], 0)
  # the same model as the subformulas written for the other parameters
  b <- fit_g(y ~ x + gas(p = 1, q = 1, time = t,
                         kappa1 ~ random(~ 1 | id),
                         pacf1 ~ random(~ 1 | id), by = id))
  expect_equal(a@loglik, b@loglik, tolerance = 1e-8)
  # a level per cluster is a random intercept in the equation
  r <- fit_g(y ~ x + random(~ 1 | id) +
               gas(p = 1, q = 1, time = t, by = ~ random(~ 1 | id)))
  expect_identical(omega_names(r), "gas.omega")
  expect_identical(r@structural[[1L]]$held, "omega")
})

test_that("without an intercept the shorthand develops the level", {
  a <- fit_g(y ~ 0 + x + gas(p = 1, q = 1, time = t,
                             by = ~ random(~ 1 | id)))
  om <- omega_names(a)
  expect_true("gas.omega.(Intercept)" %in% om)
  expect_length(grep("^gas\\.omega\\.random", om), nlevels(dp$id))
  expect_length(a@structural[[1L]]$held, 0L)
  # the filter is grouped by id, as the spelling per parameter groups it
  b <- fit_g(y ~ 0 + x + gas(p = 1, q = 1, time = t,
                             omega ~ random(~ 1 | id),
                             kappa1 ~ random(~ 1 | id),
                             pacf1 ~ random(~ 1 | id), by = id))
  expect_equal(a@loglik, b@loglik, tolerance = 1e-8)
})

test_that("a subformula written for the level is kept beside an intercept", {
  e <- fit_g(y ~ x + gas(p = 1, q = 1, time = t, omega ~ random(~ 1 | id),
                         by = id))
  om <- omega_names(e)
  expect_length(grep("^gas\\.omega\\.random", om), nlevels(dp$id))
  expect_identical(e@structural[[1L]]$held, "omega.(Intercept)")
})

test_that("a random intercept does not hold the level", {
  held_of <- function(fm) {
    d <- statmod_design(statmod_spec(fm, distributions7::gaussian1_distrib(),
                                     dp))
    statmod_structural_state(d)$held[[1L]]
  }
  # penalized columns span the constant, but the penalty identifies them
  expect_length(held_of(y ~ 0 + x + random(~ 1 | id) +
                          gas(p = 1, q = 1, time = t, by = id)), 0L)
  expect_identical(held_of(y ~ x + random(~ 1 | id) +
                             gas(p = 1, q = 1, time = t, by = id)), "omega")
  # an unpenalized factor coded without an intercept still holds it
  expect_identical(held_of(y ~ 0 + id + gas(p = 1, q = 1, time = t,
                                            by = id)), "omega")
  # the two spellings of a level per cluster are one model: with the
  # intercept the level is held and the intercept is the population value;
  # without it the level carries the population value through the
  # recursion. Under ml() they reach the same maximum. Under reml() they do
  # not, because the intercept is integrated out as a fixed effect and the
  # level is not (misure_jfp/by_omega/RISULTATI.md).
  fit_ml <- function(fm) statmod(fm, distributions7::gaussian1_distrib(), dp,
                                 outer_criterion = ml())
  a <- fit_ml(y ~ x + random(~ 1 | id) + gas(p = 1, q = 1, time = t,
                                             by = id))
  b <- fit_ml(y ~ 0 + x + random(~ 1 | id) + gas(p = 1, q = 1, time = t,
                                                 by = id))
  expect_length(b@structural[[1L]]$held, 0L)
  cb <- coef(b)$mu
  expect_gt(abs(cb[["gas.omega"]]), 0)
  expect_equal(a@loglik, b@loglik, tolerance = 1e-7)
  expect_equal(coef(a)$mu[["(Intercept)"]],
               cb[["gas.omega"]] / (1 - cb[["gas.xi1"]]), tolerance = 1e-4)
})
