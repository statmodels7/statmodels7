# reml() integrates a free structural level as the constant that it adds to
# the predictor, a column of ones in its equation, as it integrates an
# intercept: the spelling with the intercept (the level held) and the
# spelling without it (the level free) give one criterion, and the exact
# outer gradient reads the matrix with the level's column appended.

rfl_panel <- function(S = 5, Tn = 30, seed = 3) {
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
rfl_d <- rfl_panel()
rfl_gau <- distributions7::gaussian1_distrib()

# the criterion and its exact gradient with the mode refitted at each point
rfl_harness <- function(form, data) {
  fit0 <- statmod(form, rfl_gau, data, outer_criterion = NULL)
  idx <- outer_hyper_index(fit0@spec,
                           statmod_blocks(fit0@spec,
                                          statmod_design(fit0@spec)))
  refit <- function(eta) {
    hy <- eta_to_hyper(eta, idx, fit0@hyper)
    a <- fit_at_hyper(form, rfl_gau, data, hy)
    list(a = a, hy = hy)
  }
  list(eta0 = hyper_to_eta(fit0@hyper, idx), idx = idx,
       fn = function(eta) {
         r <- refit(eta)
         statmod_marginal(r$a$spec, r$a$design, r$a$coefficients, r$hy,
                          reml())$value
       },
       gr = function(eta) {
         r <- refit(eta)
         statmod_marginal_grad(r$a$spec, r$a$design, r$a$coefficients, r$hy,
                               reml(), idx)
       },
       ok = function() {
         r <- refit(hyper_to_eta(fit0@hyper, idx))
         outer_gradient_ok(r$a$spec, r$a$design, idx, reml(), 1L)
       })
}
rfl_central <- function(fn, eta, h) {
  vapply(seq_along(eta), function(j) {
    ep <- eta
    em <- eta
    ep[j] <- ep[j] + h
    em[j] <- em[j] - h
    (fn(ep) - fn(em)) / (2 * h)
  }, numeric(1))
}

test_that("the two spellings of a level per cluster give one reml() fit", {
  a <- statmod(y ~ x + random(~ 1 | id) + gas(p = 1, q = 1, time = t,
                                              by = id), rfl_gau, rfl_d)
  b <- statmod(y ~ 0 + x + random(~ 1 | id) + gas(p = 1, q = 1, time = t,
                                                  by = id), rfl_gau, rfl_d)
  expect_identical(a@structural[[1L]]$held, "omega")
  expect_length(b@structural[[1L]]$held, 0L)
  expect_equal(a@loglik, b@loglik, tolerance = 1e-7)
  expect_equal(hyper(a)$estimate, hyper(b)$estimate, tolerance = 1e-5)
  cb <- coef(b)$mu
  expect_equal(coef(a)$mu[["(Intercept)"]],
               cb[["gas.omega"]] / (1 - cb[["gas.xi1"]]), tolerance = 1e-6)
})

test_that("the two spellings agree for regime() and for a joint filter", {
  skip_on_cran()
  dr <- rfl_d
  dr$y <- dr$y + rep(rep(c(0, 2), each = 5), 15)
  ra <- statmod(y ~ x + random(~ 1 | id) + regime(k = 2, by = id, time = t),
                rfl_gau, dr)
  rb <- statmod(y ~ 0 + x + random(~ 1 | id) +
                  regime(k = 2, by = id, time = t), rfl_gau, dr)
  expect_identical(ra@structural[[1L]]$held, "level1")
  expect_length(rb@structural[[1L]]$held, 0L)
  expect_equal(ra@loglik, rb@loglik, tolerance = 1e-7)
  expect_equal(hyper(ra)$estimate, hyper(rb)$estimate, tolerance = 1e-5)

  ja <- statmod(y ~ x + random(~ 1 | id) + gas(p = 1, q = 1, time = t,
                                               by = id) |
                  sigma ~ gas(p = 1, q = 1, time = t, by = id), rfl_gau, rfl_d)
  jb <- statmod(y ~ 0 + x + random(~ 1 | id) + gas(p = 1, q = 1, time = t,
                                                   by = id) |
                  sigma ~ 0 + gas(p = 1, q = 1, time = t, by = id),
                rfl_gau, rfl_d)
  expect_length(jb@structural[[1L]]$held, 0L)
  expect_equal(ja@loglik, jb@loglik, tolerance = 1e-7)
  expect_equal(hyper(ja)$estimate, hyper(jb)$estimate, tolerance = 1e-5)
})

test_that("the exact gradient reads the level's column", {
  skip_on_cran()
  h <- rfl_harness(y ~ 0 + x + random(~ 1 | id) +
                     gas(p = 1, q = 1, time = t, by = id), rfl_d)
  expect_true(h$ok())
  eta <- h$eta0 + 0.3
  g <- h$gr(eta)
  ref3 <- rfl_central(h$fn, eta, 1e-3)
  ref2 <- rfl_central(h$fn, eta, 1e-2)
  rel3 <- max(abs(g - ref3)) / max(abs(ref3))
  rel2 <- max(abs(g - ref2)) / max(abs(ref2))
  # measured 7.9e-07 and 7.8e-05: O(h^2), so nothing is missing
  expect_lt(rel3, 1e-5)
  expect_lt(rel3, rel2 / 10)
  # the negative control: the gradient of the criterion without the level's
  # column is 24 per cent away from the criterion with it
  local_mocked_bindings(level_augmented = function(...) NULL)
  g0 <- h$gr(eta)
  expect_gt(max(abs(g0 - ref3)) / max(abs(ref3)), 0.05)
})

test_that("the exact gradient reads the level of a joint filter", {
  skip_on_cran()
  # the level of mu free, the level of sigma held by the intercept
  h <- rfl_harness(y ~ 0 + x + random(~ 1 | id) +
                     gas(p = 1, q = 1, time = t, by = id) |
                     sigma ~ gas(p = 1, q = 1, time = t, by = id), rfl_d)
  expect_true(h$ok())
  eta <- h$eta0 + 0.15
  g <- h$gr(eta)
  ref <- rfl_central(h$fn, eta, 1e-3)
  expect_lt(max(abs(g - ref)) / max(abs(ref)), 1e-4)
})
