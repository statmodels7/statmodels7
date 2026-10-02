# A multivariate Student t prior over a block of random effects.
#
# Its Hessian in the coefficients moves with them, so the exact outer
# derivatives read the prior's third and fourth response derivatives per
# block (penalties7 0.30.0, from distributions7 0.67.0). Before those existed
# the criterion had no exact gradient, the search fell to the simplex, and on
# the random slopes of sleepstudy it stopped after 10 evaluations.

mvt_panel <- function(seed = 7L, m = 30L, ni = 8L) {
  set.seed(seed)
  g <- factor(rep(seq_len(m), each = ni))
  x <- rep(seq(-1, 1, length.out = ni), m)
  b <- matrix(stats::rnorm(2 * m, 0, 0.4), m, 2)
  b[1:3, ] <- c(2.5, -2.4, 2.8, 1.5, -1.8, 1.2)
  data.frame(g = g, x = x,
             y = 1 + 0.5 * x + b[as.integer(g), 1] + b[as.integer(g), 2] * x +
               stats::rnorm(m * ni, 0, 0.5))
}

mvt_prior <- distributions7::fixed(distributions7::mvstudent_t1_distrib(2),
                                   mu1 = 0, mu2 = 0)

test_that("the exact outer gradient over a multivariate t prior is exact", {
  skip_on_cran()
  dt <- mvt_panel()
  f <- y ~ x + random(~ x | g, distrib = mvt_prior)
  fam <- distributions7::gaussian1_distrib()
  fit0 <- statmod(f, fam, dt, outer_criterion = NULL)
  spec <- fit0@spec
  design <- statmod_design(spec)
  idx <- outer_hyper_index(spec, statmod_blocks(spec, design))
  for (mt in list(reml(hessian = "observed", marginal = "none"), aic())) {
    expect_true(outer_gradient_ok(spec, design, idx, mt, 1L))
    expect_true(outer_gradient_ok(spec, design, idx, mt, 2L))
    at <- function(eta) {
      hy <- eta_to_hyper(eta, idx, fit0@hyper)
      list(hy = hy, cf = fit_at_hyper(f, fam, dt, hy,
                                      polish = TRUE)$coefficients)
    }
    crit <- function(eta) {
      p <- at(eta)
      if (identical(mt@kind, "reml")) {
        statmod_marginal(spec, design, p$cf, p$hy, mt)$value
      } else {
        -statmod_pe(spec, design, p$cf, p$hy, mt)$value
      }
    }
    eta <- hyper_to_eta(fit0@hyper, idx) + 0.3
    p <- at(eta)
    ge <- if (identical(mt@kind, "reml")) {
      statmod_marginal_grad(spec, design, p$cf, p$hy, mt, idx)
    } else {
      -statmod_pe_derivs(spec, design, p$cf, p$hy, mt, idx, 1L)$grad
    }
    h <- 1e-3
    fd <- vapply(seq_along(eta), function(j) {
      e1 <- eta; e1[j] <- e1[j] + h
      e2 <- eta; e2[j] <- e2[j] - h
      (crit(e1) - crit(e2)) / (2 * h)
    }, numeric(1))
    # Measured 5.5e-07 under reml and 2.9e-06 under aic, a hundredth of the
    # reading at h = 1e-2: the difference's own truncation.
    expect_lt(max(abs(ge - fd) / abs(fd)), 2e-5)
  }
})

test_that("a multivariate t prior is searched on the gradient and converges", {
  skip_on_cran()
  fit <- statmod(y ~ x + random(~ x | g, distrib = mvt_prior),
                 distributions7::gaussian1_distrib(), mvt_panel())
  expect_identical(class(fit@methods$search)[[1L]],
                   class(optimizers7::lbfgs())[[1L]])
  expect_identical(statmod_certificate(fit)$state, "converged")
  # three groups far out in both effects: nu is small
  expect_lt(hyper(fit)$estimate[[4L]], 10)
})

test_that("a derivative-free default search is not given a stall rule", {
  skip_on_cran()
  # regime() has no exact outer gradient, so the default search is
  # derivative-free: brent() over the one hyperparameter here, the simplex
  # over more. The resolution rule reads a CHANGE in the objective, and the
  # best point of such a search does not move over an iteration that only
  # reshapes its bracket or simplex, so the rule would end the search there.
  set.seed(3)
  n <- 300L
  z <- stats::rnorm(n)
  st <- integer(n)
  st[1L] <- 1L
  for (t in 2:n) {
    st[t] <- if (stats::runif(1) < 0.05) 3L - st[t - 1L] else st[t - 1L]
  }
  d <- data.frame(z = z,
                  y = 0.8 * z + ifelse(st == 1L, -1.5, 1.5) +
                    stats::rnorm(n, sd = 0.8))
  fit <- statmod(y ~ ridge(~ 0 + z) + regime(k = 2),
                 distributions7::gaussian1_distrib(), d,
                 outer_criterion = reml())
  s <- fit@methods$search
  expect_identical(class(s)[[1L]], class(optimizers7::brent())[[1L]])
  expect_identical(class(s@criterion)[[1L]],
                   class(optimizers7::brent()@criterion)[[1L]])
})
