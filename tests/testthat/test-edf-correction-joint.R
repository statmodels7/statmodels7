# What estimating a hyperparameter cost, priced on the vector the mode
# actually moves in. For a model carrying a filter that is the JOINT vector,
# the coefficients followed by the term's own parameters, and a penalty over
# those parameters is a column of no design -- so read on the coefficients
# alone the movement is identically zero in exactly the coordinates such a
# penalty shrinks.
#
# The checks are built so that the route this replaced FAILS them: the
# coefficient-space reading is computed here beside the shipped one and
# asserted to be exactly zero, so restoring it cannot pass.

edf_panel <- function(m = 10L, ni = 40L, seed = 31L) {
  set.seed(seed)
  u <- stats::rnorm(m, sd = 0.6)
  g <- factor(rep(seq_len(m), each = ni))
  y <- unlist(lapply(seq_len(m), function(i) {
    f <- 0
    out <- numeric(ni)
    for (t in seq_len(ni)) {
      # the observation is drawn INSIDE the loop and its residual drives the
      # recursion, which is the score for a gaussian. Driving it with the
      # noiseless predictor instead leaves f identically zero and the panel
      # has no dynamics at all.
      mu <- 0.4 + f
      out[t] <- stats::rnorm(1, mu, 1)
      f <- 0.25 * exp(u[i]) * (out[t] - mu) + 0.6 * f
    }
    out
  }))
  data.frame(y = y, g = g)
}

test_that("a penalty over a filter's own parameters is priced at all", {
  skip_on_cran()
  d <- edf_panel()
  fit <- statmod(y ~ gas(p = 1, q = 1, by = g, alpha1 ~ 1 + random(~ 1 | g)),
                 gaussian1_distrib(), d, outer_criterion = reml())
  dz <- statmod_design(fit@spec)
  idx <- outer_hyper_index(fit@spec, statmod_blocks(fit@spec, dz))

  # the premise: one hyperparameter, estimated, and the fit moved off its
  # starting value of 1 -- a fit parked at the start has no uncertainty to
  # propagate and would pass a positivity check for the wrong reason
  expect_equal(nrow(idx), 1L)
  expect_gt(abs(log(hyper(fit)$estimate[[1L]])), 0.5)

  cc <- statmod_edf_correction(fit@spec, fit@coefficients, fit@hyper, dz,
                               fit@methods$outer)
  expect_equal(cc$n_hyper, 1L)
  # measured 1.108889 on this panel against a total edf of 4.4976. The
  # assertion is on the RELATION rather than the value, a fitted quantity
  # being platform arithmetic: it is a substantial fraction of the count and
  # is nowhere near the zero the coefficient-space reading gives.
  expect_gt(cc$total, 0.5)
  expect_lt(cc$total, 2)
})

test_that("the route it replaced returns exactly zero on the same fit", {
  skip_on_cran()
  # THE INJECTION. This is the coefficient-space reading, computed here, and
  # it is what shipped before: hyper_mode_cross() skips a structural penalty
  # where the matrix spans the coefficients alone, so every column of the
  # movement is zero and so is the correction. Restoring that route makes the
  # test above fail and this one pass with the two agreeing, which is what
  # says the checks discriminate.
  d <- edf_panel()
  fit <- statmod(y ~ gas(p = 1, q = 1, by = g, alpha1 ~ 1 + random(~ 1 | g)),
                 gaussian1_distrib(), d, outer_criterion = reml())
  dz <- statmod_design(fit@spec)
  idx <- outer_hyper_index(fit@spec, statmod_blocks(fit@spec, dz))

  H <- statmod_information_at(fit@spec, fit@coefficients, dz, TRUE, "opg")
  nb <- nrow(as.matrix(H))
  cr0 <- hyper_mode_cross(fit@spec, dz, fit@coefficients, fit@hyper, idx, nb)
  expect_equal(cr0$skipped, 1L)
  expect_equal(max(abs(cr0$cross)), 0)

  nj <- nrow(as.matrix(statmod_full_information(fit@spec, fit@coefficients,
                                                dz)))
  cr1 <- hyper_mode_cross(fit@spec, dz, fit@coefficients, fit@hyper, idx, nj,
                          joint = TRUE)
  expect_equal(cr1$skipped, 0L)
  expect_gt(max(abs(cr1$cross)), 1)
})

test_that("the joint correction is the two matrices the criterion reads", {
  skip_on_cran()
  # neither matrix is assembled inside statmod_edf_correction():
  # statmod_marginal_full() is the one place K + S is built on the joint
  # vector and statmod_full_information() the one place K is. Assembling the
  # contraction here from those two must reproduce what the function returns,
  # which is what says it reads them rather than a second copy.
  d <- edf_panel()
  fit <- statmod(y ~ gas(p = 1, q = 1, by = g, alpha1 ~ 1 + random(~ 1 | g)),
                 gaussian1_distrib(), d, outer_criterion = reml())
  dz <- statmod_design(fit@spec)
  idx <- outer_hyper_index(fit@spec, statmod_blocks(fit@spec, dz))

  M <- as.matrix(statmod_marginal_full(fit@spec, dz, fit@coefficients,
                                       fit@hyper))
  Hj <- as.matrix(statmod_full_information(fit@spec, fit@coefficients, dz))
  cr <- hyper_mode_cross(fit@spec, dz, fit@coefficients, fit@hyper, idx,
                         nrow(Hj), joint = TRUE)
  Ho <- statmod_marginal_hess(fit@spec, dz, fit@coefficients, fit@hyper,
                              fit@methods$outer, idx, NULL)
  Vth <- solve(-as.matrix(Ho))
  J <- -solve(M) %*% cr$cross
  ref <- sum((J %*% Vth %*% t(J)) * t(Hj))

  cc <- statmod_edf_correction(fit@spec, fit@coefficients, fit@hyper, dz,
                               fit@methods$outer)
  expect_equal(cc$total, ref, tolerance = 1e-10)
})

test_that("a model with no filter takes the branch that was there before", {
  skip_on_cran()
  # by CONSTRUCTION and not by tolerance: statmod_marginal_full() returns
  # NULL where the design carries no structural term of the filter shape, so
  # `joint` is FALSE and the executed lines are the ones that shipped. The
  # coefficient-space assembly is computed here and must agree exactly.
  set.seed(4L)
  n <- 600L
  x <- stats::runif(n)
  g <- factor(rep(1:20, length.out = n))
  y <- sin(2 * pi * x) + stats::rnorm(20, sd = 0.5)[g] +
    stats::rnorm(n, sd = 0.4)
  d <- data.frame(y = y, x = x, g = g)
  fit <- statmod(y ~ ridge(~ x + I(x^2) + I(x^3)) + random(~ 1 | g),
                 gaussian1_distrib(), d, outer_criterion = reml())
  dz <- statmod_design(fit@spec)

  expect_null(statmod_marginal_full(fit@spec, dz, fit@coefficients,
                                    fit@hyper))

  idx <- outer_hyper_index(fit@spec, statmod_blocks(fit@spec, dz))
  H <- statmod_information_at(fit@spec, fit@coefficients, dz, TRUE, "opg")
  S <- zap_nonfinite(statmod_penalty_at(fit@spec, fit@coefficients,
                                        fit@hyper, dz, "hessian"))
  J <- -solve(H + S) %*%
    hyper_mode_cross(fit@spec, dz, fit@coefficients, fit@hyper, idx,
                     nrow(as.matrix(H)))$cross
  Ho <- statmod_marginal_hess(fit@spec, dz, fit@coefficients, fit@hyper,
                              fit@methods$outer, idx, NULL)
  ref <- sum((J %*% solve(-as.matrix(Ho)) %*% t(J)) *
               t(as.matrix(H)))

  cc <- statmod_edf_correction(fit@spec, fit@coefficients, fit@hyper, dz,
                               fit@methods$outer)
  expect_gt(cc$total, 0)
  expect_equal(cc$total, ref, tolerance = 1e-10)
})

test_that("a term of the likelihood shape is not claimed by the joint route", {
  skip_on_cran()
  # regime() contributes a mixed LIKELIHOOD rather than a predictor, so it
  # answers none of the filter rungs and carries no penalty over its own
  # parameters -- it takes no subformula. statmod_marginal_full() filters on
  # the filter shape, so the branch does not claim it.
  set.seed(7L)
  n <- 500L
  Z <- matrix(stats::rnorm(n * 4), n, 4)
  colnames(Z) <- paste0("z", 1:4)
  st <- integer(n)
  st[1L] <- 1L
  for (t in 2:n) {
    st[t] <- if (stats::runif(1) < 0.05) 3L - st[t - 1L] else st[t - 1L]
  }
  y <- as.numeric(Z %*% c(0.9, -0.6, 0, 0)) +
    ifelse(st == 1L, -1.5, 1.5) + stats::rnorm(n, sd = 0.8)
  d <- data.frame(y = y)
  d$Z <- Z
  fit <- statmod(y ~ ridge(Z) + regime(k = 2), gaussian1_distrib(), d,
                 outer_criterion = reml())
  dz <- statmod_design(fit@spec)

  expect_true(length(attr(dz, "structural")) > 0L)
  expect_null(statmod_marginal_full(fit@spec, dz, fit@coefficients,
                                    fit@hyper))
  expect_equal(length(modelterms7::term_penalties(regime(k = 2))), 0L)

  cc <- statmod_edf_correction(fit@spec, fit@coefficients, fit@hyper, dz,
                               fit@methods$outer)
  expect_equal(cc$n_hyper, 1L)
  expect_true(is.finite(cc$total))
})

test_that("the correction is the trace vcov()'s two matrices already imply", {
  skip_on_cran()
  # the identity test-inference.R asserts on a smooth --
  #   tr((V' - V_b) H) == statmod_edf_correction()$total
  # -- extended to a model carrying a filter. It was BROKEN there before this
  # change and not by a little: hyper_correction(), which
  # vcov(type = "unconditional") reads, has passed joint = TRUE since
  # 0.114.0 while the correction passed joint = FALSE, so one consumer of
  # hyper_mode_cross() counted the movement and the other counted zero.
  d <- edf_panel()
  fit <- statmod(y ~ gas(p = 1, q = 1, by = g, alpha1 ~ 1 + random(~ 1 | g)),
                 gaussian1_distrib(), d, outer_criterion = reml())
  dz <- statmod_design(fit@spec)
  K <- as.matrix(statmod_full_information(fit@spec, fit@coefficients, dz))
  C <- as.matrix(vcov(fit, type = "unconditional", readable = FALSE)) -
    as.matrix(vcov(fit, readable = FALSE))
  expect_equal(dim(C), dim(K))

  # ORDERING. vcov() orders by EQUATION and the joint vector does not: it
  # gives mu's coefficients, then the term's own parameters, then sigma's,
  # where statmod_marginal_full() gives the stacked coefficients first and
  # the term's parameters after. The two are permutations of each other, so a
  # trace taken across them is not a trace -- measured, comparing them
  # unpermuted reads 2.3043 against 1.1089, which is the ordering and not a
  # disagreement between the two routes.
  np <- vapply(dz[fit@spec@distrib@params], function(b) b$npar, integer(1))
  nb <- sum(np)
  nz <- nrow(K) - nb
  n1 <- np[[1L]]
  perm <- c(seq_len(n1), n1 + nz + seq_len(nb - n1), n1 + seq_len(nz))
  expect_setequal(perm, seq_len(nrow(K)))

  cc <- statmod_edf_correction(fit@spec, fit@coefficients, fit@hyper, dz,
                               fit@methods$outer)
  expect_equal(sum(C[perm, perm] * t(K)), cc$total, tolerance = 1e-8)
})
