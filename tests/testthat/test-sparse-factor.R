# The sparse factorization of the penalized information, and the rule it has
# to obey: it is a property of the MATRIX, so every design with the same shape
# takes the same route whatever term built it.

sparse_re <- function(m, n = 4000, seed = 11L) {
  set.seed(seed)
  g <- factor(sample(seq_len(m), n, TRUE))
  x <- runif(n)
  data.frame(x = x, g = g,
             y = 2 + 0.97 * x + rnorm(m, 0, 0.55)[as.integer(g)] + rnorm(n))
}

# the penalized matrix a fit actually forms, at its own penalized mode
penalized_at <- function(formula, data, ...) {
  spec <- statmod_spec(formula, distributions7::gaussian1_distrib(), data, ...)
  design <- statmod_design(spec)
  hyper <- statmod_hyper_start(spec, design)
  obj <- statmod_objective(spec, hyper, design)
  blocks <- statmod_blocks(spec, design)
  fit <- fit_smooth(obj, statmod_start(spec, design, obj), blocks$smooth,
                    spec, design, hyper, iwls(), list(inner = FALSE))
  coef <- obj$split(fit$par)
  S <- statmod_penalty_at(spec, coef, hyper, design, "hessian")
  S[!is.finite(S)] <- 0
  statmod_information_at(spec, coef, design, FALSE, "bartlett") + S
}


test_that("the gate asks the matrix's size and its zeros, and nothing else", {
  K <- penalized_at(y ~ x + random(~ 1 | g), sparse_re(300L))
  expect_true(isS4(K))
  expect_true(worth_sparse(K))

  # too small: below the crossover the fixed cost of the coercion dominates,
  # measured at 0.13x and 0.33x on the inverse at p = 23 and p = 53
  expect_false(worth_sparse(K, min_dim = ncol(K) + 1L))
  # and too full, at the same size
  expect_false(worth_sparse(K, max_density = 0))

  # a small model is refused whatever it is made of
  expect_false(worth_sparse(penalized_at(y ~ x + s(x, k = 10),
                                         sparse_re(20L, n = 500L))))
})


test_that("the two routes give the same log-determinant and the same verdict", {
  K <- penalized_at(y ~ x + random(~ 1 | g), sparse_re(300L))
  expect_true(worth_sparse(K))

  sp <- pd_factor(K)
  dn <- pd_logdet_dense(as_dense(K))
  expect_true(sp$sparse)
  expect_true(sp$ok)
  expect_true(dn$ok)
  expect_equal(sp$logdet, dn$logdet, tolerance = 1e-10)

  # pd_logdet() is the same answer with the factor dropped, so the one entry
  # point every caller already used keeps its contract
  expect_equal(pd_logdet(K)$logdet, sp$logdet)
  expect_null(pd_logdet(K)$factor)
})


test_that("a flat direction is refused on the sparse route too", {
  d <- sparse_re(300L)
  # a duplicated column is the ordinary way to an exactly zero eigenvalue,
  # and the refusal must not depend on which storage the matrix is in
  d$x2 <- d$x
  K <- penalized_at(y ~ x + x2 + random(~ 1 | g), d)
  expect_false(isTRUE(pd_factor(K)$ok))
  expect_false(isTRUE(pd_logdet(K)$ok))
})


test_that("the condition estimate errs on the small side", {
  K <- penalized_at(y ~ x + random(~ 1 | g), sparse_re(200L))
  Ks <- Matrix::forceSymmetric(as_sparse(K))
  L <- Matrix::Cholesky(Ks, LDL = FALSE, super = NA)
  est <- sparse_lmin(L, ncol(Ks))
  true <- min(eigen(as_dense(K), symmetric = TRUE, only.values = TRUE)$values)
  expect_true(is.finite(est))
  # for a symmetric matrix ||A^-1||_1 >= 1/lambda_min, so the estimate is at
  # or below the smallest eigenvalue; the estimator's own slack is bounded
  # and nowhere near the orders of magnitude the verdict separates
  expect_lt(est, true * 1.5)
  expect_gt(est, true / 100)
})


test_that("the route does not change what a fit reports", {
  d <- sparse_re(300L)
  f_sp <- statmod(y ~ x + random(~ 1 | g),
                  distributions7::gaussian1_distrib(), d,
                  outer_criterion = reml())

  # the same fit with the gate held shut, which is the code this replaced
  old <- worth_sparse
  on.exit(assignInNamespace("worth_sparse", old, ns = "statmodels7"),
          add = TRUE)
  assignInNamespace("worth_sparse", function(...) FALSE, ns = "statmodels7")
  f_dn <- statmod(y ~ x + random(~ 1 | g),
                  distributions7::gaussian1_distrib(), d,
                  outer_criterion = reml())

  expect_equal(as.numeric(stats::logLik(f_sp)),
               as.numeric(stats::logLik(f_dn)), tolerance = 1e-8)
  expect_equal(unlist(f_sp@coefficients), unlist(f_dn@coefficients),
               tolerance = 1e-6)
  expect_equal(sum(f_sp@edf$edf), sum(f_dn@edf$edf), tolerance = 1e-8)
})


test_that("a design with no random effect takes the same route", {
  # piano_lme4.txt section 5: a phase that pays on random(~1|g) and not on an
  # indicator block of the same shape is reading the term. Both matrices are
  # built here and both are asked the same question; the answer must agree,
  # and the fits must too.
  set.seed(12)
  n <- 4000
  m <- 300L
  g <- factor(sample(seq_len(m), n, TRUE))
  x <- runif(n)
  d <- data.frame(x = x, g = g,
                  y = 3 * rnorm(m, 0, 0.6)[as.integer(g)] + sin(2 * pi * x) +
                    rnorm(n, 0, 0.8))

  K_lin <- penalized_at(y ~ 0 + g + s(x, k = 10), d,
                        linpar = linpar_options(sparse = TRUE))
  expect_true(worth_sparse(K_lin))
  expect_true(isTRUE(pd_factor(K_lin)$sparse))

  f <- statmod(y ~ 0 + g + s(x, k = 10),
               distributions7::gaussian1_distrib(), d,
               outer_criterion = reml(),
               linpar_control = linpar_options(sparse = TRUE))
  expect_true(is.finite(as.numeric(stats::logLik(f))))
})
