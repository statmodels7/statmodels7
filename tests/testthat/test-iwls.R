# The objective and the scoring step.
#
# The gaussian with an identity link and no penalty is the case where a
# tolerance need not be chosen: the model IS least squares, so the
# coefficients must equal lm()'s and the scale the maximum likelihood one.

set.seed(1)
dd <- data.frame(x = runif(80), z = runif(80),
                 g = factor(rep(letters[1:4], 20)))
dd$y <- 1 + 2 * dd$x + rnorm(80, sd = 0.7)

fit_pieces <- function(spec, method = iwls()) {
  design <- statmod_design(spec)
  hyper <- statmod_hyper_start(spec)
  obj <- statmod_objective(spec, hyper, design)
  start <- rep(0, sum(obj$npar))
  start[1L] <- mean(spec@response)
  start[obj$npar[1L] + 1L] <- log(stats::sd(spec@response))
  res <- iwls_fit(obj, start, method, spec@n_obs,
                  function(v) iwls_pieces(spec, design, obj$split(v), hyper,
                                          method))
  list(res = res, coef = obj$split(res$par), obj = obj)
}

test_that("a gaussian fit is least squares", {
  spec <- statmod_spec(y ~ x + g, distributions7::gaussian1_distrib(), dd)
  out <- fit_pieces(spec)
  ref <- stats::lm(y ~ x + g, dd)
  expect_true(out$res$converged)
  expect_equal(out$coef$mu, unname(stats::coef(ref)), tolerance = 1e-10)
  expect_equal(exp(out$coef$sigma),
               sqrt(mean(stats::residuals(ref)^2)), tolerance = 1e-10)
  expect_equal(-out$res$value, as.numeric(stats::logLik(ref)),
               tolerance = 1e-10)
})

test_that("every decomposition reaches the same point", {
  spec <- statmod_spec(y ~ x + g, distributions7::gaussian1_distrib(), dd)
  ref <- unname(stats::coef(stats::lm(y ~ x + g, dd)))
  for (how in c("qr", "svd", "chol", "chol_crossprod")) {
    out <- fit_pieces(spec, iwls(decomposition = how))
    expect_equal(out$coef$mu, ref, tolerance = 1e-9, info = how)
  }
})

test_that("the augmented route is taken, and says when it is not", {
  spec <- statmod_spec(y ~ x, distributions7::gaussian1_distrib(), dd)
  out <- fit_pieces(spec, iwls(decomposition = "qr"))
  expect_true(all(out$res$history$route == "qr"))
  # the observed curvature away from the optimum is routinely indefinite, so
  # the square root does not exist and the step falls back, reporting it
  out2 <- fit_pieces(spec, iwls(hessian = "observed", decomposition = "qr"))
  expect_true(out2$res$converged)
  expect_true(any(grepl("chol", out2$res$history$route)))
})

test_that("the square-root design squares to the information", {
  # R'R must BE the information the assembled route computes, or the two
  # decompositions are solving different problems
  spec <- statmod_spec(y ~ x | sigma ~ z, distributions7::gaussian1_distrib(),
                       dd)
  design <- statmod_design(spec)
  coef <- list(mu = c(0.5, 1.5), sigma = c(-0.2, 0.3))
  th <- statmodels7:::statmod_eta(spec, design, coef)$theta
  L <- statmodels7:::chol_blocks(statmodels7:::info_blocks(spec, th))
  R <- statmodels7:::sqrt_design(design, L)
  I <- statmodels7:::statmod_information_at(spec, coef, design)
  expect_equal(crossprod(R), I, tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("the penalty's factor squares to its Hessian", {
  S <- crossprod(diff(diag(6), differences = 2))
  C <- statmodels7:::penalty_sqrt(S)
  expect_equal(crossprod(C), S, tolerance = 1e-10, ignore_attr = TRUE)
  # a rank-deficient penalty keeps only the rows its range needs
  expect_identical(nrow(C), 4L)
  # an indefinite one has no square root and says so
  expect_null(statmodels7:::penalty_sqrt(diag(c(1, -1))))
})

test_that("the score and the information agree with numerical differentiation", {
  skip_if_not_installed("numDeriv")
  spec <- statmod_spec(y ~ x | sigma ~ z, distributions7::gaussian1_distrib(),
                       dd)
  design <- statmod_design(spec)
  hyper <- statmod_hyper_start(spec)
  obj <- statmod_objective(spec, hyper, design)
  v <- c(0.8, 1.7, -0.4, 0.2)
  expect_equal(obj$gr(v), numDeriv::grad(obj$fn, v), tolerance = 1e-7)
  # the observed information is minus the Hessian of the log-likelihood, so
  # against the objective's own Hessian it is the objective's curvature
  I <- statmodels7:::statmod_information_at(spec, obj$split(v), design,
                                            expected = FALSE)
  expect_equal(I, numDeriv::hessian(obj$fn, v), tolerance = 1e-5,
               ignore_attr = TRUE)
})

test_that("a penalized term shrinks towards zero", {
  dd2 <- dd
  dd2$R <- matrix(rnorm(160), 80, 2)
  spec <- statmod_spec(y ~ x + ridge(R), distributions7::gaussian1_distrib(),
                       dd2)
  design <- statmod_design(spec)
  hyper <- statmod_hyper_start(spec)
  obj <- statmod_objective(spec, hyper, design)
  start <- rep(0, sum(obj$npar))
  start[1L] <- mean(dd2$y)
  start[obj$npar[1L] + 1L] <- log(stats::sd(dd2$y))
  m <- iwls()
  res <- iwls_fit(obj, start, m, spec@n_obs,
                  function(v) iwls_pieces(spec, design, obj$split(v), hyper, m))
  expect_true(res$converged)
  # the penalized block is smaller than the unpenalized fit of the same block
  cols <- design$mu$blocks[[setdiff(names(design$mu$blocks), "linpar")]]
  pen_coef <- obj$split(res$par)$mu[cols]
  unpen <- stats::lm(y ~ x + R, dd2)
  expect_lt(sum(pen_coef^2), sum(unname(stats::coef(unpen))[3:4]^2))
})

test_that("iwls states what it rejects and prints itself", {
  expect_error(iwls(decomposition = "lu"), "arg")
  expect_error(iwls(maxit = 0), "positive integer")
  expect_error(iwls(tol = -1), "positive number")
  expect_output(print(iwls()), "expected information, qr")
})
