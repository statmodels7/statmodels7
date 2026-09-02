# What a fit reports about its own uncertainty.

set.seed(4)
n <- 400
dd <- data.frame(x = runif(n), z = runif(n),
                 g = factor(rep(c("a", "b"), length.out = n)))

test_that("an unpenalized gaussian fit reproduces the closed-form variance", {
  # the reference is the linear model's own, which shares no code with the
  # penalized information assembled here
  sim <- rstatmod(y ~ x + g, distributions7::gaussian1_distrib(), dd,
                  par = list(mu = c(1, 2, -0.5), sigma = log(0.4)))$data
  fit <- statmod(y ~ x + g, distributions7::gaussian1_distrib(), sim)
  V <- vcov(fit)

  lm_fit <- stats::lm(y ~ x + g, sim)
  s2 <- sum(stats::residuals(lm_fit)^2) / n   # the MLE, not the unbiased one
  V_lm <- s2 * solve(crossprod(stats::model.matrix(lm_fit)))

  idx <- grep("^mu:", rownames(V))
  expect_equal(unname(V[idx, idx]), unname(V_lm), tolerance = 1e-6)

  # and the scale block: eta = log(sigma), so the information is 2n and the
  # standard error 1/sqrt(2n) whatever sigma is
  expect_equal(unname(sqrt(V["sigma:(Intercept)", "sigma:(Intercept)"])),
               1 / sqrt(2 * n), tolerance = 1e-6)
})

test_that("with no penalty the two conventions are the same matrix", {
  sim <- rstatmod(y ~ x, distributions7::gaussian1_distrib(), dd,
                  par = list(mu = c(1, 2), sigma = log(0.4)))$data
  fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), sim)
  expect_equal(vcov(fit, "bayesian"), vcov(fit, "frequentist"),
               tolerance = 1e-10)
})

test_that("a ridge makes the bayesian variance the larger of the two", {
  sim <- rstatmod(y ~ ridge(~ x + z), distributions7::gaussian1_distrib(), dd,
                  par = list(mu = c(1, 2, -1), sigma = log(0.5)))$data
  fit <- statmod(y ~ ridge(~ x + z), distributions7::gaussian1_distrib(), sim)
  vb <- diag(vcov(fit, "bayesian"))
  vf <- diag(vcov(fit, "frequentist"))
  expect_true(all(vb >= vf - 1e-12))
  expect_true(any(vb > vf * 1.000001))
  # V_f = V_b H V_b is the sandwich the definition states, computed here from
  # the pieces rather than read back from the function under test
  H <- -hessian(fit, expected = TRUE)
  Vb <- vcov(fit, "bayesian")
  expect_equal(Vb %*% H %*% Vb, vcov(fit, "frequentist"), tolerance = 1e-8)
})

test_that("vcov agrees with a numerical Hessian of the objective", {
  # the assembled information is a sum of crossprods; numDeriv on the whole
  # penalized objective shares none of that arithmetic
  skip_if_not_installed("numDeriv")
  sim <- rstatmod(y ~ x | sigma ~ z, distributions7::gaussian1_distrib(), dd,
                  par = list(mu = c(1, 2), sigma = c(-1, 0.8)))$data
  fit <- statmod(y ~ x | sigma ~ z, distributions7::gaussian1_distrib(), sim)
  spec <- fit@spec
  design <- statmod_design(spec)
  obj <- statmod_objective(spec, fit@hyper, design, expected = FALSE)
  b <- unlist(fit@coefficients[spec@distrib@params], use.names = FALSE)
  Hn <- numDeriv::hessian(obj$fn, b)
  expect_equal(unname(vcov(fit, expected = FALSE)), solve(Hn),
               tolerance = 1e-5)
})

test_that("a confidence interval is symmetric and covers the truth", {
  sim <- rstatmod(y ~ x, distributions7::gaussian1_distrib(), dd,
                  par = list(mu = c(1, 2), sigma = log(0.4)))$data
  fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), sim)
  ci <- confint(fit)
  expect_named(ci, c("parameter", "term", "coefficient", "estimate", "se",
                     "lower", "upper"))
  # a coefficient of a linear predictor is free on the whole line, so the
  # interval is symmetric about the estimate on the scale it is built on
  expect_equal(ci$upper - ci$estimate, ci$estimate - ci$lower,
               tolerance = 1e-12)
  # 400 observations: the truth is inside
  expect_true(ci["mu:x", "lower"] < 2 && 2 < ci["mu:x", "upper"])
  expect_true(ci["mu:(Intercept)", "lower"] < 1 &&
                1 < ci["mu:(Intercept)", "upper"])
  # the width follows the level
  wide <- confint(fit, level = 0.99)
  expect_true(all(wide$upper - wide$lower > ci$upper - ci$lower))
})

test_that("confint selects by parameter, by term and by label", {
  sim <- rstatmod(y ~ x | sigma ~ z, distributions7::gaussian1_distrib(), dd,
                  par = list(mu = c(1, 2), sigma = c(-1, 0.8)))$data
  fit <- statmod(y ~ x | sigma ~ z, distributions7::gaussian1_distrib(), sim)
  expect_identical(unique(confint(fit, "sigma")$parameter), "sigma")
  expect_identical(nrow(confint(fit, "mu:x")), 1L)
  expect_error(confint(fit, "nonesuch"), "matched nothing")
  expect_error(confint(fit, level = 1), "strictly between")
})

test_that("a coefficient a lasso set to zero has no variance", {
  # at the kink the penalty is not twice differentiable, so there is no
  # curvature to read and NA is the reading rather than a number
  sim <- dd
  sim$y <- 1 + 2 * dd$x + stats::rnorm(n, sd = 0.4)
  sim$noise1 <- stats::runif(n)
  sim$noise2 <- stats::runif(n)
  fit <- statmod(y ~ x + lasso(~ noise1 + noise2, lambda = 200),
                 distributions7::gaussian1_distrib(), sim)
  b <- fit@coefficients$mu
  V <- vcov(fit)
  zero <- which(b == 0)
  expect_gt(length(zero), 0)
  expect_true(all(is.na(diag(V)[zero])))
  # while the unpenalized ones are perfectly available
  expect_true(all(is.finite(diag(V)[grep("^mu:x$", rownames(V))])))
  s <- summary(fit)
  expect_true(any(grepl("kinked penalty", s@notes)))
})

test_that("a summary carries the tables, the criteria and the notes", {
  sim <- rstatmod(y ~ x | sigma ~ z, distributions7::gaussian1_distrib(), dd,
                  par = list(mu = c(1, 2), sigma = c(-1, 0.8)))$data
  fit <- statmod(y ~ x | sigma ~ z, distributions7::gaussian1_distrib(), sim)
  s <- summary(fit)
  expect_named(s@tables, c("mu", "sigma"))
  # one block, the parametric one, with a row per coefficient
  expect_identical(length(s@tables$mu), 1L)
  expect_identical(s@tables$mu[[1L]]$kind, "parametric")
  expect_identical(nrow(s@tables$mu[[1L]]$table), 2L)
  expect_true(all(c("estimate", "se", "statistic", "p_value", "lower",
                    "upper", "role") %in% names(s@tables$mu[[1L]]$table)))
  # the criteria are built on the same log-likelihood and df logLik reports,
  # so AIC() and the summary cannot disagree
  expect_equal(s@aic, stats::AIC(fit), tolerance = 1e-12)
  expect_equal(s@bic, stats::BIC(fit), tolerance = 1e-12)
  expect_equal(s@loglik, as.numeric(stats::logLik(fit)), tolerance = 1e-12)
  # a slope of 2 on 400 observations is not in doubt
  expect_lt(s@tables$mu[[1L]]$table$p_value[2L], 1e-10)
  expect_output(print(s), "log-likelihood")
  expect_output(print(s), "Distribution")
})

test_that("a penalized term is a block of its own, with its hyperparameter", {
  sim <- rstatmod(y ~ ridge(~ x + z), distributions7::gaussian1_distrib(), dd,
                  par = list(mu = c(1, 2, -1), sigma = log(0.5)))$data
  fit <- statmod(y ~ ridge(~ x + z), distributions7::gaussian1_distrib(),
                 sim, outer_criterion = NULL)
  s <- summary(fit)
  kinds <- vapply(s@tables$mu, `[[`, character(1), "kind")
  expect_true("penalized" %in% kinds)
  b <- s@tables$mu[[which(kinds == "penalized")]]
  # a ridge's coefficients stay interpretable, so they are shown, and the
  # hyperparameter is a row of its own marked as held rather than estimated
  expect_true(any(b$table$role == "coefficient"))
  expect_true(any(b$table$role == "fixed"))
  expect_true(all(is.na(b$table$se[b$table$role == "fixed"])))
  expect_true(any(grepl("held at the value", s@notes)))
})

test_that("a smooth shows its linear part and its edf, not its basis", {
  set.seed(12)
  n2 <- 400
  ds <- data.frame(x = runif(n2, -2, 2))
  ds$y <- sin(1.4 * ds$x) + stats::rnorm(n2, sd = 0.3)
  fit <- statmod(y ~ s(x, k = 10, lambda = 5),
                 distributions7::gaussian1_distrib(), ds)
  s <- summary(fit)
  kinds <- vapply(s@tables$mu, `[[`, character(1), "kind")
  b <- s@tables$mu[[which(kinds == "smooth")]]

  # a block of several columns, and exactly two rows: the linear component and
  # the smoothing parameter. The coordinates of the orthonormal basis say
  # nothing one at a time and what they say together is the edf.
  expect_gt(b$n_coef, 4)
  expect_identical(nrow(b$table), 2L)
  expect_lt(nrow(b$table), b$n_coef)
  expect_identical(b$table$name, c("lambda", "s(x).lin"))
  expect_identical(b$table$role, c("fixed", "coefficient"))
  expect_equal(b$table$estimate[1L], 5)
  expect_true(is.finite(b$edf) && b$edf > 1 && b$edf < b$n_coef)
  # the linear component IS an ordinary coefficient and carries inference
  expect_true(is.finite(b$table$se[2L]))

  expect_output(print(s), "s(x, k = 10, lambda = 5)", fixed = TRUE)
  expect_output(print(s), "edf")
  expect_output(print(s), "[fixed]", fixed = TRUE)
})

test_that("a random effect shows its variance parameters, not its levels", {
  set.seed(13)
  n2 <- 400
  ds <- data.frame(g = factor(sample(paste0("g", 1:15), n2, replace = TRUE)))
  u <- stats::rnorm(15, sd = 0.7)
  ds$y <- 1 + u[as.integer(ds$g)] + stats::rnorm(n2, sd = 0.4)
  fit <- statmod(y ~ random(~ 1 | g), distributions7::gaussian1_distrib(),
                 ds, outer_criterion = NULL)
  s <- summary(fit)
  kinds <- vapply(s@tables$mu, `[[`, character(1), "kind")
  b <- s@tables$mu[[which(kinds == "random")]]

  expect_identical(b$n_coef, 15L)
  # fifteen effects, and not one of them is a row: what a reader wants is the
  # distribution they came from
  expect_true(all(b$table$role == "fixed"))
  expect_gt(nrow(b$table), 0)
  expect_false(any(startsWith(b$table$name, "random.")))
  expect_output(print(s), "random(~1 | g)", fixed = TRUE)
})

test_that("a selection lists the survivors and counts the zeros", {
  set.seed(14)
  n2 <- 300
  ds <- data.frame(x = runif(n2))
  for (j in 1:5) ds[[paste0("n", j)]] <- stats::rnorm(n2)
  ds$y <- 1 + 2 * ds$x + stats::rnorm(n2, sd = 0.4)
  fit <- statmod(y ~ x + lasso(~ n1 + n2 + n3 + n4 + n5, lambda = 200),
                 distributions7::gaussian1_distrib(), ds)
  s <- summary(fit)
  kinds <- vapply(s@tables$mu, `[[`, character(1), "kind")
  b <- s@tables$mu[[which(kinds == "selection")]]

  expect_identical(b$n_coef, 5L)
  expect_gt(b$n_zero, 0)
  shown <- b$table[b$table$role == "coefficient", , drop = FALSE]
  expect_identical(nrow(shown), b$n_coef - b$n_zero)
  expect_true(all(shown$estimate != 0))
  expect_output(print(s), "at zero")
})

test_that("a fit whose information is singular says so rather than guessing", {
  # two identical columns: the data does not identify their difference, and a
  # pseudo-inverse would report a standard error for a direction it cannot see
  sim <- rstatmod(y ~ x, distributions7::gaussian1_distrib(), dd,
                  par = list(mu = c(1, 2), sigma = log(0.4)))$data
  sim$xcopy <- sim$x
  expect_error(
    vcov(statmod(y ~ x + xcopy, distributions7::gaussian1_distrib(), sim)),
    "not positive definite")
})


test_that("a correlated random effect reports variance components", {
  # its hyperparameters are the free values of a matrix parameter, and nobody
  # reads the logarithm of a Cholesky diagonal: what the prior is ABOUT is the
  # standard deviations and the correlation of the effects
  skip_on_cran()
  set.seed(202)
  m <- 40L
  ni <- 12L
  g <- factor(rep(seq_len(m), each = ni))
  x <- stats::rnorm(m * ni)
  sd1 <- 1.2
  sd2 <- 0.6
  b1 <- stats::rnorm(m, 0, sd1)
  b2 <- stats::rnorm(m, 0, sd2)
  dd <- data.frame(
    y = 2 + 0.7 * x + b1[as.integer(g)] + b2[as.integer(g)] * x +
      stats::rnorm(m * ni, 0, 0.8),
    x = x, g = g)

  fit <- statmod(y ~ x + random(~ x | g),
                 distributions7::gaussian1_distrib(), dd,
                 outer_criterion = reml())
  s <- summary(fit)
  tb <- s@tables$mu[[2L]]$table
  expect_identical(tb$name[1:3], c("sd_v1", "sd_v2", "cor_v1_v2"))

  # the quantities, against the covariance the coordinates imply, computed
  # apart from the summary
  eta <- unlist(fit@hyper$mu[[1L]])
  sig <- parameters7::param_value(
    parameters7::log_cholesky(2), eta)
  expect_equal(tb$estimate[1L], sqrt(sig[1, 1]), tolerance = 1e-8)
  expect_equal(tb$estimate[2L], sqrt(sig[2, 2]), tolerance = 1e-8)
  expect_equal(tb$estimate[3L], sig[1, 2] / sqrt(sig[1, 1] * sig[2, 2]),
               tolerance = 1e-8)
  # and near the truth they were drawn from
  expect_gt(tb$estimate[1L], 0.8)
  expect_lt(tb$estimate[1L], 1.7)
  expect_gt(tb$estimate[2L], 0.3)
  expect_lt(tb$estimate[2L], 1.0)

  # each interval is built where the quantity ranges and mapped back, so a
  # standard deviation cannot be given a negative lower end and a correlation
  # cannot be given one that leaves (-1, 1)
  expect_true(all(is.finite(tb$se[1:3])))
  expect_gt(tb$lower[1L], 0)
  expect_gt(tb$lower[2L], 0)
  expect_gt(tb$lower[3L], -1)
  expect_lt(tb$upper[3L], 1)
  expect_true(all(tb$lower[1:3] < tb$estimate[1:3]))
  expect_true(all(tb$upper[1:3] > tb$estimate[1:3]))
  # no test is printed: the null a z would report on is that a standard
  # deviation is zero, which is the edge of its range
  expect_true(all(is.na(tb$statistic[1:3])))
})


test_that("a hyperparameter that IS its quantity is reported as it stands", {
  # the one-column case has a standard deviation for a hyperparameter, so
  # there is nothing to translate and the row keeps its own name
  skip_on_cran()
  set.seed(7)
  g <- factor(rep(seq_len(20), each = 10))
  x <- stats::rnorm(200)
  u <- stats::rnorm(20, 0, 0.9)
  dd <- data.frame(y = 1 + 2 * x + u[as.integer(g)] + stats::rnorm(200, 0, 0.5),
                   x = x, g = g)
  fit <- statmod(y ~ x + random(~ 1 | g),
                 distributions7::gaussian1_distrib(), dd,
                 outer_criterion = reml())
  tb <- summary(fit)@tables$mu[[2L]]$table
  expect_identical(tb$name[1L], "sigma")
  expect_equal(tb$estimate[1L], fit@hyper$mu[[1L]][["sigma"]])
})


test_that("a hyperparameter the readable block does not describe keeps its row", {
  # a multivariate Student t is ABOUT the standard deviations and the
  # correlations of its scale matrix, and its degrees of freedom are none of
  # those: replacing the coordinate rows wholesale dropped nu from the summary
  skip_on_cran()
  set.seed(78)
  m <- 40L
  ni <- 12L
  g <- factor(rep(seq_len(m), each = ni))
  x <- stats::rnorm(m * ni)
  b <- matrix(stats::rnorm(2 * m, 0, 0.9), ncol = 2L)
  dd <- data.frame(
    y = 1 + 0.8 * x + b[as.integer(g), 1L] + b[as.integer(g), 2L] * x +
      stats::rnorm(m * ni, 0, 0.7), x = x, g = g)
  mvt <- do.call(distributions7::fixed,
                 list(distributions7::mvstudent_t1_distrib(2),
                      mu1 = 0, mu2 = 0))
  fit <- statmod(y ~ x + random(~ x | g, distrib = mvt),
                 distributions7::gaussian1_distrib(), dd,
                 outer_criterion = reml())
  tb <- summary(fit)@tables$mu[[2L]]$table
  expect_identical(tb$name[1:4],
                   c("scale_sd_v1", "scale_sd_v2", "cor_v1_v2", "nu"))
  expect_equal(tb$estimate[4L], fit@hyper$mu[[1L]][["nu"]])
})


test_that("a shape parameter is reported with a standard error and interval", {
  # the question is whether EVERY estimated hyperparameter carries one, not
  # only a scale: this prior is a log-transformed gamma, whose free parameter
  # is the variance of the gamma underneath, and the effects are drawn from it
  skip_on_cran()
  set.seed(80)
  m <- 60L
  ni <- 12L
  a <- 4
  g <- factor(rep(seq_len(m), each = ni))
  x <- stats::rnorm(m * ni)
  u <- log(stats::rgamma(m, shape = a, rate = a))
  dd <- data.frame(y = 1.5 + 0.8 * x + u[as.integer(g)] +
                     stats::rnorm(m * ni, 0, 0.7), x = x, g = g)
  # the mean is held at ONE, which is what centers the logarithm: the prior's
  # own mean is digamma(a) - log(a), within sigma2/2 of zero
  lg <- distributions7::fixed(
    distributions7::transformation(distributions7::gamma2_distrib(),
                                   distributions7::log_transform()), mu = 1)
  fit <- statmod(y ~ x + random(~ 1 | g, distrib = lg),
                 distributions7::gaussian1_distrib(), dd,
                 outer_criterion = reml())
  tb <- summary(fit)@tables$mu[[2L]]$table
  expect_identical(tb$name[1L], "sigma2")
  expect_true(is.finite(tb$se[1L]))
  expect_lt(tb$lower[1L], tb$estimate[1L])
  expect_gt(tb$upper[1L], tb$estimate[1L])
  expect_gt(tb$lower[1L], 0)
  # and near the 1/a it was drawn from
  expect_gt(tb$estimate[1L], 0.1)
  expect_lt(tb$estimate[1L], 0.5)
})


test_that("the certificate is a property of the point, not of the search", {
  # WHY IT EXISTS. The convergence flag says whether a search stopped on its
  # own rule, which is a statement about the search, and measured across
  # shapes it does not order fits by quality: on one model the default
  # reported success at a criterion of -1783.47 while the same data under
  # lbfgs() reached -1664.43 and reported failure. The certificate is read at
  # the reported point instead.
  set.seed(41)
  n <- 300
  d <- data.frame(x = runif(n))
  d$y <- sin(5 * d$x) + rnorm(n, 0, 0.3)

  fit <- statmod(y ~ s(x, k = 10), gaussian1_distrib(), d,
                 outer_criterion = reml())
  ct <- statmod_certificate(fit)
  expect_identical(ct$state, "converged")
  expect_lt(ct$gradient, 1e-2)
  expect_true(is.finite(ct$mode_error))
  expect_length(ct$boundary, 0L)
  expect_length(ct$reason, 0L)

  # IT MUST BE ABLE TO REFUSE, or "converged" says nothing. The same fit
  # against a tolerance below the gradient it actually carries.
  strict <- statmod_certificate(fit, tol = ct$gradient / 10)
  expect_identical(strict$state, "not converged")
  expect_match(strict$reason[1], "gradient")
})


test_that("the certificate declares a boundary and refuses where it cannot read", {
  # A SMOOTH ON PURE NOISE is the boundary case that is a right answer: REML
  # penalizes the wiggle away and the smoothing parameter runs to the edge.
  # The certificate says so by name rather than reporting an ordinary fit.
  set.seed(42)
  n <- 300
  d <- data.frame(x = runif(n), y = rnorm(300))
  noise <- statmod(y ~ s(x, k = 10), gaussian1_distrib(), d,
                   outer_criterion = reml())
  cn <- statmod_certificate(noise)
  expect_identical(cn$state, "boundary")
  expect_gt(length(cn$boundary), 0L)

  # AND WHERE THERE IS NO OUTER GRADIENT THE GRADIENT IS NA, always: 2p refits
  # to difference one would cost more than the fit. What the STATE is there
  # depends on why, and that is the next two tests' subject -- a model with no
  # penalty is certified by its own mode, a kinked one is refused.
  plain <- statmod(y ~ x, gaussian1_distrib(), d, outer_criterion = NULL)
  cp <- statmod_certificate(plain)
  expect_true(is.na(cp$gradient))
  expect_true(is.finite(cp$mode_error))
  expect_length(cp$boundary, 0L)
})


test_that("summary carries the certificate and prints it", {
  set.seed(43)
  n <- 200
  d <- data.frame(x = runif(n))
  d$y <- sin(4 * d$x) + rnorm(n, 0, 0.3)
  s <- summary(statmod(y ~ s(x, k = 8), gaussian1_distrib(), d,
                       outer_criterion = reml()))
  expect_false(is.null(s@certificate))
  expect_true(s@certificate$state %in%
                c("converged", "boundary", "not converged", "unknown"))
  out <- utils::capture.output(print(s))
  expect_true(any(grepl("certificate:", out, fixed = TRUE)))
})


test_that("a model with no penalty is certified by its own mode", {
  # THE LARGEST GAP THE BATTERY FOUND: seven of twenty-nine cases had no outer
  # gradient and so no state at all. Two quite different reasons hide there,
  # and only one of them is a real refusal.
  #
  # NO PENALTY AT ALL -- linpar, nl, seg, jump, jseg -- means there is no
  # hyperparameter for a gradient to be about, and the only question left is
  # whether the inner fit reached its mode. The mode error answers it.
  set.seed(51)
  n <- 200
  d <- data.frame(x = runif(n))
  d$y <- 1 + 2 * d$x + rnorm(n, 0, 0.3)
  fit <- statmod(y ~ x, gaussian1_distrib(), d, outer_criterion = NULL)
  ct <- statmod_certificate(fit)
  expect_identical(ct$state, "converged")
  expect_true(is.na(ct$gradient))
  expect_lt(ct$mode_error, mode_error_limit())
  expect_match(ct$reason[1], "no penalty")

  # AND IT MUST BE ABLE TO REFUSE HERE TOO. The limit is a property of the
  # package, so the injection moves the reading rather than the limit: a fit
  # whose mode error exceeds it is not certified. Measured, a `jump` fitted to
  # data carrying a slope and a slope change it has no term for reads 1.215,
  # against 4.2e-04 on data where the jump is the truth.
  local_mocked_bindings(mode_error_limit = function() ct$mode_error / 10)
  expect_identical(statmod_certificate(fit)$state, "not converged")
})


test_that("a kinked hyperparameter is refused rather than read", {
  # A lasso's lambda is chosen along a PATH, being the argmin over a grid
  # rather than the root of a derivative, so there is no outer gradient. And
  # the mode error is not a reading either: at a coefficient the penalty has
  # set to zero the score does not vanish but lies in the subdifferential --
  # measured on a lasso, a mode error of 4.7e-03 carried by a coordinate whose
  # coefficient is exactly 0 and whose score is -0.715.
  set.seed(52)
  n <- 200
  X <- matrix(runif(n * 4), n, 4)
  d <- data.frame(x1 = X[, 1], x2 = X[, 2], x3 = X[, 3], x4 = X[, 4])
  d$y <- 2 * d$x1 + rnorm(n, 0, 0.3)
  fit <- statmod(y ~ lasso(~ x1 + x2 + x3 + x4), gaussian1_distrib(), d)
  ct <- statmod_certificate(fit)
  expect_identical(ct$state, "unknown")
  expect_match(ct$reason[1], "kink")
  # the distinction is the point: this reason is NOT the no-penalty one
  expect_false(grepl("no penalty", ct$reason[1], fixed = TRUE))
})


test_that("every hyperparameter carries the mark its note speaks of", {
  # The note at the foot says "a hyperparameter marked REML was estimated by
  # that criterion", and no such mark was ever printed: the source went into
  # the column where a standard error would have been, which is occupied for
  # a REML-estimated one and blank only for a held or path-chosen one. It
  # goes in the NAME now, so the three kinds are marked alike and the note
  # is true. Without this the defect is invisible -- the note reads as if it
  # described something on screen.
  sim <- rstatmod(y ~ x + z, distributions7::gaussian1_distrib(), dd,
                  par = list(mu = c(1, 2, -0.5), sigma = log(0.4)))$data
  s <- summary(statmod(y ~ s(x, k = 8),
                       distributions7::gaussian1_distrib(), sim))
  expect_output(print(s), "[reml]", fixed = TRUE)
  # and the note that speaks of it is emitted
  expect_true(any(grepl("marked REML", s@notes, fixed = TRUE)))
})

test_that("the note about the point follows the certificate, not the flag", {
  # WHETHER THE POINT IS A MAXIMUM is the certificate's question and the
  # search's flag is a different one, so the note that says the estimates are
  # read where the surface is still moving is emitted from the certificate.
  # It used to read the flag, and was wrong in both directions: printed under
  # a certificate saying CONVERGED, and absent on a fit whose search met its
  # rule 0.371 log-likelihood units above its mode.
  sim <- rstatmod(y ~ x + z, distributions7::gaussian1_distrib(), dd,
                  par = list(mu = c(1, 2, -0.5), sigma = log(0.4)))$data
  for (f in list(statmod(y ~ x + z, distributions7::gaussian1_distrib(), sim),
                 statmod(y ~ s(x, k = 8),
                         distributions7::gaussian1_distrib(), sim))) {
    s <- summary(f)
    warned <- any(grepl("not certified as a maximum", s@notes, fixed = TRUE))
    ok <- !is.null(s@certificate) &&
      s@certificate$state %in% c("converged", "boundary")
    expect_identical(warned, !ok)
  }
})

test_that("a developed parameter is reported as a compartment of its own", {
  set.seed(3)
  m <- 6
  ni <- 30
  id <- factor(rep(seq_len(m), each = ni))
  x <- runif(m * ni, 0, 10)
  psi_i <- 5 + rnorm(m, 0, 0.6)
  y <- 0.3 * x + 1.5 * pmax(x - psi_i[as.integer(id)], 0) +
    rnorm(m * ni, 0, 0.5)
  dd <- data.frame(y = y, x = x, id = id)
  fit <- statmod(y ~ seg(x, psi ~ random(~ 1 | id), psi = 5),
                 gaussian1_distrib(), dd)
  s <- summary(fit)
  b <- s@tables$mu[[which(vapply(s@tables$mu, `[[`, character(1), "kind") ==
                          "breakpoint")]]
  expect_identical(length(b$components), 1L)
  cp <- b$components[[1L]]
  expect_identical(cp$name, "psi1")
  # the term's OWN table keeps what is not developed, and the developed
  # parameter's columns are gone from it
  expect_true(all(c("beta", "gamma1") %in%
                    sub("^[^.]+\\.", "", b$table$name)))
  expect_identical(nrow(b$table), 2L)
  # the compartment carries its own hyperparameter and the population value,
  # and NOT the predictions: a column of per-group numbers is what coef() is
  # for
  expect_true(any(b$components[[1L]]$table$role == "estimated"))
  expect_true("(Intercept)" %in% cp$table$name)
  expect_false(any(startsWith(cp$table$name, "random.")))
  expect_lt(nrow(cp$table), cp$n_coef)
  # and reports how many there are instead
  expect_match(cp$lines[[1L]], "^[0-9]+ predictions, sd")

  out <- paste(utils::capture.output(print(s)), collapse = "\n")
  expect_match(out, "psi1  ~ intercept + random", fixed = TRUE)
  expect_match(out, "effect sd", fixed = TRUE)
  expect_match(out, "predictions, sd", fixed = TRUE)
  expect_false(grepl("random.1", out, fixed = TRUE))
})

test_that("a term that develops nothing keeps the flat table it had", {
  set.seed(4)
  dd <- data.frame(x = runif(200, 0, 10))
  dd$y <- 0.3 * dd$x + 1.5 * pmax(dd$x - 5, 0) + rnorm(200, 0, 0.5)
  fit <- statmod(y ~ seg(x, psi = 4), gaussian1_distrib(), dd)
  s <- summary(fit)
  b <- s@tables$mu[[which(vapply(s@tables$mu, `[[`, character(1), "kind") ==
                          "breakpoint")]]
  expect_identical(b$components, list())
  expect_null(b$head)
  # the readable quantity, as before: the position and not the pair of
  # working coefficients it is read off
  expect_true("psi1" %in% b$table$name)
})

test_that("the head shows every parameter of a term at once", {
  set.seed(5)
  m <- 5
  ni <- 30
  id <- factor(rep(seq_len(m), each = ni))
  x <- runif(m * ni, 0, 10)
  psi_i <- 5 + rnorm(m, 0, 0.5)
  y <- 0.3 * x + 1.5 * pmax(x - psi_i[as.integer(id)], 0) +
    rnorm(m * ni, 0, 0.5)
  dd <- data.frame(y = y, x = x, id = id)
  fit <- statmod(y ~ seg(x, psi ~ random(~ 1 | id), psi = 5),
                 gaussian1_distrib(), dd)
  b <- summary(fit)@tables$mu[[2L]]
  expect_identical(b$head$name, "psi1")
  # it carries the population value of the development and says what
  # develops it, which is the one thing the tables below cannot show: there
  # that number is labelled by the development's intercept
  expect_true(all(is.finite(b$head$estimate)))
  expect_identical(b$head$note, "~ intercept + random")
  expect_false(any(c("beta", "gamma1") %in% b$head$name))
})

test_that("a long block is cut and says how much it hid", {
  tb <- data.frame(name = paste0("v", 1:40), estimate = 1:40 / 10,
                   se = 0.1, statistic = 1, p_value = 0.5,
                   lower = 0, upper = 1, role = "coefficient", source = "",
                   stringsAsFactors = FALSE)
  expect_identical(length(block_rows_shown(tb)), 10L)
  expect_identical(block_rows_shown(tb[1:12, ]), 1:12)
  # a hyperparameter is never among what is dropped: it governs every
  # coefficient under it
  tb$role[1L] <- "estimated"
  expect_true(1L %in% block_rows_shown(tb))
  b <- list(kind = "penalized", label = "Penalized", term = "ridge(x)",
            n_coef = 40L, edf = 5, n_zero = 0L, table = tb, head = NULL,
            components = list())
  # eleven rows are shown, the hyperparameter and ten coefficients, so
  # twenty-nine are not
  expect_output(print_block(b), "29 more, in coef()", fixed = TRUE)
})

test_that("the prefix a term repeats on every row is dropped for printing", {
  expect_identical(drop_common_prefix(c("seg.beta", "seg.gamma1")),
                   c("beta", "gamma1"))
  expect_identical(drop_common_prefix("seg.psi1"), "psi1")
  # nothing in common, nothing dropped
  expect_identical(drop_common_prefix(c("x", "gb")), c("x", "gb"))
  expect_identical(drop_common_prefix(c("a.b", "c.d")), c("a.b", "c.d"))
  # and never everything: a name is left with something to be
  expect_identical(drop_common_prefix(c("a.b", "a.b")), c("b", "b"))
  # ONE piece, not every piece they share: names that agree further along
  # keep what distinguishes them where a reader would look for it
  expect_identical(drop_common_prefix(c("ridge.r.1", "ridge.r.2")),
                   c("r.1", "r.2"))
})

gas_panel_fit <- function(q = 1L) {
  set.seed(21)
  m <- 5
  ni <- 40
  id <- factor(rep(seq_len(m), each = ni))
  tt <- rep(seq_len(ni), m)
  om <- 0.4 + rnorm(m, 0, 0.3)
  y <- numeric(m * ni)
  for (g in seq_len(m)) {
    f <- om[g] / (1 - 0.6)
    for (i in seq_len(ni)) {
      k <- (g - 1) * ni + i
      y[k] <- f + rnorm(1, 0, 1)
      f <- om[g] + 0.15 * (y[k] - f) + 0.6 * f
    }
  }
  dd <- data.frame(y = y, id = id, t = tt)
  statmod(y ~ 0 + gas(p = 1, q = q, omega ~ random(~ 1 | id), by = id,
                      time = t),
          gaussian1_distrib(), dd)
}

test_that("a structural term is a block, and its development a compartment", {
  fit <- gas_panel_fit()
  s <- summary(fit)
  kinds <- vapply(s@tables$mu, `[[`, character(1), "kind")
  b <- s@tables$mu[[which(kinds == "structural")]]
  expect_identical(length(b$components), 1L)
  cp <- b$components[[1L]]
  expect_identical(cp$name, "omega")
  # the term's own table keeps the parameters that are numbers, under the
  # names the term reports: the persistence rides a partial autocorrelation
  # and what is reported is the autoregressive coefficient
  expect_identical(b$table$name, c("alpha1", "beta1"))
  # the hyperparameter sits with the coordinates it shrinks, not in a block
  # of its own carrying nothing else
  expect_true(any(cp$table$role == "estimated"))
  expect_true("(Intercept)" %in% cp$table$name)
  expect_match(cp$lines[[1L]], "^[0-9]+ predictions, sd")

  lines <- utils::capture.output(print(s))
  out <- paste(lines, collapse = "\n")
  # a block is headed by its term and nothing else, at column zero: the
  # call also appears indented under `Call:`, which is not the heading
  expect_true(any(startsWith(lines, "gas(")))
  expect_match(out, "omega  ~ intercept + random", fixed = TRUE)
  expect_false(grepl("omega.random.1", out, fixed = TRUE))
  # nothing reported here carries a test, so the two test columns are not
  # printed at all
  expect_false(grepl("(structural: no design columns)", out, fixed = TRUE))
})

test_that("which parameter a structural quantity belongs to is read off the jacobian", {
  fit <- gas_panel_fit(q = 2L)
  st <- statmod_structural_table(fit)
  expect_true(all(c("component", "position") %in% names(st)))
  # every deviation of the developed level belongs to it
  expect_true(all(st$component[startsWith(st$name, "omega.")] == "omega"))
  # and the first autoregressive coefficient belongs to no single
  # parameter: phi_1 = rho_1 (1 - rho_2) reads two coordinates, so it is
  # filed with the term. The last one is rho_q exactly and does belong to
  # its own, which is the control saying the rule reads the Jacobian rather
  # than the name
  expect_identical(st$component[st$name == "beta1"], "")
  expect_identical(st$component[st$name == "beta2"], "pacf2")
  expect_identical(st$component[st$name == "alpha1"], "alpha1")
})

test_that("a structural term with no development is one flat block", {
  set.seed(22)
  n <- 300
  y <- numeric(n)
  f <- 1
  for (i in seq_len(n)) {
    y[i] <- f + rnorm(1, 0, 1)
    f <- 0.4 + 0.15 * (y[i] - f) + 0.6 * f
  }
  dd <- data.frame(y = y, t = seq_len(n))
  fit <- statmod(y ~ 0 + gas(p = 1, q = 1, time = t), gaussian1_distrib(), dd)
  s <- summary(fit)
  kinds <- vapply(s@tables$mu, `[[`, character(1), "kind")
  b <- s@tables$mu[[which(kinds == "structural")]]
  expect_length(b$components, 0L)
  expect_null(b$head)
  expect_identical(b$table$name, c("omega", "alpha1", "beta1"))
})
