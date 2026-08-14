# What a Student t fitted to iris found. Every check here is a defect that
# report exposed, and each one had been silent.

f3 <- Sepal.Length ~ Species | sigma ~ Species | nu ~ Species

test_that("the start is the intercept-only MLE, not zero", {
  # statmod_start read distrib_start's result by parameter name. That result
  # is a list of STARTS, each keyed by parameter, so the name matched nothing
  # and every start silently fell to zero on the link scale -- a location of 0
  # for a response centred at 5.84.
  spec <- statmod_spec(Sepal.Length ~ Species,
                       distributions7::gaussian1_distrib(), iris)
  design <- statmod_design(spec)
  obj <- statmod_objective(spec, statmod_hyper_start(spec), design)
  b <- statmod_start(spec, design, obj, NULL)
  cf <- obj$split(b)
  # the gaussian's intercept-only MLE is the sample mean and the sample sd
  expect_equal(cf$mu[1L], mean(iris$Sepal.Length), tolerance = 1e-6)
  expect_equal(exp(cf$sigma[1L]),
               sqrt(mean((iris$Sepal.Length - mean(iris$Sepal.Length))^2)),
               tolerance = 1e-4)
  # and the slopes still start at zero
  expect_equal(cf$mu[-1L], c(0, 0))
})

test_that("the same call gives the same fit", {
  # the intercept-only fit draws its own starts at random, so this inherited
  # them: the same call returned log-likelihoods of -103.49, -112.11 and
  # -111.83 on consecutive runs. The stream is pinned for the length of the
  # call and put back afterwards.
  a <- statmod(f3, distributions7::student_t1_distrib(), iris)
  b <- statmod(f3, distributions7::student_t1_distrib(), iris)
  expect_equal(a@coefficients, b@coefficients, tolerance = 1e-12)
  expect_identical(a@converged, b@converged)

  # and the caller's own stream is untouched
  set.seed(99)
  before <- stats::runif(1)
  set.seed(99)
  invisible(statmod(Sepal.Length ~ Species,
                    distributions7::gaussian1_distrib(), iris))
  expect_equal(stats::runif(1), before)
})

test_that("a fit that cannot be differentiated stops instead of erroring", {
  # `if (score < tol)` met an NA and stopped the run with "missing value where
  # TRUE/FALSE needed", naming neither the iteration nor the cause. A Student
  # t on iris reaches such a point, nu being unidentified there.
  expect_no_error(statmod(f3, distributions7::student_t1_distrib(), iris))
  expect_no_error(statmod(f3, distributions7::student_t2_distrib(), iris))
})

test_that("convergence is the inner method's, not the absence of blocks", {
  # the alternation set converged unconditionally when there was nothing to
  # alternate with, so EVERY model without a kinked penalty reported success
  # whatever the inner fit had done -- which is what let the three defects
  # above go unnoticed
  ok <- statmod(Sepal.Length ~ Species,
                distributions7::gaussian1_distrib(), iris)
  expect_true(ok@converged)

  # one iteration is not enough for this model, and the fit says so
  short <- statmod(f3, distributions7::student_t1_distrib(), iris,
                   inner_optimizer = iwls(maxit = 1L))
  expect_false(short@converged)
})

test_that("a singular information names the direction, not a guess", {
  # the message used to offer two causes -- the fit not having reached a
  # maximum, or two columns carrying the same information -- and on this fit
  # NEITHER was right: the design is full rank and the score is small. The
  # eigenvector of the smallest eigenvalue says what is actually flat.
  fit <- statmod(f3, distributions7::student_t1_distrib(), iris)
  msg <- tryCatch(vcov(fit), error = function(e) conditionMessage(e))
  skip_if(!is.character(msg), "the information happened to be definite here")
  expect_match(msg, "nu")
  expect_match(msg, "edge of its")
  expect_false(grepl("two columns of the design", msg))
})

test_that("a well-identified model is unharmed by any of it", {
  # the same three species, a family whose parameters the data do identify
  fit <- statmod(Sepal.Length ~ Species | sigma ~ Species,
                 distributions7::gaussian1_distrib(), iris)
  expect_true(fit@converged)
  # the group means, exactly
  mu <- predict(fit, "mu")
  for (s in levels(iris$Species)) {
    expect_equal(unique(mu[iris$Species == s]),
                 mean(iris$Sepal.Length[iris$Species == s]),
                 tolerance = 1e-6)
  }
  expect_true(all(is.finite(vcov(fit))))
  expect_identical(nrow(confint(fit)), 6L)
})


test_that("a fit that did not converge says where its parameters ended up", {
  # A lasso at a fixed hyperparameter, a free scale, and a design the model
  # can interpolate: fitting the coefficients shrinks the residuals, which
  # shrinks the scale, which raises the working weights, which makes the
  # penalty count for relatively less, which lets more coefficients in. At
  # 200 observations and 400 columns the scale reached 3.8e-15 and 380 of the
  # 400 coefficients survived, where the same block at a held scale kept the
  # five that were real. The fit reports failure, and now also reports the
  # number that explains it.
  set.seed(7)
  n <- 120L
  p <- 300L
  X <- scale(matrix(stats::rnorm(n * p), n, p), TRUE, FALSE)
  y <- as.numeric(X %*% c(rep(2, 5), rep(0, p - 5))) + stats::rnorm(n)
  dd <- data.frame(y = y - mean(y))
  dd$x <- X
  fit <- statmod(y ~ lasso(x, lambda = 12) - 1 | sigma ~ 1,
                 distributions7::gaussian1_distrib(), dd)
  skip_if(fit@converged, "the runaway did not happen on this platform")
  out <- paste(utils::capture.output(print(fit)), collapse = "\n")
  expect_match(out, "DID NOT CONVERGE")
  expect_match(out, "the parameters it reached")
  expect_match(out, "sigma")
  # and a fit that converged says nothing of the kind
  ok <- statmod(y ~ x1, distributions7::gaussian1_distrib(),
                data.frame(y = stats::rnorm(50), x1 = stats::rnorm(50)))
  expect_true(ok@converged)
  expect_false(grepl("the parameters it reached",
                     paste(utils::capture.output(print(ok)), collapse = "\n")))
})

test_that("an unavailable criterion at the start names its cause", {
  # A gamma model whose nl term sits inside the log link, so the fitted
  # mean is exp(a exp(-r x)): started at a = 0 the Jacobian column of the
  # rate is -a x exp(-r x) and vanishes identically, the inner fit
  # CONVERGES on that flat ridge, the penalized information has no
  # Cholesky factor at any hyperparameter, and the outer search used to
  # die inside the optimizer with "the objective is not finite at the
  # starting value", which names the point and not the cause.
  #
  # The zero start is asked for EXPLICITLY here. It used to be the default
  # and is not one any more -- a term that recomputes its block says where
  # its coefficients begin, through term_coef_start() -- so the same model
  # at the start its own `start =` names is an ordinary fit. What is
  # pinned is the message, not the old default.
  set.seed(8)
  n <- 300
  dd <- data.frame(x = runif(n, 0, 3),
                   g = factor(rep(sprintf("g%d", 1:5), length.out = n)))
  amp <- exp(rnorm(5, log(2), 0.15))
  mu <- amp[dd$g] * exp(-1.2 * dd$x)
  dd$y <- rgamma(n, shape = 20, scale = mu / 20)
  fml <- y ~ nl(~ a * exp(-r * x), a ~ ridge(~g),
                links = list(r = linkfunctions7::log_link()),
                start = list(a = 2, r = 1)) - 1
  expect_error(statmod(fml, distributions7::gamma1_distrib(), dd,
                       start = start_origin()),
               "unavailable at the starting hyperparameters")
  # held, the same model fits and reports where its parameters ended up
  f0 <- suppressWarnings(statmod(fml, distributions7::gamma1_distrib(), dd,
                                 outer_criterion = NULL))
  expect_true(isTRUE(f0@converged))
})

test_that("a smoothing parameter at 1e15 is scale separation, not singularity", {
  # REML sends lambda to ~6e15 on weak signal (counts ~0.14), which is the
  # right answer: the smooth is a straight line. The penalized information
  # then has eigenvalues from ~30 to ~6e15 -- strictly positive definite --
  # and the old min(ev) > 1e-12 * max(ev) read the separation as
  # singularity: no vcov, no edf, summary dead, on a correct fit.
  set.seed(38)
  dP <- data.frame(z = runif(400, -2, 2))
  dP$y <- rpois(400, exp(-2 + 0.8 * sin(dP$z)))
  fit <- statmod(y ~ s(z, k = 8), distributions7::poisson_distrib(), dP)
  expect_true(fit@converged)
  lam <- unlist(fit@hyper$mu)
  skip_if(lam < 1e8, "the criterion did not reach the separation regime here")
  v <- vcov(fit)
  expect_true(all(is.finite(diag(v))))
  # the shrunk directions report variances near zero, not a refusal
  expect_true(min(diag(v), na.rm = TRUE) >= 0)
  e <- fit@edf
  expect_false(is.null(e))
  expect_true(all(is.finite(e$edf)))
  # the smooth is down at its null dimension: a straight line
  expect_lt(e$edf[startsWith(e$term, "s(")], 1.6)
  out <- capture.output(print(summary(fit)))
  expect_true(any(grepl("estimated", out)))

  # and the refusal still refuses what it exists for: a genuinely flat
  # direction, two columns carrying the same information
  set.seed(1)
  dd <- data.frame(x = runif(60))
  dd$x2 <- dd$x
  dd$y <- 1 + dd$x + rnorm(60, 0.2)
  f2 <- statmod(y ~ x + x2, distributions7::gaussian1_distrib(), dd)
  expect_error(vcov(f2), "not positive definite")
})

test_that("the condition estimate stands in for the smallest eigenvalue", {
  # solve_pd() reads the smallest eigenvalue off LAPACK's estimator rather
  # than computing the spectrum. What the test above needs of the estimate
  # is that it separate the same two cases; what THIS one pins is the
  # property that makes that safe -- it errs on the small side, so the
  # verdict is conservative and never accepts a matrix the exact test would
  # refuse, and it errs by a bounded factor rather than by orders.
  est <- function(A) {
    ch <- chol(A)
    an <- max(colSums(abs(A)))
    statmodels7:::chol_rcond_cpp(ch, an) * an
  }
  set.seed(5)
  for (p in c(5L, 40L, 200L)) {
    for (cond in c(1e2, 1e8, 1e15)) {
      Q <- qr.Q(qr(matrix(stats::rnorm(p * p), p, p)))
      ev <- exp(seq(log(1), log(1 / cond), length.out = p)) * 1e3
      A <- Q %*% (ev * t(Q))
      A <- (A + t(A)) / 2
      lo <- min(eigen(A, symmetric = TRUE)$values)
      r <- est(A) / lo
      info <- sprintf("p = %d, condition %.0e", p, cond)
      expect_lt(r, 1 + 1e-8)          # never optimistic
      expect_gt(r, 1 / sqrt(p))       # and the bound the norms give
      expect_gt(r, 0.2)               # measured: 0.29 at p = 300
    }
  }

  # the inverse is the one the eigendecomposition returned, to the accuracy
  # the conditioning allows -- two routes to one matrix, so a tolerance and
  # not an identity
  set.seed(9)
  p <- 120L
  Q <- qr.Q(qr(matrix(stats::rnorm(p * p), p, p)))
  A <- Q %*% (exp(seq(log(1e4), log(1e-2), length.out = p)) * t(Q))
  A <- (A + t(A)) / 2
  e <- eigen(A, symmetric = TRUE)
  V_eig <- e$vectors %*% (t(e$vectors) / e$values)
  V_new <- statmodels7:::solve_pd(A, "x")
  expect_equal(V_new, V_eig, tolerance = 1e-6)
  expect_equal(A %*% V_new, diag(p), tolerance = 1e-8, ignore_attr = TRUE)
  # symmetric, because a variance matrix is
  expect_identical(V_new, t(V_new))
})

test_that("the stopping rule is dimensionless in the response", {
  # the score of a location equation carries the units 1/y, so the old
  # absolute threshold on score/n had a floor that CROSSED it three decades
  # down: at y ~ 1e-3 the run ended at 1.37e-6 against 1e-6, at the
  # optimum, and reported failure. Normalized by the curvature the rule
  # survives the rescaling in both directions.
  set.seed(31)
  n <- 300
  dS <- data.frame(z = runif(n, -2, 2))
  base_y <- 2 + sin(1.5 * dS$z) + rnorm(n, sd = 0.3)
  # The criterion is asked for the EXPECTED information here, not because the
  # rule under test needs it but because the default's is the observed one and
  # its plateau is what would be under test instead: at y * 1e4 the criterion
  # is flat from the optimum at 3.2e-8 down to 1e-19, Newton walks the whole
  # of it, and the fitted function is the same to five decimals while the
  # smoothing parameter reported is twelve decades past the optimum. That is a
  # property of the criterion's surface and is recorded as an open item; what
  # this test is about is the INNER rule surviving a rescaling of the response.
  fits <- lapply(c(1e-3, 1, 1e4), function(sc) {
    dS$y <- base_y * sc
    statmod(y ~ s(z, k = 10), distributions7::gaussian1_distrib(), dS,
            outer_criterion = reml("expected"))
  })
  for (f in fits) expect_true(f@converged)
  # the same fit at every scale where the outer search starts near its
  # optimum; at y * 1e4 the smoothing parameter must travel eight decades
  # from the shared start, so there the assertion is the recovery, not the
  # function to four decimals
  m1 <- predict(fits[[1]])$mu / 1e-3
  m2 <- predict(fits[[2]])$mu
  m3 <- predict(fits[[3]])$mu / 1e4
  expect_equal(m1, m2, tolerance = 1e-4)
  truth <- 2 + sin(1.5 * dS$z)
  expect_gt(cor(m3, truth), 0.99)
})
