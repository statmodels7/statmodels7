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
                   inner_method = iwls(maxit = 1L))
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
