# The information and the solver are objects, as the approximation of an
# expectation is in distributions7.

test_that("the three information objects are built, validated and printed", {
  o <- observed_info()
  e <- expected_info()
  a <- auto_info()
  expect_true(S7::S7_inherits(o, Information))
  expect_true(S7::S7_inherits(o, ObservedInfo))
  expect_true(S7::S7_inherits(e, ExpectedInfo))
  expect_true(S7::S7_inherits(a, AutoInfo))
  # the default approximation is the outer product of the scores
  expect_true(S7::S7_inherits(e@approx, distributions7::ExpectOPG))
  expect_true(S7::S7_inherits(a@approx, distributions7::ExpectOPG))
  m <- expected_info(approx = distributions7::expect_mc(nsim = 321))
  expect_identical(m@approx@nsim, 321)
  expect_identical(auto_info(approx = distributions7::expect_bartlett())@approx,
                   distributions7::expect_bartlett())
  expect_output(print(o), "Observed information")
  expect_output(print(e), "Expected information .*outer product")
  expect_output(print(m), "Monte Carlo \\(nsim = 321\\)")
  expect_output(print(a), "chosen against the family")
  # helpers
  expect_identical(c(info_label(o), info_label(e), info_label(a)),
                   c("observed", "expected", "auto"))
  expect_true(info_expected(e))
  expect_false(info_expected(o))
  expect_true(info_auto(a))
  expect_true(S7::S7_inherits(info_approx(o), distributions7::ExpectOPG))
})

test_that("a string or a wrong object for the approximation is an error", {
  expect_error(expected_info(approx = "opg"), "no longer accepted. Use expect_opg\\(\\)")
  expect_error(expected_info(approx = "mc"), "Use expect_mc\\(nsim = \\.\\.\\.\\)")
  expect_error(auto_info(approx = "bartlett"), "Use expect_bartlett\\(\\)")
  expect_error(expected_info(approx = 3), "must be an expectation object")
  expect_error(auto_info(approx = observed_info()), "must be an expectation object")
})

test_that("a string for the information is an error that names the object", {
  expect_error(iwls(information = "observed"),
               "string \"observed\" is no longer accepted. Use observed_info\\(\\)")
  expect_error(iwls("expected"), "Use expected_info\\(\\)")
  expect_error(reml(information = "auto"), "Use auto_info\\(\\)")
  expect_error(ml(information = "expected"), "Use expected_info\\(\\)")
  expect_error(aic(information = "observed"), "Use observed_info\\(\\)")
  expect_error(bic(information = "expected"), "Use expected_info\\(\\)")
  expect_error(iwls(information = 1), "must be an information object")
  expect_error(reml(information = distributions7::expect_opg()),
               "must be an information object")
  # the former arguments are gone
  expect_error(iwls(hessian = "observed"), "unused argument")
  expect_error(iwls(approx = distributions7::expect_opg()), "unused argument")
  expect_error(reml(hessian = "observed"), "unused argument")
})

test_that("aic() and bic() accept two objects and reject the one that is not a matrix", {
  expect_error(aic(information = auto_info()), "auto_info\\(\\) is not admissible")
  expect_error(bic(information = auto_info()), "auto_info\\(\\) is not admissible")
  expect_error(aic(information = "auto"), "Use observed_info\\(\\) or expected_info\\(\\)")
  expect_true(S7::S7_inherits(aic(information = expected_info())@hessian, ExpectedInfo))
  expect_true(S7::S7_inherits(bic()@hessian, ObservedInfo))
  # the defaults are those of the string interface
  expect_true(S7::S7_inherits(reml()@hessian, ObservedInfo))
  expect_true(S7::S7_inherits(ml()@hessian, ObservedInfo))
  expect_true(S7::S7_inherits(aic()@hessian, ObservedInfo))
  expect_true(S7::S7_inherits(cv()@hessian, ObservedInfo))
  expect_true(S7::S7_inherits(iwls()@hessian, AutoInfo))
})

test_that("iwls() keeps the approximation of its information", {
  b <- distributions7::expect_bartlett()
  m <- iwls(information = auto_info(approx = b))
  expect_identical(m@approx, b)
  # resolved against a family without a closed form, the observed information
  # stands and the expected one that takes the step where it cannot keeps the
  # approximation
  r <- iwls_resolve(m, pig_bare_distrib())
  expect_identical(info_label(r@hessian), "observed")
  expect_true(r@fallback)
  expect_identical(r@approx, b)
  # resolved on a family with a closed form, the expected information keeps it
  g <- iwls_resolve(m, distributions7::gaussian1_distrib())
  expect_identical(info_label(g@hessian), "expected")
  expect_identical(g@hessian@approx, b)
  # an observed method carries the default
  expect_true(S7::S7_inherits(iwls(information = observed_info())@approx,
                              distributions7::ExpectOPG))
  # the outer criterion resolves with its own approximation
  rr <- outer_resolve(reml(information = auto_info(approx = b)),
                      distributions7::student_t1_distrib())
  expect_identical(info_label(rr@hessian), "expected")
  expect_identical(rr@hessian@approx, b)
})

test_that("the approximation of a criterion is its own", {
  skip_on_cran()
  set.seed(31)
  n <- 60L
  x <- stats::rnorm(n)
  dd <- data.frame(x = x, y = distributions7::distrib_rng(
    distributions7::pig1_distrib(), n, list(mu = exp(1.2 + 0.4 * x), sigma = 0.8)))
  pig <- pig_bare_distrib()
  expect_error(statmod(y ~ x, pig, dd,
                       outer_criterion = reml(information = expected_info())),
               "pig bare")
  fb <- statmod(y ~ x, pig, dd,
                outer_criterion = reml(information = expected_info(
                  approx = distributions7::expect_bartlett())))
  fi <- statmod(y ~ x, pig, dd,
                inner_optimizer = iwls(information = auto_info(
                  approx = distributions7::expect_bartlett())),
                outer_criterion = reml(information = expected_info(
                  approx = distributions7::expect_bartlett())))
  expect_equal(unlist(coef(fb)), unlist(coef(fi)), tolerance = 1e-5)
  expect_equal(as.numeric(logLik(fb)), as.numeric(logLik(fi)), tolerance = 1e-6)
})

test_that("the solvers are objects, validated and printed", {
  expect_true(S7::S7_inherits(solve_qr(), SolveQR))
  expect_true(S7::S7_inherits(solve_svd(), SolveSVD))
  expect_true(S7::S7_inherits(solve_chol(), SolveChol))
  for (s in list(solve_qr(), solve_svd(), solve_chol())) {
    expect_true(S7::S7_inherits(s, IwlsSolver))
  }
  expect_identical(c(solver_how(solve_qr()), solver_how(solve_svd()),
                     solver_how(solve_chol())), c("qr", "svd", "chol"))
  expect_output(print(solve_qr()), "QR solver")
  expect_output(print(solve_svd()), "SVD solver")
  expect_output(print(solve_chol()), "Cholesky solver")
  expect_true(S7::S7_inherits(iwls()@solver, SolveQR))
  expect_error(iwls(solver = "svd"), "string \"svd\" is no longer accepted. Use solve_svd\\(\\)")
  expect_error(iwls(solver = 2), "must be a solver object")
  expect_error(iwls(decomposition = "qr"), "unused argument")
  expect_output(print(iwls(solver = solve_chol())), "auto information, chol")
})

test_that("every information and solver of iwls() reaches the same fit", {
  set.seed(3)
  n <- 150
  dd <- data.frame(x = runif(n))
  dd$y <- rnorm(n, 1 + sin(2 * pi * dd$x), 0.5)
  fo <- y ~ s(x, bspline_smooth(k = 8))
  g <- distributions7::gaussian1_distrib()
  ref <- statmod(fo, g, dd, outer_criterion = NULL)
  for (info in list(observed_info(), expected_info(), auto_info())) {
    for (sol in list(solve_qr(), solve_svd(), solve_chol())) {
      f <- statmod(fo, g, dd, outer_criterion = NULL,
                   inner_optimizer = iwls(information = info, solver = sol))
      expect_equal(unlist(coef(f)), unlist(coef(ref)), tolerance = 1e-6,
                   info = paste(info_label(info), solver_how(sol)))
    }
  }
  # the settled information is recorded on the fit
  expect_identical(info_label(ref@methods$smooth@hessian), "expected")
})
