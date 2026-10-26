# The information and the solver are objects. The observed and the expected
# information are defined in distributions7; auto_info() is defined here.

test_that("the three information objects are built, validated and printed", {
  o <- observed_info()
  e <- expected_info()
  a <- auto_info()
  expect_true(S7::S7_inherits(o, Information))
  expect_true(S7::S7_inherits(o, ObservedInfo))
  expect_true(S7::S7_inherits(e, ExpectedInfo))
  expect_true(S7::S7_inherits(a, AutoInfo))
  expect_true(S7::S7_inherits(a, Information))
  # the default method is NULL: the closed form where the family has one and
  # the outer product of the scores otherwise
  expect_null(e@method)
  expect_null(a@method)
  expect_identical(a@repair_floor, 1e-8)
  m <- expected_info(method = distributions7::expect_mc(nsim = 321))
  expect_identical(m@method@nsim, 321)
  sc <- distributions7::expect_integrate(identity = "score")
  expect_identical(auto_info(method = sc)@method, sc)
  expect_identical(auto_info(repair_floor = 1e-10)@repair_floor, 1e-10)
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
  expect_null(info_emethod(o))
  expect_null(info_emethod(e))
  expect_identical(info_emethod(m), m@method)
  expect_identical(info_emethod(auto_info(method = sc)), sc)
})

test_that("the classes of the observed and the expected information are those of distributions7", {
  expect_identical(Information, distributions7::Information)
  expect_identical(ObservedInfo, distributions7::ObservedInfo)
  expect_identical(ExpectedInfo, distributions7::ExpectedInfo)
  expect_false("observed_info" %in% getNamespaceExports("statmodels7"))
  expect_false("expected_info" %in% getNamespaceExports("statmodels7"))
  expect_true("auto_info" %in% getNamespaceExports("statmodels7"))
})

test_that("a wrong object for the method or the floor is an error", {
  expect_error(auto_info(method = "score"), "'method' must be NULL or an ExpectMethod object")
  expect_error(expected_info(method = 3))
  expect_error(auto_info(method = observed_info()), "'method' must be NULL or an ExpectMethod object")
  expect_error(auto_info(repair_floor = 0), "'repair_floor' must be a single positive number")
  expect_error(auto_info(repair_floor = c(1, 2)), "repair_floor")
  expect_error(auto_info(repair_floor = "a"), "repair_floor")
})

test_that("a wrong object for the information is an error", {
  expect_error(reml(information = "auto"),
               "'information' must be an information object: observed_info\\(\\), expected_info\\(\\) or auto_info\\(\\)")
  expect_error(ml(information = 1), "must be an information object")
  expect_error(reml(information = distributions7::expect_opg()),
               "must be an information object")
  expect_error(statmod(y ~ x, distributions7::gaussian1_distrib(),
                       data.frame(x = 1:5, y = c(1, 3, 2, 5, 4)),
                       inner_information = "observed"),
               "'inner_information' must be an information object")
  expect_error(bayesian_vcov(information = "observed"),
               "'information' must be an information object")
  # the former arguments are gone
  expect_error(reml(hessian = "observed"), "unused argument")
})

test_that("aic(), bic() and the variance objects reject the information that is not a matrix", {
  expect_error(aic(information = auto_info()), "auto_info\\(\\) is not admissible")
  expect_error(bic(information = auto_info()), "auto_info\\(\\) is not admissible")
  expect_error(bayesian_vcov(information = auto_info()), "auto_info\\(\\) is not admissible")
  expect_error(frequentist_vcov(information = auto_info()), "not admissible")
  expect_error(unconditional_vcov(information = auto_info()), "not admissible")
  expect_error(aic(information = "auto"), "observed_info\\(\\) or expected_info\\(\\)")
  expect_true(S7::S7_inherits(aic(information = expected_info())@hessian, ExpectedInfo))
  expect_true(S7::S7_inherits(bic()@hessian, ObservedInfo))
  expect_true(S7::S7_inherits(reml()@hessian, ObservedInfo))
  expect_true(S7::S7_inherits(ml()@hessian, ObservedInfo))
  expect_true(S7::S7_inherits(aic()@hessian, ObservedInfo))
  expect_true(S7::S7_inherits(cv()@hessian, ObservedInfo))
  # the information is no longer an argument of iwls(): it is settled by statmod()
  expect_true(S7::S7_inherits(iwls()@hessian, AutoInfo))
})

test_that("auto_info() keeps its method when it is settled", {
  b <- distributions7::expect_integrate(identity = "score")
  m <- iwls()
  # settled against a family without a closed form, the observed information
  # stands and the expected one that takes the step where it cannot keeps the
  # approximation
  r <- iwls_resolve(m, pig_bare_distrib(), auto_info(method = b))
  expect_identical(info_label(r@hessian), "observed")
  expect_identical(info_label(r@backup), "expected")
  expect_identical(r@backup@method, b)
  # settled on a family with a closed form, the expected information keeps it
  g <- iwls_resolve(m, distributions7::gaussian1_distrib(), auto_info(method = b))
  expect_identical(info_label(g@hessian), "expected")
  expect_identical(g@hessian@method, b)
  expect_null(g@backup)
  # the repair floor travels to the method
  rf <- iwls_resolve(m, distributions7::gaussian1_distrib(),
                     auto_info(repair_floor = 1e-5))
  expect_identical(rf@repair_floor, 1e-5)
  expect_identical(iwls_definite_tol(rf), 1e-5)
  # a named information leaves it to the solver's floor
  nm <- iwls_resolve(iwls(solver = solve_chol(floor = 1e-6)), NULL, observed_info())
  expect_true(is.na(nm@repair_floor))
  expect_identical(iwls_definite_tol(nm), 1e-6)
  # an observed method carries no method of the expected information
  expect_null(info_emethod(nm@hessian))
  # the outer criterion resolves with its own method
  rr <- outer_resolve(reml(information = auto_info(method = b)),
                      distributions7::student_t1_distrib())
  expect_identical(info_label(rr@hessian), "expected")
  expect_identical(rr@hessian@method, b)
})

test_that("the method of a criterion is its own", {
  skip_on_cran()
  set.seed(31)
  n <- 60L
  x <- stats::rnorm(n)
  dd <- data.frame(x = x, y = distributions7::distrib_rng(
    distributions7::pig1_distrib(), n, list(mu = exp(1.2 + 0.4 * x), sigma = 0.8)))
  pig <- pig_bare_distrib()
  sc <- distributions7::expect_integrate(identity = "score")
  fb <- statmod(y ~ x, pig, dd,
                outer_criterion = reml(information = expected_info(method = sc)))
  fi <- statmod(y ~ x, pig, dd,
                inner_information = auto_info(method = sc),
                outer_criterion = reml(information = expected_info(method = sc)))
  # the default method of the same criterion reads the outer product of the
  # scores for this family, and the fit exists
  fn <- statmod(y ~ x, pig, dd,
                outer_criterion = reml(information = expected_info()))
  expect_true(fn@converged)
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
  expect_output(print(solve_qr()), "QR solver \\(tol 1e-07, floor 1e-08\\)")
  expect_output(print(solve_svd()), "SVD solver \\(tol default, floor 1e-08\\)")
  expect_output(print(solve_svd(tol = 1e-9)), "SVD solver \\(tol 1e-09")
  expect_output(print(solve_chol()), "Cholesky solver \\(floor 1e-08\\)")
  expect_true(S7::S7_inherits(iwls()@solver, SolveQR))
  expect_error(iwls(solver = "svd"), "'solver' must be a solver object")
  expect_error(iwls(solver = 2), "must be a solver object")
  expect_error(iwls(decomposition = "qr"), "unused argument")
})

test_that("the tolerances of the solvers are arguments with the defaults that were fixed", {
  expect_identical(solve_qr()@tol, 1e-7)
  expect_identical(solve_qr()@floor, 1e-8)
  expect_null(solve_svd()@tol)
  expect_identical(solve_svd()@floor, 1e-8)
  expect_identical(solve_chol()@floor, 1e-8)
  expect_identical(solve_qr(tol = 1e-9, floor = 1e-6)@tol, 1e-9)
  expect_identical(solve_svd(tol = 1e-12)@tol, 1e-12)
  expect_identical(solve_chol(floor = 1e-4)@floor, 1e-4)
  expect_error(solve_qr(tol = 0), "'tol' must be a single positive number")
  expect_error(solve_qr(tol = NULL), "tol")
  expect_error(solve_svd(tol = -1), "tol")
  expect_error(solve_svd(tol = "a"), "tol")
  expect_error(solve_chol(floor = NA_real_), "floor")
  expect_error(solve_qr(floor = c(1e-8, 1e-6)), "floor")
  # the defaults of iwls() are those of its constants
  m <- iwls()
  expect_identical(c(m@decrement_limit, m@switch_after, m@switch_rate),
                   c(1e-6, 10, 0.5))
})

test_that("the pivoted QR tolerance decides which column is dropped", {
  set.seed(8)
  R <- matrix(stats::rnorm(60 * 3), 60, 3)
  # a column that differs from the second by 1e-9 relative is dependent at the
  # default tolerance and independent at a tighter one
  R[, 3] <- R[, 2] * (1 + 1e-9 * stats::rnorm(60))
  pcs <- list(R = R, C = matrix(0, 0, 3))
  u <- c(1, 2, 3)
  expect_length(iwls_solve(pcs, u, "qr")$dropped, 1L)
  expect_length(iwls_solve(pcs, u, "qr", tol = 1e-14)$dropped, 0L)
  expect_identical(iwls_solve(pcs, u, "qr", tol = 1e-7)$delta,
                   iwls_solve(pcs, u, "qr")$delta)
})

test_that("the cut of the singular values sets the numerical rank", {
  set.seed(9)
  R <- matrix(stats::rnorm(40 * 3), 40, 3)
  R[, 3] <- R[, 2] + 1e-6 * stats::rnorm(40)
  pcs <- list(R = R, C = matrix(0, 0, 3))
  u <- c(1, 2, 3)
  default <- iwls_solve(pcs, u, "svd")
  expect_identical(default$rank, 3L)
  loose <- iwls_solve(pcs, u, "svd", tol = 1e-3)
  expect_identical(loose$rank, 2L)
  same <- iwls_solve(pcs, u, "svd", tol = NULL)
  expect_identical(same$delta, default$delta)
})

test_that("the floor of the eigenvalues sets the step on an indefinite matrix", {
  A <- diag(c(4, -1))
  pcs <- list(R = NULL, C = NULL, A = A)
  u <- c(1, 1)
  d8 <- iwls_solve(pcs, u, "chol", floor = 1e-8)$delta
  d4 <- iwls_solve(pcs, u, "chol", floor = 1e-4)$delta
  expect_equal(d8[2], 1 / (4e-8), tolerance = 1e-6)
  expect_equal(d4[2], 1 / (4e-4), tolerance = 1e-6)
  expect_identical(iwls_solve(pcs, u, "chol")$delta, d8)
})

test_that("every information and solver reach the same fit", {
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
                   inner_information = info,
                   inner_optimizer = iwls(solver = sol))
      expect_equal(unlist(coef(f)), unlist(coef(ref)), tolerance = 1e-6,
                   info = paste(info_label(info), solver_how(sol)))
    }
  }
  # the settled information is recorded on the fit
  expect_identical(info_label(ref@methods$smooth@hessian), "expected")
  expect_identical(info_label(ref@methods$inner_information), "auto")
})

test_that("explicit tolerances equal to the defaults give the same fit bit for bit", {
  set.seed(4)
  n <- 100
  dd <- data.frame(x = runif(n))
  dd$y <- rnorm(n, sin(2 * pi * dd$x), 0.4)
  fo <- y ~ s(x, bspline_smooth(k = 8))
  g <- distributions7::gaussian1_distrib()
  a <- statmod(fo, g, dd, outer_criterion = NULL)
  b <- statmod(fo, g, dd, outer_criterion = NULL,
               inner_information = auto_info(repair_floor = 1e-8),
               inner_optimizer = iwls(solver = solve_qr(tol = 1e-7, floor = 1e-8),
                                      decrement_limit = 1e-6,
                                      switch_after = 10, switch_rate = 0.5))
  expect_identical(a@coefficients, b@coefficients)
  expect_identical(a@objective, b@objective)
  s1 <- statmod(fo, g, dd, outer_criterion = NULL,
                inner_optimizer = iwls(solver = solve_svd()))
  s2 <- statmod(fo, g, dd, outer_criterion = NULL,
                inner_optimizer = iwls(solver = solve_svd(floor = 1e-8)))
  expect_identical(s1@coefficients, s2@coefficients)
})
