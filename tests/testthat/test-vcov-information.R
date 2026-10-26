pig_data <- function(n = 200L, seed = 1L) {
  set.seed(seed)
  x <- stats::rnorm(n)
  mu <- exp(1.2 + 0.4 * x)
  # the toolkit's own generator rather than an external one: a package
  # outside the toolkit is not a dependency here, even in Suggests, and
  # naming one in a test is an unstated dependency --as-cran reports. The
  # parametrizations differ and the draws do not -- (mean, shape = mu/0.8)
  # is Var = 0.8 mu^2, which is phi = 0.8/mu where Var = phi mu^3, and both
  # run the same Michael-Schucany-Haas transformation on the same stream:
  # measured over 2e5 draws at three means, the two agree to every printed
  # digit of the mean and the variance.
  lam <- distributions7::distrib_rng(distributions7::invgauss1_distrib(), n,
                                     list(mu = mu, phi = 0.8 / mu))
  data.frame(y = stats::rpois(n, lam), x = x)
}

test_that("the standard errors are the observed ones for a smooth family, whatever the inner information", {
  dd <- pig_data(120L)
  fam <- list(pig_bare_distrib(), distributions7::pig1_distrib())
  for (d in fam) {
    f_def <- statmod(y ~ x, d, dd)
    f_exp <- statmod(y ~ x, d, dd, inner_information = expected_info())
    f_obs <- statmod(y ~ x, d, dd, inner_information = observed_info())
    ## the observed information inverted at each fit's own estimates
    for (f in list(f_def, f_exp, f_obs)) {
      expect_identical(vcov(f), vcov(f, vcov_type = bayesian_vcov(
        information = observed_info())))
    }
  }
  # one that writes its expected information out has the same default
  set.seed(2)
  gd <- data.frame(x = stats::runif(120))
  gd$y <- 1 + 2 * gd$x + stats::rnorm(120, sd = 0.4)
  f_g <- statmod(y ~ x, distributions7::gaussian1_distrib(), gd)
  expect_identical(info_label(f_g@methods$smooth@hessian), "expected")
  expect_identical(vcov(f_g), vcov(f_g, vcov_type = bayesian_vcov(
    information = observed_info())))
})

test_that("the observed standard errors are the inverse of the observed Hessian", {
  dd <- pig_data(150L)
  fit <- statmod(y ~ x, distributions7::pig1_distrib(), dd,
                 outer_criterion = NULL)
  H <- -hessian(fit, information = observed_info())
  expect_equal(unname(vcov(fit)), unname(solve(H)), tolerance = 1e-8)
  Ie <- -hessian(fit, information = expected_info())
  expect_equal(unname(vcov(fit, vcov_type = bayesian_vcov(
    information = expected_info()))), unname(solve(Ie)), tolerance = 1e-8)
  expect_false(isTRUE(all.equal(vcov(fit), vcov(fit, vcov_type = bayesian_vcov(
    information = expected_info())))))
})

test_that("each analytic convention takes the information, and the abstract class does not", {
  for (ctor in list(bayesian_vcov, frequentist_vcov, unconditional_vcov)) {
    n <- ctor()
    o <- ctor(information = observed_info())
    e <- ctor(information = expected_info())
    expect_null(n@information)
    expect_true(S7::S7_inherits(o@information, ObservedInfo))
    expect_true(S7::S7_inherits(e@information, ExpectedInfo))
    expect_true(S7::S7_inherits(e, vcov_method))
  }
  expect_identical(names(vcov_method@properties), "name")
})

test_that("the information reaches the family through the method", {
  # The two methods differ, which is what says the argument is not being
  # swallowed. The integration is the dear route, so this runs small. The
  # family must approximate its expected information for the method to mean
  # anything: pig1 computes its own exactly since distributions7 0.65.0.
  dd <- pig_data(60L)
  fit <- statmod(y ~ x, pig_bare_distrib(), dd)
  V_opg <- vcov(fit, vcov_type = bayesian_vcov(
    information = expected_info(method = distributions7::expect_opg())))
  V_sc <- vcov(fit, vcov_type = bayesian_vcov(
    information = expected_info(method = distributions7::expect_integrate(
      identity = "score"))))
  expect_false(isTRUE(all.equal(V_opg, V_sc)))
  expect_identical(dimnames(V_opg), dimnames(V_sc))
  # the default method of this family is the outer product of the scores
  V_null <- vcov(fit, vcov_type = bayesian_vcov(information = expected_info()))
  expect_identical(V_null, V_opg)
  # a bad value is rejected by the constructor
  expect_error(expected_info(method = "nonesuch"), "ExpectMethod object")
  # and the arguments that moved into the object are not arguments of the methods
  expect_error(vcov(fit, expected = TRUE), "vcov\\(\\) takes no argument 'expected'")
  expect_error(vcov(fit, method = distributions7::expect_opg()),
               "takes no argument 'method'")
  expect_error(confint(fit, expected = TRUE), "takes no argument 'expected'")
  expect_error(summary(fit, expected = TRUE), "takes no argument 'expected'")
})

test_that("summary and confint read the information of vcov_type", {
  dd <- pig_data(150L)
  fit <- statmod(y ~ x, distributions7::pig1_distrib(), dd)
  ve <- bayesian_vcov(information = expected_info())
  s_def <- summary(fit)
  s_obs <- summary(fit, vcov_type = bayesian_vcov(information = observed_info()))
  s_exp <- summary(fit, vcov_type = ve)

  # a StatmodSummary's @tables is keyed by distribution parameter, and each
  # entry is a LIST OF BLOCKS (one per term, the shape term_block_kind()
  # names), not a flat data.frame: the numbers are one level further in, at
  # each block's own $table.
  col <- function(s, which) unlist(lapply(s@tables, function(blocks)
    lapply(blocks, function(b) b$table[[which]])), use.names = FALSE)
  se <- function(s) col(s, "se")
  est <- function(s) col(s, "estimate")
  expect_equal(se(s_def), se(s_obs))
  expect_false(isTRUE(all.equal(se(s_def), se(s_exp))))
  # the ESTIMATES do not move: only which curvature is read
  expect_equal(est(s_obs), est(s_exp))
  ci_e <- confint(fit, vcov_type = ve)
  expect_equal(ci_e$se, sqrt(diag(vcov(fit, vcov_type = ve))), ignore_attr = TRUE)
})

test_that("the closed form and the outer product of the scores reach the same fit", {
  # WHAT THE DEFAULT IS ALLOWED TO CHANGE: the step, not the answer. The score
  # is exact, so a scoring iteration driven by either matrix has the same fixed
  # point. Small n, because the integration sums over the support per row.
  dd <- pig_data(60L)
  f_cf <- statmod(y ~ x, distributions7::pig1_distrib(), dd,
                  inner_information = expected_info())
  f_opg <- statmod(y ~ x, distributions7::pig1_distrib(), dd,
                   inner_information = expected_info(
                     method = distributions7::expect_opg()))
  f_sc <- statmod(y ~ x, distributions7::pig1_distrib(), dd,
                  inner_information = expected_info(
                    method = distributions7::expect_integrate(identity = "score")))
  expect_equal(unlist(coef(f_cf)), unlist(coef(f_sc)), tolerance = 1e-4)
  expect_equal(unlist(coef(f_cf)), unlist(coef(f_opg)), tolerance = 1e-4)
  expect_equal(as.numeric(logLik(f_cf)), as.numeric(logLik(f_sc)),
               tolerance = 1e-6)
  expect_equal(as.numeric(logLik(f_cf)), as.numeric(logLik(f_opg)),
               tolerance = 1e-6)
})

test_that("a family with a closed form fits alike under an explicit method", {
  set.seed(7)
  gd <- data.frame(x = stats::runif(150))
  gd$y <- 1 + 2 * gd$x + stats::rnorm(150, sd = 0.4)
  a <- statmod(y ~ x, distributions7::gaussian1_distrib(), gd)
  # an explicit method is used also by a family that has the closed form
  b <- statmod(y ~ x, distributions7::gaussian1_distrib(), gd,
               inner_information = auto_info(
                 method = distributions7::expect_integrate(identity = "score")))
  expect_equal(unlist(coef(a)), unlist(coef(b)), tolerance = 1e-6)
  expect_equal(vcov(a), vcov(b), tolerance = 1e-6)
})

test_that("hessian() takes the information as an object", {
  dd <- pig_data(80L)
  fit <- statmod(y ~ x, distributions7::pig1_distrib(), dd)
  expect_identical(hessian(fit), hessian(fit, information = observed_info()))
  expect_false(isTRUE(all.equal(hessian(fit),
                                hessian(fit, information = expected_info()))))
  expect_error(hessian(fit, information = auto_info()), "not admissible")
  expect_error(hessian(fit, information = "expected"), "information object")
  expect_error(hessian(fit, expected = TRUE), "takes no argument 'expected'")
})

test_that("the default information of a variance object is resolved by the family", {
  # the rule of distributions7: observed where every parameter is twice
  # differentiable, expected (the whole matrix) otherwise
  d7 <- asNamespace("distributions7")
  for (d in list(distributions7::gaussian1_distrib(),
                 distributions7::gamma1_distrib(),
                 distributions7::pig1_distrib(),
                 distributions7::laplace_distrib(),
                 distributions7::laplace2_distrib(),
                 distributions7::enet_distrib(),
                 pig_bare_distrib())) {
    mine <- default_vcov_information(d)
    theirs <- d7$default_information(d)
    expect_identical(class(mine), class(theirs), label = d@distrib_name)
  }
  expect_true(S7::S7_inherits(default_vcov_information(distributions7::gaussian1_distrib()),
                              ObservedInfo))
  expect_true(S7::S7_inherits(default_vcov_information(distributions7::laplace_distrib()),
                              ExpectedInfo))
  # a variance object that names no information, or one, resolves accordingly
  expect_true(S7::S7_inherits(
    vcov_information(bayesian_vcov(), distributions7::laplace_distrib()), ExpectedInfo))
  expect_true(S7::S7_inherits(
    vcov_information(bayesian_vcov(information = observed_info()),
                     distributions7::laplace_distrib()), ObservedInfo))
  expect_true(S7::S7_inherits(
    vcov_information(unconditional_vcov(), distributions7::gamma1_distrib()),
    ObservedInfo))
})

test_that("the Laplace and the Laplace2 fits have finite standard errors by default", {
  set.seed(9)
  n <- 150
  dd <- data.frame(x = stats::runif(n))
  dd$y <- 1 + 2 * dd$x + stats::rexp(n) * sample(c(-1, 1), n, TRUE)
  for (d in list(distributions7::laplace_distrib(), distributions7::laplace2_distrib())) {
    fit <- statmod(y ~ x, d, dd)
    for (ty in list(bayesian_vcov(), frequentist_vcov(), unconditional_vcov())) {
      se <- sqrt(diag(vcov(fit, vcov_type = ty)))
      expect_true(all(is.finite(se)), label = paste(d@distrib_name, ty@name))
    }
    ci <- confint(fit)
    expect_true(all(is.finite(c(ci$lower, ci$upper))), label = d@distrib_name)
    s <- summary(fit)
    # the matrix is named when it is not the observed one
    expect_output(print(s), "bayesian variance on the expected information")
    expect_true(S7::S7_inherits(s@information, ExpectedInfo))
    # and it is the expected information that was inverted
    H <- -hessian(fit, information = expected_info())
    expect_equal(unname(vcov(fit)), unname(solve(H)), tolerance = 1e-8)
  }
  # naming the observed information is the other choice
  fit <- statmod(y ~ x, distributions7::laplace_distrib(), dd)
  expect_error(suppressWarnings(vcov(fit, vcov_type = bayesian_vcov(
    information = observed_info()))))
  expect_false(anyNA(vcov(fit)))
})

test_that("the summary of a smooth family names no information", {
  dd <- pig_data(120L)
  fit <- statmod(y ~ x, distributions7::pig1_distrib(), dd)
  expect_true(S7::S7_inherits(summary(fit)@information, ObservedInfo))
  expect_output(print(summary(fit)), "bayesian variance\n")
  out <- utils::capture.output(print(summary(
    fit, vcov_type = bayesian_vcov(information = expected_info()))))
  expect_true(any(grepl("bayesian variance on the expected information", out)))
  # the printed variance object names an information only when one was given
  expect_identical(utils::capture.output(print(bayesian_vcov())), "<vcov> bayesian")
  expect_output(print(frequentist_vcov(information = expected_info())),
                "frequentist [(]expected information[)]")
})
