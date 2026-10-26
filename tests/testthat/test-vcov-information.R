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

test_that("the standard errors are the observed ones for every fit, whatever the inner information", {
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
  H <- -hessian(fit, expected = FALSE)
  expect_equal(unname(vcov(fit)), unname(solve(H)), tolerance = 1e-8)
  Ie <- -hessian(fit, expected = TRUE)
  expect_equal(unname(vcov(fit, vcov_type = bayesian_vcov(
    information = expected_info()))), unname(solve(Ie)), tolerance = 1e-8)
  expect_false(isTRUE(all.equal(vcov(fit), vcov(fit, vcov_type = bayesian_vcov(
    information = expected_info())))))
})

test_that("each analytic convention takes the information, and the abstract class does not", {
  for (ctor in list(bayesian_vcov, frequentist_vcov, unconditional_vcov)) {
    o <- ctor()
    e <- ctor(information = expected_info())
    expect_true(S7::S7_inherits(o@information, ObservedInfo))
    expect_true(S7::S7_inherits(e@information, ExpectedInfo))
    expect_true(S7::S7_inherits(e, vcov_method))
  }
  expect_identical(names(vcov_method@properties), "name")
})

test_that("the information reaches the family through the approximation", {
  # The two readings of the identity differ, which is what says the argument
  # is not being swallowed. bartlett is the dear route, so this runs small.
  # The family must approximate its expected information for approx to mean
  # anything: pig1 computes its own exactly since distributions7 0.65.0.
  dd <- pig_data(60L)
  fit <- statmod(y ~ x, pig_bare_distrib(), dd)
  V_opg <- vcov(fit, vcov_type = bayesian_vcov(
    information = expected_info(approx = distributions7::expect_opg())))
  V_bar <- vcov(fit, vcov_type = bayesian_vcov(
    information = expected_info(approx = distributions7::expect_bartlett())))
  expect_false(isTRUE(all.equal(V_opg, V_bar)))
  expect_identical(dimnames(V_opg), dimnames(V_bar))
  # a bad value is rejected by the constructor
  expect_error(expected_info(approx = "nonesuch"), "ExpectApprox object")
  # and the arguments that moved into the object are gone from the methods
  expect_error(vcov(fit, expected = TRUE), "vcov\\(\\) takes no argument 'expected'")
  expect_error(vcov(fit, approx = distributions7::expect_opg()),
               "takes no argument 'approx'")
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

test_that("the cheap and the dear expected information reach the same fit", {
  # WHAT THE DEFAULT IS ALLOWED TO CHANGE: the step, not the answer. The score
  # is exact, so a scoring iteration driven by either matrix has the same fixed
  # point. Small n, because the bartlett route sums over the support per row.
  dd <- pig_data(60L)
  f_opg <- statmod(y ~ x, distributions7::pig1_distrib(), dd,
                   inner_information = expected_info())
  f_bar <- statmod(y ~ x, distributions7::pig1_distrib(), dd,
                   inner_information = expected_info(
                     approx = distributions7::expect_bartlett()))
  expect_equal(unlist(coef(f_opg)), unlist(coef(f_bar)), tolerance = 1e-4)
  expect_equal(as.numeric(logLik(f_opg)), as.numeric(logLik(f_bar)),
               tolerance = 1e-6)
})

test_that("a family with a closed form fits identically under either approx", {
  set.seed(7)
  gd <- data.frame(x = stats::runif(150))
  gd$y <- 1 + 2 * gd$x + stats::rnorm(150, sd = 0.4)
  a <- statmod(y ~ x, distributions7::gaussian1_distrib(), gd)
  b <- statmod(y ~ x, distributions7::gaussian1_distrib(), gd,
               inner_information = auto_info(
                 approx = distributions7::expect_bartlett()))
  expect_equal(unlist(coef(a)), unlist(coef(b)))
  expect_equal(vcov(a), vcov(b))
})
