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

test_that("fit_expected asks the family as well as the fit", {
  dd <- pig_data(120L)
  # a family with no closed expected information: the step may use an
  # approximation of it, the report must not
  f_pig <- statmod(y ~ x, distributions7::pig1_distrib(), dd)
  expect_identical(f_pig@methods$smooth@hessian, "expected")
  expect_false(statmodels7:::fit_expected(f_pig))

  # one that writes it out keeps it
  set.seed(2)
  gd <- data.frame(x = stats::runif(120))
  gd$y <- 1 + 2 * gd$x + stats::rnorm(120, sd = 0.4)
  f_g <- statmod(y ~ x, distributions7::gaussian1_distrib(), gd)
  expect_true(statmodels7:::fit_expected(f_g))

  # and an explicit "observed" is followed whatever the family
  f_obs <- statmod(y ~ x, distributions7::gaussian1_distrib(), gd,
                   inner = iwls(hessian = "observed"))
  expect_false(statmodels7:::fit_expected(f_obs))
})

test_that("the default report is the observed information where the expected one is approximated", {
  dd <- pig_data(150L)
  fit <- statmod(y ~ x, distributions7::pig1_distrib(), dd)
  expect_equal(vcov(fit), vcov(fit, expected = FALSE))
  expect_false(isTRUE(all.equal(vcov(fit), vcov(fit, expected = TRUE))))
})

test_that("vcov takes approx and it reaches the family", {
  # The two readings of the identity differ, which is what says the argument
  # is not being swallowed. bartlett is the dear route, so this runs small.
  dd <- pig_data(60L)
  fit <- statmod(y ~ x, distributions7::pig1_distrib(), dd)
  V_opg <- vcov(fit, expected = TRUE, approx = "opg")
  V_bar <- vcov(fit, expected = TRUE, approx = "bartlett")
  expect_false(isTRUE(all.equal(V_opg, V_bar)))
  expect_identical(dimnames(V_opg), dimnames(V_bar))
  # a bad value is refused by name rather than passed on
  expect_error(vcov(fit, approx = "nonesuch"))
})

test_that("summary passes expected and approx through to vcov", {
  dd <- pig_data(150L)
  fit <- statmod(y ~ x, distributions7::pig1_distrib(), dd)
  s_def <- summary(fit)
  s_obs <- summary(fit, expected = FALSE)
  s_exp <- summary(fit, expected = TRUE)

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
})

test_that("the cheap and the dear expected information reach the same fit", {
  # WHAT THE DEFAULT IS ALLOWED TO CHANGE: the step, not the answer. The score
  # is exact, so a scoring iteration driven by either matrix has the same fixed
  # point. Small n, because the bartlett route sums over the support per row.
  dd <- pig_data(60L)
  f_opg <- statmod(y ~ x, distributions7::pig1_distrib(), dd)
  f_bar <- statmod(y ~ x, distributions7::pig1_distrib(), dd,
                   inner = iwls(approx = "bartlett"))
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
               inner = iwls(approx = "bartlett"))
  expect_equal(unlist(coef(a)), unlist(coef(b)))
  expect_equal(vcov(a), vcov(b))
})
