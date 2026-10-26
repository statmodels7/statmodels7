# A criterion on the expected information of a family without a closed form.
#
# distributions7::expected_info() with method = NULL reads the closed form of
# the family where it has one and the outer product of the scores at the data
# otherwise. The criterion built on the second reads the same method object for
# the matrix and for the derivatives that its exact gradient and Hessian need,
# so it is fitted like any other. The family these tests use is
# pig_bare_distrib(), from helper-pigbare.R.

test_that("the method of the derivatives follows the family", {
  # a family that writes its expected information out keeps NULL
  expect_null(resolve_emethod(distributions7::gamma1_distrib(), NULL))
  expect_null(resolve_emethod(distributions7::negbin2_distrib(), NULL))
  expect_null(resolve_emethod(distributions7::pig1_distrib(), NULL))
  # a truncated family registers its quadrature and its derivative
  expect_null(resolve_emethod(
    distributions7::truncated(distributions7::gaussian1_distrib(), lower = -1), NULL))
  # a family with a density only is read as the outer product of the scores
  expect_true(S7::S7_inherits(resolve_emethod(pig_bare_distrib(), NULL),
                              distributions7::ExpectOPG))
  expect_true(S7::S7_inherits(
    resolve_emethod(distributions7::fixed(pig_bare_distrib(), sigma = 2), NULL),
    distributions7::ExpectOPG))
  # an object is returned as given, for any family
  for (m in list(distributions7::expect_opg(),
                 distributions7::expect_integrate(identity = "score"),
                 distributions7::expect_mc(nsim = 50))) {
    expect_identical(resolve_emethod(pig_bare_distrib(), m), m)
    expect_identical(resolve_emethod(distributions7::gamma1_distrib(), m), m)
  }
  # and the caller's random stream is left where it was
  set.seed(3)
  before <- .Random.seed
  resolve_emethod(pig_bare_distrib(), NULL)
  expect_identical(.Random.seed, before)
})

test_that("the derivative probes follow the method", {
  pig <- pig_bare_distrib()
  expect_false(expected_deriv_ok(pig, NULL))
  expect_true(expected_deriv_ok(pig, distributions7::expect_opg()))
  expect_true(expected_deriv2_ok(pig, distributions7::expect_opg()))
  expect_true(expected_deriv_ok(distributions7::gamma1_distrib(), NULL))
})

test_that("a criterion on the expected information of a density-only family is fitted", {
  skip_on_cran()
  set.seed(11)
  n <- 120
  dd <- data.frame(x = runif(n))
  mu <- exp(0.8 + 0.5 * sin(2 * pi * dd$x))
  dd$y <- distributions7::distrib_rng(distributions7::pig1_distrib(), n,
                                      list(mu = mu, sigma = 0.5))
  fo <- y ~ s(x, bspline_smooth(k = 8))
  fit <- statmod(fo, pig_bare_distrib(), dd,
                 outer_criterion = reml(information = expected_info()))
  ct <- statmod_certificate(fit)
  expect_identical(ct$curvature, "analytic")
  expect_identical(ct$state, "converged")
  # the outer product of the scores is the method that the default selects,
  # so naming it gives the same fit
  fo2 <- statmod(fo, pig_bare_distrib(), dd,
                 outer_criterion = reml(information = expected_info(
                   method = distributions7::expect_opg())))
  expect_identical(fit@hyper, fo2@hyper)
  expect_identical(fit@coefficients, fo2@coefficients)
  # a model with nothing for the criterion to estimate is fitted as well
  expect_no_error(statmod(y ~ x, pig_bare_distrib(), dd,
                          outer_criterion = reml(information = expected_info(),
                                                 marginal = "none")))
  # while the shipped family, with its closed form, is fitted on the expected
  # route with an analytic certificate
  fe <- statmod(fo, distributions7::pig1_distrib(), dd,
                outer_criterion = reml(information = expected_info()))
  ce <- statmod_certificate(fe)
  expect_identical(ce$curvature, "analytic")
  expect_identical(ce$state, "converged")
})
