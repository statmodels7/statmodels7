# A criterion on the expected information needs an expectation.
#
# A family that does not write its expected information out returns, under the
# default approx = expect_opg(), the outer product of the scores at the data, which
# depends on the response. A criterion built on it is neither the Laplace
# approximation nor its Fisher variant and has no exact outer derivatives, so
# statmod() rejects it and names the observed route, which is exact. No
# shipped family has done so since distributions7 0.65.0, when the two
# Poisson-inverse Gaussians gained an exact one; the family these tests use is
# pig_bare_distrib(), from helper-pigbare.R.

test_that("the probe finds the families whose expected information reads the data", {
  expect_true(expected_is_opg(pig_bare_distrib(), distributions7::expect_opg()))
  # the Poisson-inverse Gaussians compute theirs exactly and are not probed
  for (d in list(distributions7::pig1_distrib(), distributions7::pig2_distrib())) {
    expect_false(expected_is_opg(d, distributions7::expect_opg()), label = d@distrib_name)
  }
  # the location-scale families compute theirs by quadrature and are not probed
  for (d in list(distributions7::skewnormal1_distrib(),
                 distributions7::skewnormal2_distrib(),
                 distributions7::skewt_distrib(),
                 distributions7::pseudohuber_distrib())) {
    expect_false(expected_is_opg(d, distributions7::expect_opg()), label = d@distrib_name)
  }
  # a wrapper inherits it from its parent
  expect_true(expected_is_opg(
    distributions7::fixed(pig_bare_distrib(), sigma = 2), distributions7::expect_opg()))
  # a family that writes it out is not probed, and a truncated family, which
  # does not write it out either, returns a quadrature that does not move
  # with the response
  expect_false(expected_is_opg(distributions7::gamma1_distrib(), distributions7::expect_opg()))
  expect_false(expected_is_opg(distributions7::negbin2_distrib(), distributions7::expect_opg()))
  expect_false(expected_is_opg(
    distributions7::truncated(distributions7::gaussian1_distrib(), lower = -1),
    distributions7::expect_opg()))
  # an expectation asked for by name is an expectation
  expect_false(expected_is_opg(pig_bare_distrib(), distributions7::expect_bartlett()))
  # and the caller's random stream is left where it was
  set.seed(3); before <- .Random.seed
  expected_is_opg(pig_bare_distrib(), distributions7::expect_opg())
  expect_identical(.Random.seed, before)
})


test_that("the refusal names the family and the remedy, and only where it applies", {
  pig <- pig_bare_distrib()
  for (m in list(reml(information = expected_info()),
                 ml(information = expected_info()),
                 aic(information = expected_info()),
                 bic(information = expected_info()))) {
    expect_error(assert_criterion_information(pig, m),
                 "outer product of the scores")
  }
  expect_error(assert_criterion_information(pig, reml(information = expected_info())),
               "reml(information = observed_info())", fixed = TRUE)
  expect_error(assert_criterion_information(pig, aic(information = expected_info())),
               "aic(information = observed_info())", fixed = TRUE)
  expect_error(assert_criterion_information(pig, reml(information = expected_info())),
               "expected_info(approx = expect_bartlett())", fixed = TRUE)
  # not on the observed information, not on a family that writes it out, not
  # on an expectation asked for by name, and not where there is no criterion
  expect_null(assert_criterion_information(pig, reml(information = observed_info())))
  expect_null(assert_criterion_information(distributions7::gamma1_distrib(),
                                           reml(information = expected_info())))
  expect_null(assert_criterion_information(
    pig, reml(information = expected_info(approx = distributions7::expect_bartlett()))))
  expect_null(assert_criterion_information(pig, NULL))
})


test_that("statmod() rejects it at the start and the observed route is exact", {
  skip_on_cran()
  set.seed(11)
  n <- 400
  dd <- data.frame(x = runif(n))
  mu <- exp(0.8 + 0.5 * sin(2 * pi * dd$x))
  dd$y <- distributions7::distrib_rng(distributions7::pig1_distrib(), n,
                                      list(mu = mu, sigma = 0.5))
  fo <- y ~ s(x, bspline_smooth(k = 8))
  expect_error(statmod(fo, pig_bare_distrib(), dd,
                       outer_criterion = reml(information = expected_info())),
               "pig bare")
  fit <- statmod(fo, pig_bare_distrib(), dd,
                 outer_criterion = reml(information = observed_info()))
  ct <- statmod_certificate(fit)
  expect_identical(ct$curvature, "analytic")
  expect_identical(ct$state, "converged")
  # a model with nothing for the criterion to estimate is not refused: the
  # criterion is inert there and the check runs after it is dropped
  expect_no_error(statmod(y ~ x, pig_bare_distrib(), dd,
                          outer_criterion = reml(information = expected_info(),
                                                 marginal = "none")))
  # while by default it estimates sigma, and is refused like any other
  expect_error(statmod(y ~ x, pig_bare_distrib(), dd,
                       outer_criterion = reml(information = expected_info())), "pig bare")
  # and the shipped family, exact since distributions7 0.65.0, is fitted on
  # the expected route with an analytic certificate
  fe <- statmod(fo, distributions7::pig1_distrib(), dd,
                outer_criterion = reml(information = expected_info()))
  ce <- statmod_certificate(fe)
  expect_identical(ce$curvature, "analytic")
  expect_identical(ce$state, "converged")
})
