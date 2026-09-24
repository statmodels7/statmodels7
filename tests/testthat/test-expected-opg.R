# A criterion on the expected information needs an expectation.
#
# Six shipped families do not write their expected information out, and under
# the default iwls(approx = "opg") what they return is the outer product of
# the scores at the data, which depends on the response. A criterion built on
# it is neither the Laplace approximation nor its Fisher variant and has no
# exact outer derivatives, so statmod() rejects it and names the observed
# route, which is exact on all six.

test_that("the probe finds the families whose expected information reads the data", {
  six <- list(distributions7::pig1_distrib(), distributions7::pig2_distrib(),
              distributions7::skewnormal1_distrib(),
              distributions7::skewnormal2_distrib(),
              distributions7::skewt_distrib(), distributions7::pseudohuber_distrib())
  for (d in six) expect_true(expected_is_opg(d, "opg"), label = d@distrib_name)
  # a wrapper inherits it from its parent
  expect_true(expected_is_opg(
    distributions7::fixed(distributions7::skewt_distrib(), nu = 6), "opg"))
  expect_true(expected_is_opg(
    distributions7::zero_inflated(distributions7::pig1_distrib()), "opg"))
  # a family that writes it out is not probed, and a truncated family, which
  # does not write it out either, returns a quadrature that does not move
  # with the response
  expect_false(expected_is_opg(distributions7::gamma1_distrib(), "opg"))
  expect_false(expected_is_opg(distributions7::negbin2_distrib(), "opg"))
  expect_false(expected_is_opg(
    distributions7::truncated(distributions7::gaussian1_distrib(), lower = -1),
    "opg"))
  # an expectation asked for by name is an expectation
  expect_false(expected_is_opg(distributions7::pig1_distrib(), "bartlett"))
  # and the caller's random stream is left where it was
  set.seed(3); before <- .Random.seed
  expected_is_opg(distributions7::skewt_distrib(), "opg")
  expect_identical(.Random.seed, before)
})


test_that("the refusal names the family and the remedy, and only where it applies", {
  pig <- distributions7::pig1_distrib()
  for (m in list(reml("expected"), ml("expected"), aic(hessian = "expected"),
                 bic(hessian = "expected"))) {
    expect_error(assert_criterion_information(pig, m, "opg"),
                 "outer product of the scores")
  }
  expect_error(assert_criterion_information(pig, reml("expected"), "opg"),
               "reml\\(hessian = \"observed\"\\)")
  expect_error(assert_criterion_information(pig, aic(hessian = "expected"), "opg"),
               "aic\\(hessian = \"observed\"\\)")
  # not on the observed information, not on a family that writes it out, not
  # on an expectation asked for by name, and not where there is no criterion
  expect_null(assert_criterion_information(pig, reml("observed"), "opg"))
  expect_null(assert_criterion_information(distributions7::gamma1_distrib(),
                                           reml("expected"), "opg"))
  expect_null(assert_criterion_information(pig, reml("expected"), "bartlett"))
  expect_null(assert_criterion_information(pig, NULL, "opg"))
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
  expect_error(statmod(fo, distributions7::pig1_distrib(), dd,
                       outer_criterion = reml("expected")),
               "poisson-inverse gaussian")
  fit <- statmod(fo, distributions7::pig1_distrib(), dd,
                 outer_criterion = reml("observed"))
  ct <- statmod_certificate(fit)
  expect_identical(ct$curvature, "analytic")
  expect_identical(ct$state, "converged")
  # a model with nothing for the criterion to estimate is not refused: the
  # criterion is inert there and the check runs after it is dropped
  expect_no_error(statmod(y ~ x, distributions7::pig1_distrib(), dd,
                          outer_criterion = reml("expected")))
})
