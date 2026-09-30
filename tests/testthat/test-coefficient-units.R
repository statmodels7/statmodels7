test_that("a coefficient on the criterion is read in no particular units", {
  skip_on_cran()
  # The restricted likelihood of dist ~ speed has its maximum at the variance
  # lm() reports, whatever link or scale carries it.
  s2 <- summary(stats::lm(dist ~ speed, data = cars))$sigma^2

  # sigma^2 on the identity link: an estimate of 236.5 with a standard error
  # of 48.3, started at 651 where the criterion is not concave in it
  v <- statmod(dist ~ speed, distrib = gaussian2_distrib(
                 link_sigma2 = identity_link()), data = cars)
  expect_equal(unname(coef(v)$sigma2), s2, tolerance = 1e-6)
  expect_lt(nrow(v@history$outer), 20)
  cv <- statmod_certificate(v)
  expect_identical(cv$state, "converged")
  expect_length(cv$boundary, 0)

  # sigma on the identity link in two sets of units: the same fit, and the
  # same certificate
  cc <- cars
  cc$dist100 <- 100 * cc$dist
  a <- statmod(dist ~ speed, distrib = gaussian1_distrib(
                 link_sigma = identity_link()), data = cc)
  b <- statmod(dist100 ~ speed, distrib = gaussian1_distrib(
                 link_sigma = identity_link()), data = cc)
  expect_equal(unname(coef(b)$sigma) / 100, unname(coef(a)$sigma),
               tolerance = 1e-5)
  expect_lt(nrow(b@history$outer), 20)
  expect_identical(statmod_certificate(a)$state, "converged")
  expect_identical(statmod_certificate(b)$state, "converged")
  expect_length(statmod_certificate(b)$boundary, 0)
})
