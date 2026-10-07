test_that("the default marginal criterion uses the expected information for a t", {
  # Giovanni, 2026-10-07: the observed information of a Student t is not
  # positive definite at an outlier, and the Laplace determinant built on it
  # failed (GAGurine with nu ~ Age: 91.5 s, not converged)
  expect_identical(reml()@hessian, "auto")
  expect_identical(ml()@hessian, "auto")
  t1 <- distributions7::student_t1_distrib()
  expect_identical(outer_resolve(reml(), t1)@hessian, "expected")
  expect_identical(outer_resolve(ml(), distributions7::skewt_distrib())@hessian,
                   "expected")
  expect_identical(outer_resolve(reml(), distributions7::cauchy_distrib())@hessian,
                   "expected")
  # through a wrapper
  expect_identical(outer_resolve(reml(), distributions7::fixed(t1, nu = 4))@hessian,
                   "expected")
  expect_identical(outer_resolve(reml(), distributions7::gaussian1_distrib())@hessian,
                   "observed")
  # a value named by the caller is kept
  expect_identical(outer_resolve(reml(hessian = "observed"), t1)@hessian,
                   "observed")
  # the fit records the information it used
  set.seed(5)
  dd <- data.frame(x = runif(80))
  dd$y <- 1 + dd$x + 0.5 * rt(80, df = 3)
  fit <- statmod(y ~ x, distrib = t1, data = dd)
  expect_identical(fit@methods$outer@hessian, "expected")
  fg <- statmod(y ~ x, distrib = distributions7::gaussian1_distrib(), data = dd)
  expect_identical(fg@methods$outer@hessian, "observed")
})
