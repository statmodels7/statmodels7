test_that("auto_info() uses the expected information for a t", {
  # the default is the observed information (Giovanni, 2026-10-07, evening);
  # "auto" settles it against the family: the observed information of a
  # Student t is not positive definite at an outlier
  expect_identical(info_label(reml()@hessian), "observed")
  expect_identical(info_label(ml()@hessian), "observed")
  t1 <- distributions7::student_t1_distrib()
  auto <- reml(information = auto_info())
  expect_identical(info_label(outer_resolve(auto, t1)@hessian), "expected")
  expect_identical(info_label(outer_resolve(ml(information = auto_info()),
                                 distributions7::skewt_distrib())@hessian),
                   "expected")
  expect_identical(info_label(outer_resolve(auto, distributions7::cauchy_distrib())@hessian),
                   "expected")
  # through a wrapper
  expect_identical(info_label(outer_resolve(auto, distributions7::fixed(t1, nu = 4))@hessian),
                   "expected")
  expect_identical(info_label(outer_resolve(auto, distributions7::gaussian1_distrib())@hessian),
                   "observed")
  # a value named by the caller is kept
  expect_identical(info_label(outer_resolve(reml(information = observed_info()), t1)@hessian),
                   "observed")
  # the fit records the information it used
  set.seed(5)
  dd <- data.frame(x = runif(80))
  dd$y <- 1 + dd$x + 0.5 * rt(80, df = 3)
  fit <- statmod(y ~ x, distrib = t1, data = dd,
                 outer_criterion = reml(information = auto_info()))
  expect_identical(info_label(fit@methods$outer@hessian), "expected")
  f0 <- statmod(y ~ x, distrib = t1, data = dd)
  expect_identical(info_label(f0@methods$outer@hessian), "observed")
  fg <- statmod(y ~ x, distrib = distributions7::gaussian1_distrib(), data = dd)
  expect_identical(info_label(fg@methods$outer@hessian), "observed")
})
