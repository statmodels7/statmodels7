# The information on which the effective degrees of freedom are read:
# 1. aic() and bic() read the count on their own information;
# 2. fit@edf reads the information of reml() or ml(), and the inner
#    information only where there is no outer criterion;
# 3. summary() reads the count on the information of its vcov_type;
# 4. a note when the observed information used for the count is indefinite.

edf_data <- function(n = 200, seed = 7) {
  set.seed(seed)
  d <- data.frame(x = runif(n))
  d$y <- sin(2 * pi * d$x) + rnorm(n, sd = 0.4)
  d
}

edf_fit <- function(oc, d = edf_data(), inner = auto_info()) {
  statmod(y ~ s(x, smoother = basis7::bspline_smooth(k = 12)),
          distributions7::gaussian1_distrib(), d, outer_criterion = oc,
          inner_information = inner)
}

# the count read directly on a stated information
edf_on <- function(fit, expected) {
  statmod_edf(fit@spec, fit@coefficients, statmod_design(fit@spec),
              fit@hyper, expected = expected)
}

smooth_edf <- function(e) e$edf[grepl("^s\\(", e$term)]

test_that("fit@edf reads the information of reml() and ml()", {
  skip_on_cran()
  for (oc in list(reml(), ml())) {
    fo <- edf_fit(oc)
    expect_equal(fo@edf$edf, edf_on(fo, FALSE)$edf, tolerance = 1e-12)
    expect_true(is(fo@methods$edf_information, "distributions7::ObservedInfo") ||
                  !info_expected(fo@methods$edf_information))
    fe <- edf_fit(if (oc@kind == "reml") reml(information = expected_info())
                  else ml(information = expected_info()))
    expect_equal(fe@edf$edf, edf_on(fe, TRUE)$edf, tolerance = 1e-12)
    expect_true(info_expected(fe@methods$edf_information))
    expect_gt(abs(smooth_edf(fo@edf) - smooth_edf(fe@edf)), 1e-6)
  }
})

test_that("the outer information wins over the inner one under reml()", {
  skip_on_cran()
  f <- edf_fit(reml(), inner = distributions7::expected_info())
  expect_equal(f@edf$edf, edf_on(f, FALSE)$edf, tolerance = 1e-12)
})

test_that("a fit without an outer criterion reads the inner information", {
  skip_on_cran()
  fe <- edf_fit(NULL, inner = distributions7::expected_info())
  expect_equal(fe@edf$edf, edf_on(fe, TRUE)$edf, tolerance = 1e-12)
  fo <- edf_fit(NULL, inner = distributions7::observed_info())
  expect_equal(fo@edf$edf, edf_on(fo, FALSE)$edf, tolerance = 1e-12)
  # auto_info() is settled against the family, as the inner fit settles it
  fa <- edf_fit(NULL)
  cfg <- inner_settings(inner_bind(iwls(), fa@spec@distrib, auto_info()))
  expect_equal(fa@edf$edf, edf_on(fa, cfg$expected)$edf, tolerance = 1e-12)
})

test_that("aic() and bic() keep the count on their own information", {
  skip_on_cran()
  for (oc in list(aic(information = expected_info()),
                  bic(information = expected_info()))) {
    f <- edf_fit(oc)
    expect_equal(f@edf$edf, edf_on(f, TRUE)$edf, tolerance = 1e-12)
    expect_true(f@methods$edf_priced)
    # the summary keeps it, whatever its vcov_type
    expect_identical(summary(f)@edf, f@edf)
    expect_identical(summary(f, vcov_type = bayesian_vcov(
      information = expected_info()))@edf, f@edf)
  }
  fa <- edf_fit(aic())
  expect_equal(fa@edf$edf, edf_on(fa, FALSE)$edf, tolerance = 1e-12)
  expect_identical(summary(fa)@edf, fa@edf)
})

test_that("summary() reports the count of fit@edf under all defaults", {
  skip_on_cran()
  for (oc in list(reml(), ml())) {
    f <- edf_fit(oc)
    expect_identical(summary(f)@edf, f@edf)
  }
  fe <- edf_fit(reml(information = expected_info()))
  # the standard errors default to the observed information, so the count
  # follows them and leaves the one stored in the fit
  s <- summary(fe)
  expect_equal(s@edf$edf, edf_on(fe, FALSE)$edf, tolerance = 1e-12)
  expect_equal(fe@edf$edf, edf_on(fe, TRUE)$edf, tolerance = 1e-12)
})

test_that("summary() recomputes the count on the information of vcov_type", {
  skip_on_cran()
  f <- edf_fit(reml())
  before <- f@edf
  s <- summary(f, vcov_type = bayesian_vcov(information = expected_info()))
  expect_equal(s@edf$edf, edf_on(f, TRUE)$edf, tolerance = 1e-12)
  expect_gt(abs(smooth_edf(s@edf) - smooth_edf(f@edf)), 1e-6)
  expect_identical(f@edf, before)
  # the other conventions read their information too
  s2 <- summary(f, vcov_type = frequentist_vcov(information = expected_info()))
  expect_equal(s2@edf$edf, s@edf$edf, tolerance = 1e-12)
  # and the count printed is the one of the summary
  expect_equal(sum(s@edf$edf), sum(edf_on(f, TRUE)$edf), tolerance = 1e-12)
})

test_that("a NULL information of vcov_type is resolved by the family", {
  skip_on_cran()
  # gaussian: observed, so nothing is recomputed
  f <- edf_fit(reml())
  expect_identical(summary(f, vcov_type = bayesian_vcov())@edf, f@edf)
  expect_false(info_expected(vcov_information(bayesian_vcov(),
                                              f@spec@distrib)))
})

test_that("edf_info_indefinite() separates indefinite from semidefinite", {
  H <- diag(c(4, 2, 1))
  expect_false(edf_info_indefinite(H))
  # semidefinite with an exactly flat direction is not indefinite
  v <- c(1, 1, 0) / sqrt(2)
  expect_false(edf_info_indefinite(diag(3) - tcrossprod(v) * 1))
  expect_true(edf_info_indefinite(diag(c(4, 2, -1e-3))))
  expect_false(edf_info_indefinite(matrix(NA_real_, 2, 2)))
})

test_that("the share of a direction of negative curvature is negative", {
  H <- diag(c(4, 2, -0.5))
  S <- diag(c(0, 1, 1))
  d <- diag(solve(H + S) %*% H)
  expect_lt(d[3], 0)
  expect_true(all(d[1:2] >= 0 & d[1:2] <= 1))
})

test_that("summary() notes an indefinite observed information", {
  skip_on_cran()
  set.seed(3)
  n <- 80
  d <- data.frame(x = runif(n))
  d$y <- sin(2 * pi * d$x) + 0.3 * rt(n, 2)
  f <- statmod(y ~ s(x, smoother = basis7::bspline_smooth(k = 12)),
               distributions7::student_t1_distrib(), d,
               outer_criterion = reml(information = observed_info()))
  H <- as.matrix(statmod_information_at(f@spec, f@coefficients,
                                        statmod_design(f@spec),
                                        expected = FALSE))
  expect_true(edf_info_indefinite(H))
  pat <- "observed information is indefinite"
  so <- summary(f, vcov_type = bayesian_vcov(information = observed_info()))
  expect_true(any(grepl(pat, so@notes)))
  expect_true(any(grepl("negative", so@notes)))
  # the default variance object reads the observed information of this family
  expect_true(any(grepl(pat, summary(f)@notes)))
  # on the expected information there is no note
  se2 <- summary(f, vcov_type = bayesian_vcov(information = expected_info()))
  expect_false(any(grepl(pat, se2@notes)))
  # a definite observed information carries none
  g <- edf_fit(reml())
  expect_false(any(grepl(pat, summary(g)@notes)))
})
