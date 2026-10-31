sep_bern <- distributions7::bernoulli_distrib()

test_that("complete separation is found along the separating direction", {
  d <- data.frame(x = 1:10, y = as.integer(1:10 > 5))
  f <- suppressWarnings(statmod(y ~ x, sep_bern, d))
  s <- separation_diagnosis(f)
  expect_true(s$found)
  expect_false(s$descent)
  dir <- s$directions[[1]]
  expect_identical(dir$kind, "recession")
  expect_setequal(names(dir$coefficients), c("mu:(Intercept)", "mu:x"))
  # the direction separates: the predictor rises where y = 1
  eta <- dir$coefficients[["mu:(Intercept)"]] + dir$coefficients[["mu:x"]] * d$x
  expect_true(all(eta[d$y == 1] > 0) && all(eta[d$y == 0] < 0))
})

test_that("a level whose counts are all zero is found, and only that level", {
  set.seed(1)
  d <- data.frame(g = factor(rep(c("a", "b", "c", "d"), each = 15)),
                  x = rnorm(60))
  d$y <- rpois(60, exp(1 + 0.3 * d$x))
  d$y[d$g == "d"] <- 0L
  f <- statmod(y ~ g + x, distributions7::poisson_distrib(), d)
  s <- separation_diagnosis(f)
  expect_true(s$found)
  expect_named(s$directions[[1]]$coefficients, "mu:gd")
  expect_lt(s$directions[[1]]$coefficients[["mu:gd"]], 0)
  expect_identical(unname(s$directions[[1]]$n_moved), 15L)
  ct <- statmod_certificate(f)
  expect_true(ct$separation$found)
  expect_match(ct$separation$message, "mu:gd moves toward minus infinity")
  # an observation with a count of one in that level makes the estimate exist
  d$y[d$g == "d"][3] <- 1L
  f1 <- statmod(y ~ g + x, distributions7::poisson_distrib(), d)
  expect_false(separation_diagnosis(f1)$found)
})

test_that("a group without zeros is found in the zero-inflation equation", {
  set.seed(3)
  d <- data.frame(g = factor(rep(c("A", "B"), each = 80)), x = rnorm(160))
  d$y <- ifelse(runif(160) < 0.3, 0L, rpois(160, exp(1.2 + 0.3 * d$x)))
  d$y[d$g == "B" & d$y == 0] <- 1L
  f <- suppressWarnings(statmod(
    y ~ x | zi ~ g, distributions7::zero_inflated(distributions7::poisson_distrib()),
    d))
  s <- separation_diagnosis(f)
  expect_true(s$found)
  expect_true(any(vapply(s$directions, function(z)
    identical(names(z$coefficients), "zi:gB"), logical(1))))
  # zi:gB is estimated on the outer criterion, whose search decides where it
  # stops; the diagnosis must not depend on that point
  for (z in c(-6, -10, -12, -16, -25)) {
    f2 <- f
    cf <- f2@coefficients
    cf$zi[2] <- z
    f2@coefficients <- cf
    s2 <- separation_diagnosis(f2)
    expect_true(s2$found, label = sprintf("zi:gB at %g", z))
    expect_true(any(vapply(s2$directions, function(w)
      identical(names(w$coefficients), "zi:gB"), logical(1))),
      label = sprintf("zi:gB named at %g", z))
  }
})

test_that("separation beside a random effect is found in the fixed effect", {
  set.seed(4)
  d <- data.frame(id = factor(rep(1:20, each = 8)), x = rnorm(160),
                  z = rbinom(160, 1, 0.2))
  d$y <- rbinom(160, 1, plogis(-0.3 + 0.7 * d$x + rnorm(20, 0, 0.8)[d$id]))
  d$y[d$z == 1] <- 1L
  f <- suppressWarnings(statmod(y ~ x + z + random(~ 1 | id), sep_bern, d))
  s <- separation_diagnosis(f)
  expect_true(s$found)
  expect_named(s$directions[[1]]$coefficients, "mu:z")
})

test_that("ordinary fits, near separation and penalized terms are not reported", {
  set.seed(11)
  d <- data.frame(x = rnorm(200), z = rnorm(200))
  d$y <- rbinom(200, 1, plogis(-0.5 + d$x - 0.5 * d$z))
  s <- separation_diagnosis(statmod(y ~ x + z, sep_bern, d))
  expect_false(s$found)
  expect_false(s$descent)
  # one overlapping pair: large but finite coefficients
  d2 <- data.frame(x = 1:10, y = as.integer(c(0, 0, 0, 0, 1, 0, 1, 1, 1, 1)))
  expect_false(separation_diagnosis(statmod(y ~ x, sep_bern, d2))$found)
  # a group of ones under a random intercept has a finite mode
  set.seed(19)
  d3 <- data.frame(id = factor(rep(1:20, each = 8)), x = rnorm(160))
  d3$y <- rbinom(160, 1, plogis(-0.3 + 0.7 * d3$x + rnorm(20, 0, 0.8)[d3$id]))
  d3$y[d3$id == "3"] <- 1L
  s3 <- separation_diagnosis(statmod(y ~ x + random(~ 1 | id), sep_bern, d3))
  expect_false(s3$found)
  expect_false(s3$descent)
})

test_that("a recession makes the certificate's state boundary and names the coefficients", {
  d <- data.frame(x = 1:10, y = as.integer(1:10 > 5))
  f <- suppressWarnings(statmod(y ~ x, sep_bern, d))
  expect_identical(certificate_core(f, 1e-2, 2e-3, 8)$state, "converged")
  ct <- statmod_certificate(f)
  expect_true(ct$separation$found)
  expect_length(ct$separation$message, length(ct$separation$directions))
  expect_identical(ct$state, "boundary")
  expect_setequal(ct$boundary_coef, c("mu:(Intercept)", "mu:x"))
  expect_setequal(ct$boundary, c("mu/(Intercept)", "mu/x"))
  expect_length(ct$boundary_key, 0L)
  expect_true(any(grepl("the estimate does not exist", ct$reason)))
})

test_that("summary() leaves the separated coefficient without inference and says why", {
  set.seed(1)
  d <- data.frame(g = factor(rep(c("a", "b", "c", "d"), each = 15)),
                  x = rnorm(60))
  d$y <- rpois(60, exp(1 + 0.3 * d$x))
  d$y[d$g == "d"] <- 0L
  f <- statmod(y ~ g + x, distributions7::poisson_distrib(), d)
  s <- summary(f)
  tb <- s@tables$mu[[1L]]$table
  row <- tb[tb$name == "gd", ]
  expect_true(is.na(row$se) && is.na(row$lower) && is.na(row$p_value))
  expect_false(anyNA(tb[tb$name == "x", "se"]))
  expect_true(any(grepl("mu:gd carries no standard error", s@notes)))
  expect_false(any(grepl("estimated on the outer criterion", s@notes)))
  out <- utils::capture.output(print(s))
  expect_false(any(grepl("outer   not available", out)))
})

test_that("a direction of descent leaves the state alone", {
  f <- suppressWarnings(statmod(
    Sepal.Length ~ Species | sigma ~ Species | nu ~ Species,
    distributions7::student_t1_distrib(), iris,
    inner_optimizer = iwls(maxit = 1L), outer_criterion = reml(marginal = "none")))
  ct <- statmod_certificate(f)
  expect_true(ct$separation$descent)
  expect_identical(ct$state, certificate_core(f, 1e-2, 2e-3, 8)$state)
})
