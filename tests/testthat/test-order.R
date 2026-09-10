## The layer reads the order of differentiability distributions7 declares and
## refuses by name where a consumer cannot run.  Lot 2 of piano_nonsmooth.txt.

lap_data <- function(n = 200L, seed = 1L) {
  set.seed(seed)
  m <- 20L; ni <- n %/% m
  id <- factor(rep(seq_len(m), each = ni))
  x  <- stats::rnorm(m * ni)
  b  <- stats::rnorm(m, 0, 0.8)
  u  <- stats::runif(m * ni) - 0.5
  data.frame(y = 1 + 0.5 * x + b[as.integer(id)] -
               0.7 * sign(u) * log(1 - 2 * abs(u)),
             x = x, id = id)
}

test_that("a criterion refuses a kinked family by name, with the remedy", {
  dd <- lap_data()
  msg <- tryCatch(statmod(y ~ x + random(~ 1 | id),
                          distrib = distributions7::laplace_distrib(),
                          data = dd),
                  error = conditionMessage)
  expect_type(msg, "character")
  # the four things the message must name
  expect_match(msg, "'laplace'",  fixed = TRUE)   # the family
  expect_match(msg, "'mu'",       fixed = TRUE)   # the parameter
  expect_match(msg, "2 derivatives", fixed = TRUE) # the order it needs
  expect_match(msg, "0 of them",  fixed = TRUE)   # the order there is
  expect_match(msg, "outer_criterion = NULL", fixed = TRUE)  # the remedy
  # and NOT the message it replaces, which spoke of something else
  expect_false(grepl("the inner fit did not converge there", msg, fixed = TRUE))
})

test_that("every criterion gives the same refusal, since all five read the curvature", {
  dd <- lap_data()
  crits <- list(reml(hessian = "observed"), reml(hessian = "expected"),
                ml(), aic(), bic())
  seen <- vapply(crits, function(cr) {
    tryCatch({
      statmod(y ~ x + random(~ 1 | id),
              distrib = distributions7::laplace_distrib(),
              data = dd, outer_criterion = cr)
      NA_character_
    }, error = conditionMessage)
  }, character(1))
  expect_false(anyNA(seen))
  for (s in seen) {
    expect_match(s, "'laplace'", fixed = TRUE)
    expect_match(s, "'mu'", fixed = TRUE)
    expect_match(s, "derivatives of the log-density", fixed = TRUE)
  }
  # they differ only in which criterion they name, and the name appears
  # TWICE in the message, so every occurrence has to go before comparing
  stripped <- gsub("the [A-Z]+ criterion", "<crit>", seen)
  expect_length(unique(stripped), 1L)
})

test_that("a fit with no criterion is NOT refused: the point estimate is available", {
  dd <- lap_data()
  f <- statmod(y ~ x + random(~ 1 | id),
               distrib = distributions7::laplace_distrib(),
               data = dd, outer_criterion = NULL)
  expect_true(is.finite(as.numeric(logLik(f))))
  # what the refusal protects is the CRITERION, which reads a curvature that
  # is not there; the coefficients are a different question and 1c of
  # piano_nonsmooth.txt measures them good to ~1e-4 of the exact LAD fit.
})

test_that("a smooth family is untouched", {
  dd <- lap_data()
  expect_silent(order_available(distributions7::gaussian1_distrib(), 5L))
  expect_true(order_available(distributions7::gaussian1_distrib(), 5L))
  expect_null(order_shortfall(distributions7::gaussian1_distrib(), 4L, "x"))
  f <- statmod(y ~ x + random(~ 1 | id),
               distrib = distributions7::gaussian1_distrib(), data = dd)
  expect_true(is.finite(as.numeric(logLik(f))))
})

test_that("INJECTION: the order gates each consumer at its own rung", {
  # order_available() is the one predicate every consumer reads, so the
  # injection is done on it: a family reporting m = 2 must satisfy the
  # criterion (which needs 2) and not the exact outer gradient (which needs
  # 3), and one reporting m = 1 must satisfy neither.
  d <- distributions7::gaussian1_distrib()
  fake_order <- function(m) function(distrib, ...) c(mu = m, sigma = Inf)

  for (m in c(1, 2, 3, 4)) {
    testthat::local_mocked_bindings(params_order = fake_order(m),
                                    .package = "distributions7")
    expect_identical(order_available(d, 2L), m >= 2, info = paste("m =", m))
    expect_identical(order_available(d, 3L), m >= 3, info = paste("m =", m))
    expect_identical(order_available(d, 4L), m >= 4, info = paste("m =", m))
    # the criterion refuses at m < 2 and not above
    why <- order_shortfall(d, 2L, "the REML criterion")
    expect_identical(is.null(why), m >= 2, info = paste("m =", m))
  }
})

test_that("INJECTION: at m = 2 the criterion runs and the exact gradient does not", {
  dd <- lap_data()
  # the control: with the real family both are available
  f0 <- statmod(y ~ x + random(~ 1 | id),
                distrib = distributions7::gaussian1_distrib(), data = dd)
  expect_true(is.finite(as.numeric(logLik(f0))))

  testthat::local_mocked_bindings(
    params_order = function(distrib, ...) c(mu = 2, sigma = Inf),
    .package = "distributions7")
  # IWLS and the criterion need 2 and are not refused
  expect_true(order_available(distributions7::gaussian1_distrib(), 2L))
  expect_null(order_shortfall(distributions7::gaussian1_distrib(), 2L, "the REML criterion"))
  # the exact outer gradient needs 3 and is
  expect_false(order_available(distributions7::gaussian1_distrib(), 3L))
  f <- statmod(y ~ x + random(~ 1 | id),
               distrib = distributions7::gaussian1_distrib(), data = dd)
  expect_true(is.finite(as.numeric(logLik(f))))
})
