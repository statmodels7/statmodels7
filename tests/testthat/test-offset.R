# An offset written in the formula, which is how R has always written one.

off_data <- function(n = 1500, seed = 7) {
  set.seed(seed)
  d <- data.frame(x = runif(n), g = factor(sample(20, n, TRUE)),
                  pop = runif(n, 100, 5000))
  d$lp <- log(d$pop)
  d$y <- rpois(n, exp(-7 + 0.9 * d$x + d$lp))
  d
}


test_that("offset() in the formula is the offsets argument", {
  d <- off_data()
  a <- statmod(y ~ x + offset(lp), distributions7::poisson_distrib(), d,
               outer_criterion = NULL)
  b <- statmod(y ~ x, distributions7::poisson_distrib(), d,
               offsets = list(mu = d$lp), outer_criterion = NULL)
  expect_equal(a@coefficients$mu, b@coefficients$mu)
  expect_equal(as.numeric(stats::logLik(a)), as.numeric(stats::logLik(b)))

  # and it is not the model without one, which is what used to be fitted:
  # model.matrix() drops an offset term, so the equation lost it entirely
  c0 <- statmod(y ~ x, distributions7::poisson_distrib(), d,
                outer_criterion = NULL)
  expect_gt(abs(a@coefficients$mu[1] - c0@coefficients$mu[1]), 1)
})


test_that("the offset enters the predictor and not the design", {
  d <- off_data()
  a <- statmod(y ~ x + offset(lp), distributions7::poisson_distrib(), d,
               outer_criterion = NULL)
  c0 <- statmod(y ~ x, distributions7::poisson_distrib(), d,
                outer_criterion = NULL)
  # no column of its own: an offset is a coefficient known to be one
  expect_equal(statmod_design(a@spec)$mu$npar,
               statmod_design(c0@spec)$mu$npar)
  expect_length(a@spec@offsets$mu, nrow(d))
  expect_equal(a@spec@offsets$mu, d$lp)
})


test_that("the offset survives prediction on new data", {
  d <- off_data()
  a <- statmod(y ~ x + offset(lp), distributions7::poisson_distrib(), d,
               outer_criterion = NULL)
  nd <- d[1:25, ]
  expect_equal(predict(a, "mu")[1:25], predict(a, "mu", newdata = nd))

  # the EXPRESSION is re-evaluated, so a different population gives a
  # different rate: carrying the offset as numbers could not do this, which
  # is why prediction used to return the predictor of a model without one
  nd2 <- nd
  nd2$lp <- nd2$lp + log(2)
  expect_equal(predict(a, "mu", newdata = nd2) / predict(a, "mu", newdata = nd),
               rep(2, nrow(nd)), tolerance = 1e-10)
})


test_that("several offsets are summed, and any equation may carry one", {
  d <- off_data()
  d$o2 <- rnorm(nrow(d), 0, 0.01)
  a <- statmod(y ~ x + offset(lp) + offset(o2),
               distributions7::poisson_distrib(), d, outer_criterion = NULL)
  b <- statmod(y ~ x, distributions7::poisson_distrib(), d,
               offsets = list(mu = d$lp + d$o2), outer_criterion = NULL)
  expect_equal(a@coefficients$mu, b@coefficients$mu)

  d$s <- rnorm(nrow(d))
  g <- statmod(y ~ x + offset(lp) | sigma ~ offset(s),
               distributions7::gaussian1_distrib(), d, outer_criterion = NULL)
  expect_equal(g@spec@offsets$sigma, d$s)
  expect_equal(g@spec@offsets$mu, d$lp)
})


test_that("an equation that is only an offset keeps its intercept", {
  d <- off_data()
  a <- statmod(y ~ offset(lp), distributions7::poisson_distrib(), d,
               outer_criterion = NULL)
  expect_equal(statmod_design(a@spec)$mu$npar, 1L)
  expect_equal(a@spec@offsets$mu, d$lp)
})


test_that("an offset that cannot be used says so", {
  d <- off_data()
  expect_error(statmod(y ~ x + offset(nope), distributions7::poisson_distrib(),
                       d, outer_criterion = NULL),
               "could not")
  # a factor is not an offset
  expect_error(statmod(y ~ x + offset(g), distributions7::poisson_distrib(),
                       d, outer_criterion = NULL),
               "must be numeric")
  # and neither is a column with a hole in it: an offset enters the predictor
  # as it stands, so there is no fitted answer at a missing value
  d$bad <- d$lp
  d$bad[3] <- NA
  expect_error(statmod(y ~ x + offset(bad), distributions7::poisson_distrib(),
                       d, outer_criterion = NULL),
               "finite")
})


test_that("an offset buried inside another term is refused", {
  d <- off_data()
  d$z <- runif(nrow(d))
  # every one of these used to FIT, with the term's own model.matrix dropping
  # the offset: the block had the columns of the model without it and the
  # intercept came back 566 times out
  expect_error(statmod(y ~ x + ridge(~ z + offset(lp)),
                       distributions7::poisson_distrib(), d,
                       outer_criterion = NULL),
               "inside another term")
  expect_error(statmod(y ~ x + random(~ 1 + offset(lp) | g),
                       distributions7::poisson_distrib(), d,
                       outer_criterion = NULL),
               "inside another term")
  expect_error(statmod(y ~ linpar(~ x + offset(lp)),
                       distributions7::poisson_distrib(), d,
                       outer_criterion = NULL),
               "inside another term")
  # the message says where it belongs
  expect_error(statmod(y ~ x + ridge(~ z + offset(lp)),
                       distributions7::poisson_distrib(), d,
                       outer_criterion = NULL),
               "term of the equation")

  # and the spelling that works is not caught
  expect_no_error(statmod(y ~ x + offset(lp) + ridge(~ z),
                          distributions7::poisson_distrib(), d,
                          outer_criterion = NULL))
})


test_that("reject_nested_offsets reads calls and not names", {
  # a column called `offset` used as a covariate is not an offset call
  expect_silent(reject_nested_offsets(~ x + offset, "mu"))
  expect_silent(reject_nested_offsets(~ x + z, "mu"))
  # buried at any depth, however spelled
  expect_error(reject_nested_offsets(~ x + f(g(offset(a))), "mu"),
               "inside another term")
  expect_error(reject_nested_offsets(~ ridge(~ z + stats::offset(a)), "mu"),
               "inside another term")
})


test_that("split_offsets takes only a top-level additive term", {
  f <- split_offsets(~ x + offset(a) + z)
  expect_equal(deparse(f$formula[[2]]), "x + z")
  expect_length(f$offsets, 1L)
  expect_equal(deparse(f$offsets[[1]]), "a")

  # stats::offset() is the same term written out
  f2 <- split_offsets(~ x + stats::offset(a))
  expect_length(f2$offsets, 1L)

  # nothing to take
  f3 <- split_offsets(~ x + z)
  expect_length(f3$offsets, 0L)
  expect_equal(deparse(f3$formula[[2]]), "x + z")
})
