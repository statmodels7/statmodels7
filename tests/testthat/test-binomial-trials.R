# cbind(successes, failures), and the number of trials that goes with it.
#
# Aggregated binomial data are written cbind(y, n - y) ~ x in R. Before
# split_binomial_matrix() the matrix went down as a vector of twice the length
# and the run died on "Parameter dimension mismatch. All parameters should
# have length 1 or 400", a message about a length nobody had asked for.

set.seed(1)
n <- 200
dd <- data.frame(x = rnorm(n), g = factor(sample(letters[1:3], n, TRUE)))
dd$sz   <- sample(5:40, n, TRUE)
dd$y    <- rbinom(n, dd$sz, plogis(-0.3 + 0.8 * dd$x))
dd$fail <- dd$sz - dd$y

test_that("cbind on a binomial family agrees with glm", {
  rif <- glm(cbind(y, fail) ~ x + g, family = binomial, data = dd)
  s <- statmod(cbind(y, fail) ~ x + g, distributions7::binomial_distrib(),
               data = dd, inner_optimizer = iwls(tol = 1e-10))
  expect_equal(unname(coef(s)$mu), unname(coef(rif)), tolerance = 1e-8)
  expect_equal(as.numeric(logLik(s)), as.numeric(logLik(rif)),
               tolerance = 1e-8)
})

test_that("cbind and size = are the same model, written twice", {
  a <- statmod(cbind(y, fail) ~ x + g, distributions7::binomial_distrib(),
               data = dd, inner_optimizer = iwls(tol = 1e-10))
  b <- statmod(y ~ x + g, distributions7::binomial_distrib(size = dd$sz),
               data = dd, inner_optimizer = iwls(tol = 1e-10))
  expect_equal(coef(a)$mu, coef(b)$mu, tolerance = 1e-10)
  expect_identical(as.numeric(a@spec@distrib@size), as.numeric(dd$sz))
})

test_that("the trials come from the response, so they follow other rows", {
  # this is the whole reason cbind() is worth having: the row sums are
  # recomputed wherever the expression is evaluated, where a `size` given to
  # the family is a vector fixed when the fit was written
  s <- statmod(cbind(y, fail) ~ x, distributions7::binomial_distrib(), dd)
  nd <- dd[1:20, ]
  spec <- statmod_respec(s@spec, nd)
  expect_identical(spec@n_obs, 20L)
  expect_identical(as.numeric(spec@distrib@size), as.numeric(nd$sz))
  expect_identical(as.numeric(spec@response), as.numeric(nd$y))
})

test_that("a size that cannot follow the rows is refused, not recycled", {
  # the negative control is at the density: 200 sizes against 20 observations
  # return 200 log-densities and signal nothing, so the count of terms comes
  # from the family rather than from the data
  muto <- distributions7::distrib_pdf(
    distributions7::binomial_distrib(size = dd$sz), dd$y[1:20],
    list(mu = rep(0.4, 20)), log = TRUE)
  expect_length(muto, 200L)

  s <- statmod(y ~ x, distributions7::binomial_distrib(size = dd$sz), dd)
  expect_error(statmod_respec(s@spec, dd[1:20, ]),
               "200 numbers of trials and these are 20 rows")
  # the same fit at its own rows is untouched
  expect_identical(statmod_respec(s@spec, dd)@n_obs, 200L)
})

test_that("cbind needs a family that HAS a number of trials", {
  # bernoulli is one trial by construction: there is nothing to write the
  # successes and failures onto, and the message says which family to use
  expect_error(
    statmod(cbind(y, fail) ~ x, distributions7::bernoulli_distrib(), dd),
    "no\n  number of trials")
  expect_error(
    statmod(cbind(y, fail) ~ x, distributions7::bernoulli_distrib(), dd),
    "binomial_distrib")
  expect_error(
    statmod(cbind(y, fail) ~ x, distributions7::gaussian1_distrib(), dd),
    "no\n  number of trials")
})

test_that("a matrix that is not two columns is refused", {
  dd$z <- rnorm(n)
  expect_error(
    statmod(cbind(y, fail, z) ~ x, distributions7::binomial_distrib(), dd),
    "matrix with 3 columns")
})

test_that("the trials given twice must agree", {
  expect_error(
    statmod(cbind(y, fail) ~ x,
            distributions7::binomial_distrib(size = 99), dd),
    "given twice and the two disagree")
  # and where they do agree there is nothing to choose
  d2 <- dd
  d2$sz <- 10L
  d2$y <- rbinom(n, 10L, 0.4)
  d2$fail <- 10L - d2$y
  expect_no_error(
    statmod(cbind(y, fail) ~ x,
            distributions7::binomial_distrib(size = 10), d2))
})

test_that("a multivariate family keeps its matrix response", {
  # a matrix is the ordinary response there, one observation per row, and the
  # split must not touch it
  fam <- distributions7::mvgaussian1_distrib(2)
  Y <- matrix(rnorm(40), 20, 2)
  out <- statmodels7:::split_binomial_matrix(Y, fam)
  expect_identical(out$response, Y)
  expect_identical(out$distrib, fam)
})
