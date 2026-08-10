# The specification: equations to built terms to a design.

set.seed(1)
dd <- data.frame(
  y = rnorm(40), x = runif(40), z = runif(40),
  g = factor(rep(letters[1:4], 10))
)
dd$R <- matrix(rnorm(80), 40, 2)

test_that("every parameter gets terms, in the family's order", {
  spec <- statmod_spec(y ~ x + g | sigma ~ z,
                       distributions7::gaussian1_distrib(), dd)
  expect_named(spec@terms, c("mu", "sigma"))
  expect_identical(spec@n_obs, 40L)
  expect_identical(spec@response, dd$y)
  expect_true(all(spec@weights == 1))
  d <- statmod_design(spec)
  # mu carries an intercept, x and the three contrasts of g
  expect_identical(d$mu$npar, 5L)
  expect_identical(d$sigma$npar, 2L)
  expect_identical(nrow(d$mu$X), 40L)
})

test_that("a parameter with no equation gets an intercept only", {
  spec <- statmod_spec(y ~ x, distributions7::gaussian1_distrib(), dd)
  d <- statmod_design(spec)
  expect_identical(d$sigma$npar, 1L)
  expect_identical(d$sigma$coef_names, "(Intercept)")
})

test_that("a penalized term is built and carries its penalty", {
  spec <- statmod_spec(y ~ x + ridge(R), distributions7::gaussian1_distrib(), dd)
  tms <- spec@terms$mu
  expect_length(tms, 2L)
  expect_true("linpar" %in% names(tms))
  pen <- lapply(tms, modelterms7::term_penalty)
  expect_true(any(!vapply(pen, is.null, logical(1))))
  d <- statmod_design(spec)
  # the blocks say which columns each term owns, which is what the penalty
  # assembly needs
  expect_identical(sum(lengths(d$mu$blocks)), d$mu$npar)
})

test_that("the design's blocks and coefficient names line up", {
  spec <- statmod_spec(y ~ x + g + ridge(R),
                       distributions7::gaussian1_distrib(), dd)
  d <- statmod_design(spec)
  expect_length(d$mu$coef_names, d$mu$npar)
  expect_identical(ncol(d$mu$X), d$mu$npar)
  expect_identical(sort(unlist(d$mu$blocks, use.names = FALSE)),
                   seq_len(d$mu$npar))
})

test_that("a three-parameter family is served", {
  spec <- statmod_spec(y ~ x | sigma ~ z | nu ~ 1,
                       distributions7::student_t1_distrib(), dd)
  expect_named(spec@terms, c("mu", "sigma", "nu"))
  d <- statmod_design(spec)
  expect_identical(vapply(d, function(z) z$npar, integer(1)),
                   c(mu = 2L, sigma = 2L, nu = 1L))
})

test_that("prior weights are taken as given and not normalized", {
  w <- runif(40, 0.5, 2)
  spec <- statmod_spec(y ~ x, distributions7::gaussian1_distrib(), dd,
                       weights = w)
  expect_identical(spec@weights, w)
  # the sum is whatever it is; normalizing would turn the log-likelihood into
  # a mean and rescale every standard error
  expect_false(isTRUE(all.equal(sum(spec@weights), 1)))
})

test_that("what the specification refuses", {
  expect_error(statmod_spec(y ~ x, distributions7::gaussian1_distrib(),
                            list(y = 1)), "must be a data frame")
  expect_error(statmod_spec(y ~ x, "gaussian", dd), "distributions7")
  expect_error(statmod_spec(y ~ x, distributions7::gaussian1_distrib(), dd,
                            weights = 1:3), "length 3")
  expect_error(statmod_spec(y ~ x, distributions7::gaussian1_distrib(), dd,
                            weights = rep(-1, 40)), "non-negative")
  expect_error(statmod_spec(y ~ x, distributions7::gaussian1_distrib(), dd,
                            offsets = list(wrong = 1)), "not a parameter")
})

test_that("an offset is recycled to the sample size", {
  spec <- statmod_spec(y ~ x, distributions7::gaussian1_distrib(), dd,
                       offsets = list(mu = 2))
  expect_length(spec@offsets$mu, 40L)
  expect_true(all(spec@offsets$mu == 2))
  expect_null(spec@offsets$sigma)
})

test_that("our terms win over an attached package's", {
  # the interpreter runs with modelterms7 in front, so a name another package
  # exports cannot change what a formula means
  shim <- statmodels7:::terms_first(globalenv())
  local({
    s <- function(...) stop("this must never be called")
    spec <- statmod_spec(y ~ s(x, k = 5), distributions7::gaussian1_distrib(),
                         dd)
    expect_true(any(vapply(spec@terms$mu,
                           function(tm) S7::S7_inherits(tm, modelterms7::SmoothTerm),
                           logical(1))))
  })
})
