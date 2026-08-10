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

test_that("a design carrying a tensor smooth and an intercept has full rank", {
  # te() contains the constant and so does its penalty's null space, so before
  # the centering constraint this design was rank deficient by exactly one --
  # 25 of 26 columns, smallest singular value 0, condition 8.0e15 -- and chol()
  # accepted the penalized information anyway, so vcov(), confint() and the
  # outer criterion returned numbers computed on a singular matrix
  set.seed(24)
  dt <- data.frame(a = runif(300, -1, 1), b = runif(300, -1, 1))
  dt$y <- dt$a^2 + dt$b + stats::rnorm(300, sd = 0.3)
  spec <- statmod_spec(y ~ te(a, b, k = 5),
                       distributions7::gaussian1_distrib(), dt)
  des <- statmod_design(spec)
  X <- des$mu$X
  expect_identical(qr(X)$rank, ncol(X))
  sv <- svd(X)$d
  expect_lt(sv[1] / sv[length(sv)], 1e6)

  # and the penalized information is definite by a margin, not by the luck of
  # rounding: chol() succeeding is not the assertion, the eigenvalue is
  S <- statmod_penalty_at(spec, lapply(des, function(d) numeric(d$npar)),
                          statmod_hyper_start(spec), design = des,
                          what = "hessian")
  mu <- seq_len(des$mu$npar)
  ev <- eigen(crossprod(X) + S[mu, mu], symmetric = TRUE,
              only.values = TRUE)$values
  expect_gt(min(ev), sqrt(.Machine$double.eps) * max(ev))
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

# A term the fitting scheme does not cover is rejected at specification time,
# where the term can be named, rather than several frames down in the design.

test_that("a term whose block moves with its coefficients is rejected", {
  for (call in list(quote(seg(x)), quote(jump(x)), quote(jseg(x)))) {
    fm <- stats::as.formula(paste("y ~", deparse(call)))
    expect_error(
      statmod_spec(fm, distributions7::gaussian1_distrib(), dd),
      "block depends on its own coefficients")
  }
})

test_that("nl() is rejected for the same reason", {
  expect_error(
    statmod_spec(y ~ nl(~ a * exp(b * x), params = c(a = 1, b = 0.1)),
                 distributions7::gaussian1_distrib(), dd),
    "block depends on its own coefficients")
})

test_that("a structural term is rejected, and says which shape it is", {
  tt <- data.frame(y = rnorm(60), t = 1:60)
  for (call in list(quote(gas(1, 1)), quote(regime(2)))) {
    fm <- stats::as.formula(paste("y ~", deparse(call)))
    expect_error(statmod_spec(fm, distributions7::gaussian1_distrib(), tt),
                 "structural term rewrites the likelihood")
  }
})

test_that("the message names the term and its parameter", {
  err <- tryCatch(statmod_spec(y ~ x | sigma ~ seg(z),
                               distributions7::gaussian1_distrib(), dd),
                  error = conditionMessage)
  expect_match(err, "'seg(z)'", fixed = TRUE)
  expect_match(err, "'sigma'", fixed = TRUE)
})

test_that("every offending equation is reported, not just the first", {
  err <- tryCatch(statmod_spec(y ~ seg(x) | sigma ~ jump(z),
                               distributions7::gaussian1_distrib(), dd),
                  error = conditionMessage)
  expect_match(err, "seg(x)", fixed = TRUE)
  expect_match(err, "jump(z)", fixed = TRUE)
})

test_that("the terms the scheme does cover are not rejected", {
  # the guard reads a property off the term, so widening it would show here
  for (fm in list(y ~ x, y ~ s(x), y ~ te(x, z), y ~ s(x, by = g),
                  y ~ x + random(~ 1 | g), y ~ ridge(R), y ~ lasso(R),
                  y ~ x | sigma ~ z)) {
    expect_no_error(statmod_spec(fm, distributions7::gaussian1_distrib(), dd))
  }
})

test_that("the predicate reads the method's owner, not a list of classes", {
  built <- function(tm) modelterms7::term_build(tm, dd)
  expect_true(statmodels7:::refreshes_own_block(built(seg(x))))
  expect_true(statmodels7:::refreshes_own_block(
    built(nl(~ a * exp(b * x), params = c(a = 1, b = 0.1)))))
  # these inherit term_refresh's identity method on model_term
  expect_false(statmodels7:::refreshes_own_block(built(linpar(~x))))
  expect_false(statmodels7:::refreshes_own_block(built(s(x))))
  expect_false(statmodels7:::refreshes_own_block(built(ridge(R))))
  expect_false(statmodels7:::refreshes_own_block(built(random(~ 1 | g))))
})
