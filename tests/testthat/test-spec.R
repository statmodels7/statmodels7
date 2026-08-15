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

test_that("a term whose block moves with its coefficients is accepted", {
  # it was rejected until the alternation learned to refresh a block; what
  # the design now carries is the list of terms it has to refresh
  for (call in list(quote(seg(x)), quote(jump(x)), quote(jseg(x)),
                    quote(nl(~ a * exp(b * x), params = c("a", "b"))))) {
    fm <- stats::as.formula(paste("y ~", deparse(call)))
    spec <- expect_no_error(
      statmod_spec(fm, distributions7::gaussian1_distrib(), dd))
    expect_length(statmod_refreshable(spec), 1L)
    expect_length(attr(statmod_design(spec), "refresh"), 1L)
  }
  # and an ordinary model carries none, so it reaches exactly the same
  # arithmetic as before
  plain <- statmod_spec(y ~ x + s(z), distributions7::gaussian1_distrib(), dd)
  expect_length(statmod_refreshable(plain), 0L)
  expect_null(attr(statmod_design(plain), "refresh"))
})

test_that("the refreshable terms are found in every equation", {
  spec <- statmod_spec(y ~ seg(x) | sigma ~ jump(z),
                       distributions7::gaussian1_distrib(), dd)
  rf <- statmod_refreshable(spec)
  expect_length(rf, 2L)
  expect_identical(vapply(rf, function(r) r$param, character(1)),
                   c("mu", "sigma"))
  expect_identical(vapply(rf, function(r) r$term, character(1)),
                   c("seg(x)", "jump(z)"))
})

test_that("a structural term is routed by the shape it implements", {
  tt <- data.frame(y = rnorm(60), t = 1:60)
  # both shapes are fitted, and which one a term is is read off the methods
  # it registers rather than from a list of class names
  spec <- statmod_spec(y ~ gas(1, 1, time = t) - 1,
                       distributions7::gaussian1_distrib(), tt)
  expect_identical(statmod_structural(spec)[[1L]]$kind, "filter")
  s2 <- statmod_spec(y ~ regime(2, time = t) - 1,
                     distributions7::gaussian1_distrib(), tt)
  expect_identical(statmod_structural(s2)[[1L]]$kind, "loglik")
})

test_that("the message names the term and its parameter", {
  # a structural term implementing NEITHER shape of the contract is what
  # remains unfittable, and the error says which term and which equation
  Hollow <- S7::new_class("Hollow", parent = modelterms7::structural_term)
  build <- modelterms7::term_build
  S7::method(build, Hollow) <- function(term, data, ...) term
  tt <- data.frame(y = rnorm(60), z = runif(60))
  hollow <- function() Hollow(label = "hollow")
  err <- tryCatch(statmod_spec(y ~ z | sigma ~ hollow(),
                               distributions7::gaussian1_distrib(), tt),
                  error = conditionMessage)
  expect_match(err, "'hollow()'", fixed = TRUE)
  expect_match(err, "'sigma'", fixed = TRUE)
  expect_match(err, "neither shape of the contract", fixed = TRUE)
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

test_that("linpar_options reaches the implicit parametric block", {
  # the IMPLICIT linpar -- the bare covariates collapsed into one term -- is
  # never written by the caller, so this is the only place its arguments can
  # come from
  set.seed(61)
  n <- 400L
  m <- 60L
  d <- data.frame(g = factor(sample.int(m, n, TRUE)), z = stats::rnorm(n))
  d$y <- stats::rnorm(n)

  dense <- statmod_spec(y ~ 0 + g + z, distributions7::gaussian1_distrib(),
                        d, NULL, NULL)
  sp <- statmod_spec(y ~ 0 + g + z, distributions7::gaussian1_distrib(),
                     d, NULL, NULL, linpar = linpar_options(sparse = TRUE))
  Xd <- statmod_design(dense)$mu$X
  Xs <- statmod_design(sp)$mu$X
  expect_true(is.matrix(Xd))
  expect_true(methods::is(Xs, "sparseMatrix"))
  # the same design, only stored differently
  expect_equal(unname(as.matrix(Xs)), unname(as.matrix(Xd)))
  expect_lt(as.numeric(utils::object.size(Xs)),
            as.numeric(utils::object.size(Xd)) / 3)

  # and the SPECIFICATION carries it, which is what lets a rebuild reproduce
  # the storage rather than quietly densifying
  expect_true(isTRUE(sp@linpar$sparse))
  fold <- statmod_spec(sp@formula, sp@distrib, d[1:200, , drop = FALSE],
                       NULL, NULL, linpar = sp@linpar)
  expect_true(methods::is(statmod_design(fold)$mu$X, "sparseMatrix"))

  # the contrasts reach it too, and are a different coding of the same model
  ct <- statmod_spec(y ~ g, distributions7::gaussian1_distrib(), d, NULL, NULL,
                     linpar = linpar_options(contrasts = list(g = "contr.sum")))
  pl <- statmod_spec(y ~ g, distributions7::gaussian1_distrib(), d, NULL, NULL)
  Xc <- statmod_design(ct)$mu$X
  Xp <- statmod_design(pl)$mu$X
  expect_identical(ncol(Xc), ncol(Xp))
  expect_false(isTRUE(all.equal(unname(Xc), unname(Xp))))

  expect_error(linpar_options(sparse = 1), "TRUE or FALSE")
  expect_error(linpar_options(contrasts = "contr.sum"), "named list")
})

test_that("a sparse parametric block fits to the same answer", {
  set.seed(62)
  n <- 500L
  m <- 80L
  d <- data.frame(g = factor(sample.int(m, n, TRUE)), z = stats::rnorm(n))
  b <- stats::rnorm(m, sd = 0.6)
  d$y <- b[as.integer(d$g)] + 0.8 * d$z + stats::rnorm(n, sd = 0.5)

  a <- statmod(y ~ 0 + g + z, distributions7::gaussian1_distrib(), d)
  s <- statmod(y ~ 0 + g + z, distributions7::gaussian1_distrib(), d,
               linpar_control = linpar_options(sparse = TRUE))
  # the storage is a storage: the fit is the same model, to the digit
  expect_equal(a@loglik, s@loglik)
  expect_equal(unname(unlist(a@coefficients)), unname(unlist(s@coefficients)))
  expect_true(methods::is(statmod_design(s@spec)$mu$X, "sparseMatrix"))
  expect_error(statmod(y ~ z, distributions7::gaussian1_distrib(), d,
                       linpar_control = "sparse"),
               "linpar_options()", fixed = TRUE)
})
