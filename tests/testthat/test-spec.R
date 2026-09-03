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

  expect_error(linpar_options(sparse = 1), "TRUE, FALSE, or NULL")
  # NULL is the default and reaches linpar(), which settles the storage from
  # the size of the design rather than being told
  expect_null(linpar_options()$sparse)
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

test_that("a smooth per level fits the same model in either storage", {
  # the storage is a storage, and the penalty that avoids assembling
  # I (x) P is the same penalty: everything a reader looks at must agree
  skip_on_cran()
  set.seed(71)
  n <- 800L
  m <- 25L
  d <- data.frame(x = stats::runif(n), g = factor(sample.int(m, n, TRUE)))
  b <- stats::rnorm(m, sd = 0.5)
  d$y <- b[as.integer(d$g)] + sin(3 * d$x) + stats::rnorm(n, sd = 0.4)

  # the dense side is asked for EXPLICITLY: left NULL the storage is settled
  # from the size of the block, and 800 rows by 25 levels by a basis of eight
  # is 160000 cells, past the threshold, so the default settles sparse here
  a <- statmod(y ~ s(x, k = 8, by = g, sparse = FALSE),
               distributions7::gaussian1_distrib(), d)
  s <- statmod(y ~ s(x, k = 8, by = g, sparse = TRUE),
               distributions7::gaussian1_distrib(), d)

  expect_true(methods::is(statmod_design(s@spec)$mu$X, "sparseMatrix"))
  expect_true(is.matrix(statmod_design(a@spec)$mu$X))
  expect_equal(s@loglik, a@loglik)
  expect_equal(s@hyper$mu[[1L]][[1L]], a@hyper$mu[[1L]][[1L]])
  expect_equal(s@edf$edf, a@edf$edf)
  expect_equal(unname(unlist(s@coefficients)), unname(unlist(a@coefficients)))
  expect_equal(sqrt(diag(as.matrix(vcov(s)))),
               sqrt(diag(as.matrix(vcov(a)))))
  expect_equal(unlist(predict(s)), unlist(predict(a)))

  # the penalty is sparse and was NEVER assembled, yet its rank is the rank
  # the assembled matrix has -- which is the claim the shortcut rests on
  # the smooth, not the intercept's linpar, which comes first
  tm <- Filter(function(z) S7::S7_inherits(z, modelterms7::SmoothTerm),
               s@spec@terms$mu)[[1L]]
  pen <- modelterms7::term_penalty(tm)
  npar <- modelterms7::term_npar(tm)
  S <- penalties7::penalty_hessian(pen, numeric(npar), list(lambda = 1))
  expect_true(methods::is(S, "sparseMatrix"))
  ev <- eigen(as.matrix(S), symmetric = TRUE, only.values = TRUE)$values
  expect_identical(penalties7::penalty_rank(pen),
                   sum(ev > 1e-10 * max(ev)))
  # m times the block's, the eigenvalues of I (x) P being P's repeated
  expect_identical(penalties7::penalty_rank(pen) %% m, 0L)
})


# A covariance label is recorded and not yet fittable ----------------------

test_that("a labelled random effect becomes one covariance class", {
  # the terms carrying a label are collected across the equations into one
  # class, whose penalty covers their stacked columns
  spec <- statmod_spec(y ~ random(~ 1 | u | g) | sigma ~ random(~ 1 | u | g),
                       distributions7::gaussian1_distrib(), dd)
  des <- statmod_design(spec)
  us <- statmod_penalized(spec, des)
  expect_length(us, 1L)
  u <- us[[1L]]
  expect_identical(u$key, "u | g")
  expect_identical(u$params, c("mu", "sigma"))
  expect_identical(u$class$dim, 2L)
  expect_identical(u$class$m, nlevels(dd$g))
  # the class's index reaches BOTH equations
  expect_true(any(u$index <= des$mu$npar))
  expect_true(any(u$index > des$mu$npar))
  expect_length(u$index, 2L * nlevels(dd$g))
  # and it is group-major over the union: each group's two columns adjacent
  expect_identical(u$index[1:2],
                   c(u$pieces[[1L]]$index[1L], u$pieces[[2L]]$index[1L]))
})

test_that("term_tags_deep reaches a label inside a subformula", {
  # the labelled term is not one of the equation's terms: it develops the
  # break-point of a seg(), so the walk has to go through term_components().
  # Under an ADDITIVE parent that label is fitted; the walk is what tells the
  # two parents apart, so it is tested on its own.
  tm <- modelterms7::term_build(
    modelterms7::seg(x, psi ~ modelterms7::random(~ 1 | u | g)), dd)
  lab <- term_tags_deep(tm)
  expect_identical(unname(lab), "u")
  expect_match(names(lab), "psi1", fixed = TRUE)
  # and nothing under an unlabelled one
  expect_length(term_tags_deep(modelterms7::term_build(
    modelterms7::seg(x, psi ~ modelterms7::random(~ 1 | g)), dd)), 0L)
})

test_that("an unlabelled model is untouched", {
  # the guard reads a property off the term, so a formula that carried no
  # label before still fits, and the walk costs nothing where there is none
  expect_length(term_tags_deep(
    modelterms7::term_build(modelterms7::random(~ 1 | g), dd)), 0L)
  expect_length(term_tags_deep(
    modelterms7::term_build(modelterms7::seg(x, psi ~ random(~ 1 | g)), dd)), 0L)
  expect_no_error(statmod_spec(y ~ x + random(~ 1 | g),
                               distributions7::gaussian1_distrib(), dd))
})


test_that("a class's index is the interleaving its penalty reads", {
  # penalties7 reshapes a blockwise penalty's argument by row, so the index
  # must list, for each group in turn, that group's columns from every member.
  # Getting it wrong is not a rounding: the value is a different number.
  spec <- statmod_spec(y ~ random(~ 1 | u | g) | sigma ~ random(~ 1 | u | g),
                       distributions7::gaussian1_distrib(), dd)
  des <- statmod_design(spec)
  u <- statmod_penalized(spec, des)[[1L]]
  m <- u$class$m

  # the joint penalty against the form written out by hand, at a covariance
  # built from the chart itself so the two describe the same matrix exactly
  dr <- parameters7::dr_prod(2)
  v <- c(log(0.8), log(1.3), 0.4)
  names(v) <- dr@free_names
  Sig <- parameters7::param_value(dr, v)
  th <- as.list(stats::setNames(v, u$penalty@params))
  set.seed(4)
  b <- stats::rnorm(length(u$index))
  B <- matrix(b, ncol = 2, byrow = TRUE)
  manual <- 0.5 * sum(diag(B %*% solve(Sig) %*% t(B))) +
    (m / 2) * log(det(Sig)) + m * log(2 * pi)
  expect_equal(penalties7::penalty_value(u$penalty, b, th), manual,
               tolerance = 1e-12)

  # the negative control: swapping the two members inside each group is a
  # different number, so the test could fail
  bad <- as.numeric(t(B[, 2:1]))
  expect_false(isTRUE(all.equal(penalties7::penalty_value(u$penalty, bad, th),
                                manual, tolerance = 1e-6)))
})

test_that("a class of one member is today's random effect exactly", {
  # a single-member class of one column has no correlation to carry, so its
  # default prior is the centered univariate gaussian an unlabelled term
  # builds, and the two fits are the same fit
  set.seed(202)
  m <- 25; ni <- 8
  d2 <- data.frame(id = factor(rep(seq_len(m), each = ni)),
                   x = stats::rnorm(m * ni))
  d2$y <- 1 + 0.5 * d2$x + stats::rnorm(m, sd = 0.7)[d2$id] +
    stats::rnorm(m * ni, sd = 0.5)
  a <- statmod(y ~ x + random(~ 1 | id),
               distributions7::gaussian1_distrib(), d2, outer_criterion = reml())
  b <- statmod(y ~ x + random(~ 1 | u | id),
               distributions7::gaussian1_distrib(), d2, outer_criterion = reml())
  expect_identical(as.numeric(logLik(a)), as.numeric(logLik(b)))
  expect_identical(unlist(coef(a), use.names = FALSE),
                   unlist(coef(b), use.names = FALSE))
  expect_identical(hyper(a)$estimate, hyper(b)$estimate)
  # what differs is only what it is filed under
  expect_identical(hyper(b)$term, "u | id")
})

test_that("a shared block is refused where it would mean nothing", {
  d2 <- dd
  d2$g2 <- factor(rep(c("p", "q"), length.out = nrow(dd)))
  # one label, two groupings: the effects are indexed by different things
  expect_error(
    statmod_spec(y ~ random(~ 1 | u | g) | sigma ~ random(~ 1 | u | g2),
                 distributions7::gaussian1_distrib(), d2),
    "more than one grouping")
  # the joint prior is one object, so it is named once or not at all
  mv <- do.call(distributions7::fixed,
                list(distributions7::mvgaussian1_distrib(2), mu1 = 0, mu2 = 0))
  expect_error(
    statmod_spec(y ~ random(~ 1 | u | g, distrib = mv) +
                   random(~ 0 + x | u | g, distrib = mv),
                 distributions7::gaussian1_distrib(), dd),
    "named on one term or on none")
  # and its dimension is the class's total
  mv3 <- do.call(distributions7::fixed,
                 list(distributions7::mvgaussian1_distrib(3),
                      mu1 = 0, mu2 = 0, mu3 = 0))
  expect_error(
    statmod_spec(y ~ random(~ 1 | u | g, distrib = mv3) +
                   random(~ 0 + x | u | g),
                 distributions7::gaussian1_distrib(), dd),
    "collects 2")
})

test_that("a label inside an additive subformula is fitted, not refused", {
  # what L4 changed: under an additive parent the labelled effect occupies
  # columns of that term's block, so it sits in the same vector as one written
  # in an equation and the class reaches it with no new machinery
  spec <- statmod_spec(y ~ seg(x, psi ~ random(~ 1 | u | g)),
                       distributions7::gaussian1_distrib(), dd)
  u <- Filter(function(z) !is.null(z$pieces),
              statmod_penalized(spec, statmod_design(spec)))
  expect_length(u, 1L)
  expect_identical(u[[1L]]$key, "u | g")
  expect_false(is.null(u[[1L]]$pieces[[1L]]$within))
})


# What a fit reports about a covariance class ------------------------------

test_that("hyper() names every equation a class spans", {
  set.seed(404)
  m <- 20; ni <- 10
  d2 <- data.frame(id = factor(rep(seq_len(m), each = ni)))
  u <- matrix(stats::rnorm(2 * m, sd = 0.6), m, 2)
  d2$y <- stats::rnorm(m * ni, mean = 1 + u[d2$id, 1],
                       sd = exp(-0.7 + u[d2$id, 2]))
  fit <- statmod(y ~ random(~ 1 | b | id) | sigma ~ random(~ 1 | b | id),
                 distributions7::gaussian1_distrib(), d2,
                 outer_criterion = reml())
  h <- hyper(fit)
  # a class has no single parameter: `param` on its unit is the first
  # member's, which is a convention the store is keyed by and not a fact
  expect_true(all(h$parameter == "mu, sigma"))
  expect_true(all(h$term == "b | id"))
  expect_identical(nrow(h), 3L)

  # the edf total is the trace of the assembled smoother, cross block and all
  des <- statmod_design(fit@spec)
  H <- as.matrix(statmod_information_at(fit@spec, fit@coefficients, des,
                                        TRUE, "bartlett"))
  S <- as.matrix(statmod_penalty_at(fit@spec, fit@coefficients, fit@hyper,
                                    des, "hessian"))
  expect_equal(sum(fit@edf$edf), sum(diag(solve(H + S) %*% H)),
               tolerance = 1e-8)
})

test_that("a labelled random effect is reported as a random effect", {
  set.seed(505)
  m <- 15; ni <- 10
  d2 <- data.frame(id = factor(rep(seq_len(m), each = ni)))
  d2$y <- stats::rnorm(m * ni) + stats::rnorm(m, sd = 0.6)[d2$id]
  fit <- statmod(y ~ random(~ 1 | b | id), distributions7::gaussian1_distrib(),
                 d2, outer_criterion = reml())
  # it declares no penalty, so read after the penalties it came back
  # "parametric" and its grouping indicators were printed one per line
  expect_identical(term_block_kind(fit@spec@terms$mu[["random(~1 | b | id)"]]),
                   "random")
  out <- utils::capture.output(print(summary(fit)))
  expect_true(any(grepl("^random", out)))
  expect_false(any(grepl("random.1 ", out, fixed = TRUE)))
  # and the class's hyperparameter is reported, which reading the term's own
  # penalties could not do
  expect_true(any(grepl("sigma", out)))
})

test_that("a shared block is announced in a note, and a private one is not", {
  set.seed(606)
  m <- 20; ni <- 10
  d2 <- data.frame(id = factor(rep(seq_len(m), each = ni)))
  u <- matrix(stats::rnorm(2 * m, sd = 0.6), m, 2)
  d2$y <- stats::rnorm(m * ni, mean = 1 + u[d2$id, 1],
                       sd = exp(-0.7 + u[d2$id, 2]))
  shared <- statmod(y ~ random(~ 1 | b | id) | sigma ~ random(~ 1 | b | id),
                    distributions7::gaussian1_distrib(), d2,
                    outer_criterion = reml())
  nt <- summary(shared)@notes
  expect_true(any(grepl("covariance block 'b | id' is shared", nt, fixed = TRUE)))
  expect_true(any(grepl("'mu'", nt, fixed = TRUE) & grepl("'sigma'", nt)))
  # the hyperparameters are printed ONCE, under the first member
  out <- utils::capture.output(print(summary(shared)))
  expect_identical(sum(grepl("cor_v1_v2", out)), 1L)
  # a class of one member shares nothing and gets no note
  alone <- statmod(y ~ random(~ 1 | b | id),
                   distributions7::gaussian1_distrib(), d2,
                   outer_criterion = reml())
  expect_false(any(grepl("is shared", summary(alone)@notes)))
})


# A label written inside a subformula --------------------------------------

test_that("a label inside an additive term's subformula joins the class", {
  set.seed(88)
  m <- 8; ni <- 15
  d2 <- data.frame(id = factor(rep(seq_len(m), each = ni)),
                   x = rep(seq(-3, 3, length.out = ni), m))
  d2$y <- stats::rnorm(m * ni)
  spec <- statmod_spec(y ~ seg(x, psi ~ random(~ 1 | u | id)) |
                         sigma ~ random(~ 1 | u | id),
                       distributions7::gaussian1_distrib(), d2)
  des <- statmod_design(spec)
  u <- Filter(function(z) !is.null(z$pieces), statmod_penalized(spec, des))
  expect_length(u, 1L)
  u <- u[[1L]]
  expect_identical(u$class$dim, 2L)

  # the sub-term's piece is PART of its parent's block, named by `within`,
  # where the equation-level one is the whole of its term
  sub <- Filter(function(z) !is.null(z$within), u$pieces)
  top <- Filter(function(z) is.null(z$within), u$pieces)
  expect_length(sub, 1L)
  expect_length(top, 1L)
  expect_identical(sub[[1L]]$cols,
                   des$mu$blocks[[sub[[1L]]$term]][sub[[1L]]$within])
  expect_length(sub[[1L]]$cols, m)

  # and the two are interleaved group by group, as the penalty reads them
  expect_identical(u$index[1:2],
                   c(u$pieces[[1L]]$index[1L], u$pieces[[2L]]$index[1L]))
  expect_length(u$index, 2L * m)
})

test_that("a label under a structural term is refused for its own reason", {
  set.seed(89)
  m <- 6; ni <- 20
  d2 <- data.frame(id = factor(rep(seq_len(m), each = ni)))
  d2$y <- stats::rnorm(m * ni)
  err <- tryCatch(statmod_spec(y ~ gas(p = 1, q = 1,
                                       omega ~ random(~ 1 | u | id), by = id),
                               distributions7::gaussian1_distrib(), d2),
                  error = conditionMessage)
  # its coefficients are the term's own parameters and contribute no design
  # column, so the class's index has nowhere to point
  expect_match(err, "structural term", fixed = TRUE)
  expect_match(err, "covariance label 'u'", fixed = TRUE)
  expect_match(err, "under an ordinary term is fitted", fixed = TRUE)
})

test_that("a parent keeps its own penalty beside a labelled sub-term", {
  set.seed(90)
  m <- 8; ni <- 15
  d2 <- data.frame(id = factor(rep(seq_len(m), each = ni)),
                   x = rep(seq(-3, 3, length.out = ni), m),
                   g = factor(rep(c("a", "b"), length.out = m * ni)))
  d2$y <- stats::rnorm(m * ni)
  spec <- statmod_spec(
    y ~ seg(x, gamma1 ~ 0 + ridge(~ g), psi ~ random(~ 1 | u | id)) |
      sigma ~ random(~ 1 | u | id),
    distributions7::gaussian1_distrib(), d2)
  us <- statmod_penalized(spec, statmod_design(spec))
  # two units: the term's own ridge and the class, and the class does not
  # swallow the first
  expect_length(us, 2L)
  expect_length(Filter(function(z) is.null(z$pieces), us), 1L)
  expect_length(Filter(function(z) !is.null(z$pieces), us), 1L)
})
