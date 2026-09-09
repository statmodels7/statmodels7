# The degrees of freedom of a model carrying a filter -----------------------
#
# The effective count is the trace of the model's smoother over the vector the
# model estimates, which where a filter is present is the coefficients AND the
# term's own free parameters. The rule this replaced -- one degree of freedom
# per free parameter of a structural term, on the reading that those are
# estimated and unpenalized -- is the same trace restricted to the coordinates
# no penalty covers, and the tests below pin both halves of that sentence: the
# unpenalized coordinates still cost one, and a coordinate a prior shrinks
# costs what the shrinkage leaves.

panel <- function(seed = 51L, m = 6L, ni = 30L) {
  set.seed(seed)
  g <- factor(rep(seq_len(m), each = ni))
  al <- 0.25 * exp(stats::rnorm(m, 0, 0.5))
  bi <- stats::rnorm(m, 0, 0.6)
  y <- unlist(lapply(seq_len(m), function(i) {
    f <- 0
    out <- numeric(ni)
    for (t in seq_len(ni)) {
      mu <- bi[i] + f
      out[t] <- stats::rnorm(1, mu, 1)
      f <- al[i] * (out[t] - mu) + 0.6 * f
    }
    out
  }))
  data.frame(g = g, y = y)
}

# the hyperparameters are HELD: what is under test is the reading, not the
# search, and holding them keeps every fit here to a couple of seconds
fit_panel <- function(form, dd) {
  statmod(form, distributions7::gaussian1_distrib(), dd,
          outer_criterion = NULL)
}

test_that("an unpenalized parameter of a filter still costs exactly one", {
  dd <- panel()
  fit <- fit_panel(y ~ gas(p = 1, q = 1, by = g,
                           alpha1 ~ 1 + random(~ 1 | g)), dd)
  design <- statmod_design(fit@spec)
  js <- joint_smoother_diag(fit@spec, fit@coefficients, design, fit@hyper)
  expect_false(is.null(js))

  sst <- statmod_structural_state(design)
  key <- names(sst$zeta)[[1L]]
  free <- setdiff(names(sst$zeta[[key]]), sst$held[[key]])
  pen <- structural_range_cols(fit@spec, design, key, free)
  unp <- setdiff(seq_along(free), pen)
  expect_true(length(unp) > 0L)
  expect_true(length(pen) > 0L)

  # S has a zero column outside its own support, so F = I - M^-1 S has a unit
  # diagonal there whatever the penalty does. A TOLERANCE and not an identity:
  # the diagonal is read off a computed inverse rather than off that
  # cancellation, so the last bits are the platform's.
  expect_equal(unname(js$zeta[unp]), rep(1, length(unp)), tolerance = 1e-10)

  # THE NEGATIVE CONTROL, without which the test above is satisfied by any
  # route that returns ones: a coordinate the prior shrinks costs less than
  # one, and the old rule counted it as a whole.
  expect_lt(sum(js$zeta[pen]), length(pen) - 0.5)
  expect_true(all(js$zeta[pen] < 0.9))
})

test_that("a structural term's count is the joint trace over its parameters", {
  dd <- panel()
  fit <- fit_panel(y ~ gas(p = 1, q = 1, by = g,
                           alpha1 ~ 1 + random(~ 1 | g)), dd)
  design <- statmod_design(fit@spec)
  sst <- statmod_structural_state(design)
  key <- names(sst$zeta)[[1L]]
  free <- setdiff(names(sst$zeta[[key]]), sst$held[[key]])
  nb <- sum(vapply(design[fit@spec@distrib@params], function(d) d$npar,
                   integer(1)))

  # the reference goes through a plain solve rather than solve_pd, so the two
  # sides share the matrices and not the arithmetic that inverts them
  K <- statmod_full_information(fit@spec, fit@coefficients, design)
  M <- statmod_marginal_full(fit@spec, design, fit@coefficients, fit@hyper)
  d <- diag(solve(as_dense(M), as_dense(K)))

  row <- fit@edf[fit@edf$term != "linpar", ]
  expect_equal(nrow(row), 1L)
  expect_equal(row$edf, sum(d[nb + seq_along(free)]), tolerance = 1e-8)
  expect_equal(sum(fit@edf$edf), sum(d), tolerance = 1e-8)

  # and it is strictly less than the count it replaced, which is what says the
  # change is doing something on this model
  expect_lt(row$edf, length(free) - 0.5)
})

test_that("a model with no filter keeps the coefficient-only reading", {
  set.seed(11)
  n <- 300
  dd <- data.frame(x = stats::runif(n, -3, 3))
  dd$y <- stats::rnorm(n, sin(dd$x), 0.5)
  fit <- statmod(y ~ s(x, bspline_smooth(k = 10)), distributions7::gaussian1_distrib(), dd,
                 outer_criterion = reml())
  design <- statmod_design(fit@spec)
  # there is no joint vector to read, so the route is not even offered
  # NULL, and not the `failed` marker: there is no joint vector here at all,
  # which is a different answer from one that could not be read
  expect_null(joint_smoother_diag(fit@spec, fit@coefficients, design,
                                  fit@hyper))

  H <- statmod_information_at(fit@spec, fit@coefficients, design, TRUE, "opg")
  S <- zap_nonfinite(statmod_penalty_at(fit@spec, fit@coefficients, fit@hyper,
                                        design, "hessian"))
  d <- diag(solve(as_dense(H + S), as_dense(H)))
  expect_equal(sum(fit@edf$edf), sum(d), tolerance = 1e-8)
})

test_that("a term that mixes over states keeps one apiece, and is right to", {
  # A regime's contribution is a likelihood mixed over latent states, so
  # statmod_marginal_full() has no determinant for it and the joint reading is
  # not available. It costs nothing: regime() takes no subformula, so its own
  # parameters carry no penalty, and one apiece is exactly what the joint
  # trace would say for an unpenalized coordinate.
  # written out here rather than borrowed: a test file is sourced on its own,
  # so a global defined in another one is not in scope
  P <- matrix(c(0.93, 0.07, 0.06, 0.94), 2, 2, byrow = TRUE)
  set.seed(21)
  n <- 300
  st <- integer(n)
  st[1L] <- 1L
  for (i in 2:n) st[i] <- sample.int(2L, 1L, prob = P[st[i - 1L], ])
  dd <- data.frame(t = seq_len(n), y = stats::rnorm(n, c(0, 3)[st], 1),
                   x = stats::rnorm(n))
  fit <- statmod(y ~ x + regime(2, time = t),
                 distributions7::gaussian1_distrib(), dd)
  design <- statmod_design(fit@spec)
  # NULL, and not the `failed` marker: there is no joint vector here at all,
  # which is a different answer from one that could not be read
  expect_null(joint_smoother_diag(fit@spec, fit@coefficients, design,
                                  fit@hyper))
  zn <- modelterms7::term_params(fit@spec@terms$mu[[2L]])
  row <- fit@edf[fit@edf$term != "linpar" & fit@edf$parameter == "mu", ]
  expect_equal(row$edf, as.numeric(length(zn) - 1L))
})

test_that("a class counts the coordinates it collects, not its members", {
  dd <- panel()
  fit <- fit_panel(y ~ gas(p = 1, q = 1, by = g,
                           alpha1 ~ 1 + random(~ 1 | u | g),
                           omega ~ 1 + random(~ 1 | u | g)), dd)
  s <- summary(fit)
  expect_length(s@classes, 1L)
  ed <- as.numeric(sub("^edf ", "",
                       grep("^edf ", s@classes[[1L]]$bits, value = TRUE)))
  expect_length(ed, 1L)

  design <- statmod_design(fit@spec)
  u <- Filter(function(q) !is.null(q$pieces) && length(q$pieces) > 1L,
              statmod_penalized(fit@spec, design))[[1L]]
  js <- joint_smoother_diag(fit@spec, fit@coefficients, design, fit@hyper)
  pos <- unit_joint_positions(u, fit@spec, design)
  expect_length(pos, 2L * nlevels(dd$g))
  expect_equal(ed, sum(c(js$beta, js$zeta)[pos]), tolerance = 5e-3)

  # THE CLASS DOES NOT COLLECT THE WHOLE TERM: the level and the persistence
  # are the term's and not the class's, so the class's count is strictly the
  # smaller of the two. Printed as the sum of the members' rows it was the
  # term's own number.
  term_row <- fit@edf$edf[fit@edf$term != "linpar"]
  expect_length(term_row, 1L)
  expect_lt(ed, term_row - 0.5)
})

test_that("a count that cannot be read is missing, not the old rule", {
  # WHERE THE JOINT MATRIX CANNOT BE READ the one-apiece rule is a different
  # quantity, and returning it under the same name reported 20.00 on a model
  # whose penalized count is 10.40. Measured, the refusal happens where that
  # matrix is INDEFINITE -- smallest eigenvalue -7.1e-05, chol() refusing it
  # too -- so there is nothing to relax and the count is simply absent.
  #
  # Driven through a stub rather than through a data set that happens to reach
  # such a point: which panels do is a property of the sample, and a test that
  # depended on one would be a test of the sample.
  dd <- panel()
  fit <- fit_panel(y ~ gas(p = 1, q = 1, by = g,
                           alpha1 ~ 1 + random(~ 1 | g)), dd)
  design <- statmod_design(fit@spec)
  orig <- joint_smoother_diag
  on.exit(utils::assignInNamespace("joint_smoother_diag", orig,
                                   ns = "statmodels7"), add = TRUE)

  # the reading as it stands, which the stub must differ from
  ok <- statmod_edf(fit@spec, fit@coefficients, design, fit@hyper)
  row_ok <- ok$edf[ok$term != "linpar"]
  expect_length(row_ok, 1L)
  expect_true(is.finite(row_ok))

  utils::assignInNamespace(
    "joint_smoother_diag",
    function(spec, coef, design, hyper) list(beta = NULL, zeta = NULL,
                                             failed = TRUE),
    ns = "statmodels7")
  # the value is captured by assignment: in the third edition expect_warning()
  # returns the CONDITION and not the value of its argument
  bad <- NULL
  expect_warning(
    bad <- statmod_edf(fit@spec, fit@coefficients, design, fit@hyper),
    "cannot be read")
  row_bad <- bad$edf[bad$term != "linpar"]
  expect_length(row_bad, 1L)
  expect_true(is.na(row_bad))
  # and it is NOT the count of the term's free parameters, which is what the
  # rule it replaces returned
  zn <- modelterms7::term_params(fit@spec@terms$mu[["gas(p = 1, q = 1, by = g, alpha1 ~ 1 + random(~1 | g))"]])
  expect_false(isTRUE(all.equal(row_bad, as.numeric(length(zn) - 1L))))
  # the ordinary rows still say what they can
  expect_true(all(is.finite(bad$edf[bad$term == "linpar"])))
})


test_that("logLik's fallback for a missing count is an upper bound", {
  # logLik() replaces a count it could not obtain by the term's number of
  # COLUMNS, on the stated reading that this is an upper bound and that
  # dropping it would flatter every criterion. A structural term has no
  # columns, so that bound is zero and says the opposite of what it means:
  # measured, a filter carrying eighteen free parameters reported a df of 2.
  dd <- panel()
  fit <- fit_panel(y ~ gas(p = 1, q = 1, by = g,
                           alpha1 ~ 1 + random(~ 1 | g)), dd)
  nm <- fit@edf$term[fit@edf$term != "linpar"]
  expect_length(nm, 1L)
  expect_identical(fit@edf$coefficients[fit@edf$term == nm], 0L)

  design <- statmod_design(fit@spec)
  sp <- statmod_structural_par(fit@spec, design)
  free <- length(sp[[nm]]$parameter) - length(sp[[nm]]$held)
  expect_gt(free, 1L)

  broken <- fit
  broken@edf$edf[broken@edf$term == nm] <- NA_real_
  df <- attr(stats::logLik(broken), "df")
  # the free parameters plus the two intercepts, and NOT the two alone
  expect_equal(df, free + 2)
  expect_gt(df, 2)
})


test_that("a class of ordinary members counts the same either way", {
  # THE CONTROL for the test above. Where every member is an ordinary term
  # whose block the class covers entirely, the trace over the class's own
  # coordinates IS the members' rows added up, so the reading that separates
  # the two above must not separate them here.
  set.seed(707)
  m <- 12L
  ni <- 10L
  g <- factor(rep(seq_len(m), each = ni))
  b <- stats::rnorm(m, 0, 1.2)
  u <- stats::rnorm(m, 0, 0.35)
  dd <- data.frame(g = g)
  dd$y <- stats::rnorm(m * ni, mean = 1 + b[as.integer(g)],
                       sd = exp(-0.3 + u[as.integer(g)]))
  fit <- statmod(y ~ random(~ 1 | a | g) | sigma ~ random(~ 1 | a | g),
                 distributions7::gaussian1_distrib(), dd,
                 outer_criterion = reml())
  s <- summary(fit)
  ed <- as.numeric(sub("^edf ", "",
                       grep("^edf ", s@classes[[1L]]$bits, value = TRUE)))
  members <- sum(fit@edf$edf[fit@edf$term == "random(~1 | a | g)"])
  expect_equal(ed, members, tolerance = 5e-3)
})
