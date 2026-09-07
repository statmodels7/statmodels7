# A covariance class inside a structural term ------------------------------
#
# A filter's own parameters are estimated beside the coefficients, so a
# labelled random effect developing one of them carries a prior like any
# other. What differs is only where its coordinates are: among the term's
# own parameters, which contribute no design column, rather than in the
# stacked coefficient vector.

# a panel whose loading varies by group, which is the model the case exists
# for: alpha_i lognormal around 0.25, a persistence of 0.6 and a random
# intercept per group
panel <- function(seed = 51L, m = 8L, ni = 40L) {
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

test_that("a class of one structural member is the unlabelled fit exactly", {
  # THE CONTROL THAT SAYS THE ADDRESSING IS RIGHT. A label collecting one
  # member describes the same prior over the same coordinates as no label at
  # all, so every number has to come back the same -- and it does so through
  # the whole of the new route: a class is assembled, its penalty built from
  # the class's own default, and its value read at interleaved positions in
  # the term's parameter vector.
  dd <- panel()
  one <- statmod(y ~ gas(p = 1, q = 1, by = g,
                         alpha1 ~ 1 + random(~ 1 | u | g)),
                 distributions7::gaussian1_distrib(), dd,
                 outer_criterion = reml())
  bare <- statmod(y ~ gas(p = 1, q = 1, by = g,
                          alpha1 ~ 1 + random(~ 1 | g)),
                  distributions7::gaussian1_distrib(), dd,
                  outer_criterion = reml())
  expect_equal(as.numeric(stats::logLik(one)),
               as.numeric(stats::logLik(bare)), tolerance = 1e-12)
  expect_equal(one@coefficients, bare@coefficients, tolerance = 1e-12)
  expect_equal(sum(one@edf$edf), sum(bare@edf$edf), tolerance = 1e-12)
  expect_equal(unname(one@structural[[1L]]$parameter),
               unname(bare@structural[[1L]]$parameter), tolerance = 1e-12)
  # the hyperparameter is the same number under a different key: the class's,
  # since the prior belongs to the class and not to the member
  expect_equal(hyper(one)$estimate, hyper(bare)$estimate, tolerance = 1e-12)
  expect_identical(hyper(one)$term, "u | g")
})

test_that("the penalty a one-member class builds is the member's own", {
  # the mechanism behind the control above, checked directly so that a
  # failure there says which half moved
  dd <- panel()
  bare <- modelterms7::term_build(
    modelterms7::gas(p = 1, q = 1, by = g, alpha1 ~ 1 + random(~ 1 | g)), dd)
  lab <- modelterms7::term_build(
    modelterms7::gas(p = 1, q = 1, by = g, alpha1 ~ 1 + random(~ 1 | u | g)),
    dd)
  # a LABELLED sub-term declares no penalty of its own: the class carries it
  expect_length(modelterms7::term_penalties(lab), 0L)
  own <- modelterms7::term_penalties(bare)[[1L]]$penalty
  cls <- statmod_classes(list(mu = list(T = lab)))[[1L]]$penalty
  b <- stats::rnorm(nlevels(dd$g))
  th <- list(sigma = 0.7)
  expect_identical(penalties7::penalty_value(own, b, th),
                   penalties7::penalty_value(cls, b, th))
  expect_identical(penalties7::penalty_gradient(own, b, th),
                   penalties7::penalty_gradient(cls, b, th))
})

test_that("two developed parameters of one filter share a covariance", {
  # the case the work exists for: the level and the loading of one filter
  # vary by group and the two deviations are correlated
  dd <- panel()
  fit <- statmod(y ~ gas(p = 1, q = 1, by = g,
                         alpha1 ~ 1 + random(~ 1 | u | g),
                         omega ~ 1 + random(~ 1 | u | g)),
                 distributions7::gaussian1_distrib(), dd,
                 outer_criterion = reml())
  des <- statmod_design(fit@spec)
  u <- Filter(function(z) !is.null(z$pieces),
              statmod_penalized(fit@spec, des))[[1L]]
  expect_true(u$structural)
  expect_identical(u$class$dim, 2L)
  expect_length(u$cols, 2L * nlevels(dd$g))
  # INTERLEAVED GROUP BY GROUP, which is the order the prior reads its
  # argument in: for each group, that group's coordinate from each member
  expect_identical(u$cols[1:2],
                   c(u$pieces[[1L]]$cols[[1L]], u$pieces[[2L]]$cols[[1L]]))
  # the prior is over the two coordinates of one group, repeated over groups
  expect_identical(as.integer(u$penalty@block), 2L)

  # and the report names each coordinate for the parameter it develops, not
  # for the equation alone: both are effects inside one term of `mu`
  cd <- class_coords(u)
  expect_identical(cd$column,
                   c("omega:(Intercept)", "alpha1:(Intercept)"))
  s <- summary(fit)
  expect_length(s@classes, 1L)
  # the class's own count is the term's, taken ONCE: both members are
  # developments of the same term, and a sum over the pieces read its row
  # twice
  bits <- s@classes[[1L]]$bits
  ed <- as.numeric(sub("^edf ", "", grep("^edf ", bits, value = TRUE)))
  expect_lte(ed, sum(fit@edf$edf))
})

test_that("two classes in one filter accumulate into one vector", {
  # Both are addressed in the SAME term's parameters, and the accumulator is
  # built once per term, so the second must add to it rather than reset it.
  dd <- panel(m = 8L, ni = 20L)
  spec <- statmod_spec(y ~ 0 + gas(p = 1, q = 1, by = g,
                                   omega ~ 1 + random(~ 1 | u | g),
                                   alpha1 ~ 1 + random(~ 1 | v | g)),
                       distributions7::gaussian1_distrib(), dd)
  des <- statmod_design(spec)
  us <- Filter(function(z) !is.null(z$pieces), statmod_penalized(spec, des))
  expect_length(us, 2L)
  expect_true(all(vapply(us, function(z) isTRUE(z$structural), TRUE)))
  # two classes, disjoint coordinates, one term
  expect_length(intersect(us[[1L]]$cols, us[[2L]]$cols), 0L)
  expect_identical(us[[1L]]$term, us[[2L]]$term)

  # displaced by hand, because the START puts one class's coordinates at
  # exactly zero -- where a centered gaussian prior's gradient IS zero -- and
  # a check read there would pass with the second class doing all the work
  sst <- statmod_structural_state(des)
  v <- sst$zeta[[1L]]
  v[us[[1L]]$cols] <- seq(-0.4, 0.4, length.out = length(us[[1L]]$cols))
  v[us[[2L]]$cols] <- seq(0.3, -0.3, length.out = length(us[[2L]]$cols))
  sst$zeta[[1L]] <- v
  sst$key <- NULL
  sst$value <- NULL
  hy <- statmod_hyper_start(spec, des)
  g <- statmod_structural_penalty(spec, des, hy, "gradient")[[1L]]
  h <- statmod_structural_penalty(spec, des, hy, "hessian")[[1L]]
  expect_true(all(g[us[[1L]]$cols] != 0))
  expect_true(all(g[us[[2L]]$cols] != 0))
  expect_identical(sort(unname(which(rowSums(abs(h)) > 0))),
                   sort(c(us[[1L]]$cols, us[[2L]]$cols)))
  # at sigma = 1 a centered gaussian's gradient is the coordinate itself
  expect_equal(unname(g[us[[1L]]$cols]), v[us[[1L]]$cols],
               ignore_attr = TRUE, tolerance = 1e-12)
})

test_that("a model with no label is untouched by the class route", {
  # the negative control: the unit a class produces is new, but the two sites
  # a class reaches -- the penalty's value and the structural gradient -- are
  # on the path of every fit
  set.seed(7)
  n <- 300L
  g <- factor(rep(seq_len(15L), each = 20L))
  x <- stats::rnorm(n)
  b <- stats::rnorm(15L, 0, 0.5)
  d2 <- data.frame(y = 1 + 0.8 * x + b[as.integer(g)] + stats::rnorm(n, 0, 0.7),
                   x = x, g = g)
  fit <- statmod(y ~ x + random(~ 1 | g),
                 distributions7::gaussian1_distrib(), d2,
                 outer_criterion = reml())
  us <- statmod_penalized(fit@spec, statmod_design(fit@spec))
  expect_length(Filter(function(z) !is.null(z$pieces), us), 0L)
  expect_false(any(vapply(us, function(z) isTRUE(z$structural), TRUE)))
  expect_true(fit@converged)
})
