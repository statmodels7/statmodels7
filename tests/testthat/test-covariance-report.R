# What a summary says about a covariance shared between terms -------------
#
# The chart's coordinates are numbered by the multivariate family, which is a
# law on R^d and has no model to name them after. Which equation, which term
# and which within-group column each one is, is this layer's answer, and
# without it a printed correlation between `v1` and `v3` says nothing.

test_that("a coordinate is named for the equation and the column it is", {
  set.seed(11)
  m <- 12L
  ni <- 6L
  dd <- data.frame(g = factor(rep(seq_len(m), each = ni)),
                   x = stats::rnorm(m * ni))
  dd$y <- stats::rnorm(m * ni)

  spec <- statmod_spec(y ~ x + random(~ 1 + x | a | g) | sigma ~ random(~ 1 | a | g),
                       distributions7::gaussian1_distrib(), dd)
  des <- statmod_design(spec)
  u <- Filter(function(z) !is.null(z$pieces), statmod_penalized(spec, des))
  expect_length(u, 1L)
  cd <- class_coords(u[[1L]])

  # THE ORDER IS THE CLASS'S OWN: the pieces in the order the equations were
  # walked, each contributing its within-group columns in order, which is how
  # class_index() interleaves them
  expect_identical(cd$param, c("mu", "mu", "sigma"))
  expect_identical(cd$column, c("(Intercept)", "x", "(Intercept)"))
  expect_identical(cd$label,
                   c("mu:(Intercept)", "mu:x", "sigma:(Intercept)"))
})

test_that("the label carries the term where the equation is not enough", {
  set.seed(12)
  m <- 12L
  ni <- 6L
  dd <- data.frame(g = factor(rep(seq_len(m), each = ni)),
                   x = stats::rnorm(m * ni))
  dd$y <- stats::rnorm(m * ni)

  # two terms of ONE equation, whose columns are told apart by name
  spec <- statmod_spec(y ~ random(~ 1 | a | g) + random(~ 0 + x | a | g),
                       distributions7::gaussian1_distrib(), dd)
  des <- statmod_design(spec)
  u <- Filter(function(z) !is.null(z$pieces), statmod_penalized(spec, des))
  cd <- class_coords(u[[1L]])
  expect_identical(cd$label, c("mu:(Intercept)", "mu:x"))
  expect_identical(cd$term,
                   c("random(~1 | a | g)", "random(~0 + x | a | g)"))

  # and where they are NOT, the term goes into the label rather than two
  # coordinates being given one name
  cd2 <- coord_labels_of(data.frame(
    param = c("mu", "mu"), term = c("t1", "t2"), column = c("x", "x"),
    stringsAsFactors = FALSE))
  expect_identical(cd2, c("mu:t1:x", "mu:t2:x"))

  # one term alone needs no equation in front of its columns
  cd3 <- coord_labels_of(data.frame(
    param = c("mu", "mu"), term = c("t1", "t1"),
    column = c("(Intercept)", "x"), stringsAsFactors = FALSE))
  expect_identical(cd3, c("(Intercept)", "x"))
})

test_that("a label in a subformula names what the effect is an effect on", {
  set.seed(13)
  m <- 8L
  ni <- 12L
  d2 <- data.frame(id = factor(rep(seq_len(m), each = ni)),
                   x = rep(seq(-3, 3, length.out = ni), m))
  d2$y <- stats::rnorm(m * ni)
  spec <- statmod_spec(y ~ seg(x, psi ~ random(~ 1 | u | id)) |
                         sigma ~ random(~ 1 | u | id),
                       distributions7::gaussian1_distrib(), d2)
  des <- statmod_design(spec)
  u <- Filter(function(z) !is.null(z$pieces), statmod_penalized(spec, des))
  cd <- class_coords(u[[1L]])
  # THE FIRST COORDINATE IS NOT AN EFFECT ON THE MEAN. It is the random
  # intercept of the break-point, and `mu:(Intercept)` would be read as the
  # mean's own, which is a different quantity in the same equation.
  expect_identical(cd$label, c("mu:psi1:(Intercept)", "sigma:(Intercept)"))
  expect_identical(cd$term[[1L]], "seg(x, psi ~ random(~1 | u | id))")
})


test_that("only the shape the family names by position is rewritten", {
  lab <- c("mu:(Intercept)", "sigma:z")
  expect_identical(
    readable_coord_names(c("sd_v1", "sd_v2", "cor_v1_v2"), lab),
    c("sd[mu:(Intercept)]", "sd[sigma:z]",
      "cor[mu:(Intercept), sigma:z]"))
  # a prefix of its own is kept, whatever it is
  expect_identical(readable_coord_names("scale_sd_v2", lab),
                   "scale_sd[sigma:z]")
  # AND WHAT DOES NOT MATCH IS LEFT ALONE, which is what keeps a reading of
  # another kind from being renamed into a wrong one: a name with no index,
  # and one whose index no label answers for
  expect_identical(readable_coord_names(c("nu", "sd_v3", "lambda"), lab),
                   c("nu", "sd_v3", "lambda"))
  expect_identical(readable_coord_names("cor_v1_v2", character(0)),
                   "cor_v1_v2")
})


test_that("a term whose only penalized effect is labelled is not parametric", {
  # its own penalties are empty -- the class carries them -- and read after
  # the penalties it came back "parametric": twelve predictions printed one
  # per line with a z and a p beside them, the term's compartments gone
  set.seed(14)
  m <- 8L
  ni <- 10L
  dd <- data.frame(id = factor(rep(seq_len(m), each = ni)),
                   x = rep(seq(0.2, 4, length.out = ni), m))
  dd$y <- 5 * exp(-0.8 * dd$x) + stats::rnorm(m * ni, 0, 0.2)
  tm <- modelterms7::term_build(
    modelterms7::nl(~ a * exp(-r * x), a ~ 1 + random(~ 1 | u | id)), dd)
  expect_length(modelterms7::term_penalties(tm), 0L)
  expect_identical(term_block_kind(tm), "penalized")
  # and a term with neither a penalty nor a label under it still is
  expect_identical(term_block_kind(modelterms7::term_build(
    modelterms7::linpar(~ x), dd)), "parametric")
})


test_that("a flat direction the coefficients do not span is said to be one", {
  # the matrix spans a structural term's own parameters as well, and only the
  # coefficients are named: the list came out empty and the message read
  # "carried by: , "
  A <- diag(c(1, 1, 1e-14))
  expect_match(flat_directions(A, c("mu:x", "mu:z", "")),
               "none of the coefficients carry", fixed = TRUE)
  expect_false(grepl("carried by", flat_directions(A, c("mu:x", "mu:z", ""))))
  # where it is spanned, the coefficients are named as before
  expect_match(flat_directions(diag(c(1e-14, 1, 1)), c("mu:x", "mu:z", "")),
               "carried by: mu:x", fixed = TRUE)
})


test_that("a covariance inside a subformula is named from its own columns", {
  skip_on_cran()
  set.seed(909)
  m <- 12L
  ni <- 14L
  dd <- data.frame(id = factor(rep(seq_len(m), each = ni)),
                   x = rep(seq(0.2, 4, length.out = ni), m),
                   z = stats::rnorm(m * ni))
  a_i <- 5 + stats::rnorm(m, 0, 1)
  dd$y <- a_i[as.integer(dd$id)] * exp(-0.8 * dd$x) +
    stats::rnorm(m * ni, 0, 0.4)
  fit <- statmod(y ~ 0 + nl(~ a * exp(-r * x), a ~ 1 + random(~ 1 + z | id)),
                 distributions7::gaussian1_distrib(), dd,
                 outer_criterion = reml())
  b <- summary(fit)@tables$mu[[1L]]
  # THE COLUMNS ARE THE SUB-TERM'S. Asked of the parent, `nl()` has no
  # grouping and the coordinates stayed numbered.
  tb <- b$components[[1L]]$table
  expect_identical(tb$name[1:3],
                   c("sd[(Intercept)]", "sd[z]", "cor[(Intercept), z]"))
})


test_that("a shared covariance is reported once, ahead of the equations", {
  skip_on_cran()
  set.seed(707)
  m <- 25L
  ni <- 10L
  g <- factor(rep(seq_len(m), each = ni))
  # THE TWO SCALES ARE FAR APART, which is what makes the labelling
  # falsifiable: an order read backwards puts 0.35 where 1.2 belongs
  sd_mu <- 1.2
  sd_sg <- 0.35
  b <- stats::rnorm(m, 0, sd_mu)
  u <- stats::rnorm(m, 0, sd_sg)
  dd <- data.frame(g = g)
  dd$y <- stats::rnorm(m * ni, mean = 1 + b[as.integer(g)],
                       sd = exp(-0.3 + u[as.integer(g)]))
  fit <- statmod(y ~ random(~ 1 | a | g) | sigma ~ random(~ 1 | a | g),
                 distributions7::gaussian1_distrib(), dd,
                 outer_criterion = reml())
  s <- summary(fit)

  expect_length(s@classes, 1L)
  b1 <- s@classes[[1L]]
  expect_identical(b1$term, "a | g")
  expect_identical(b1$coords$label,
                   c("mu:(Intercept)", "sigma:(Intercept)"))
  expect_identical(b1$table$name,
                   c("sd[mu:(Intercept)]", "sd[sigma:(Intercept)]",
                     "cor[mu:(Intercept), sigma:(Intercept)]"))

  # THE LABELS ARE CHECKED AGAINST THE DATA and not against our own
  # bookkeeping: the mean's effects were drawn three times as wide as the
  # scale's, so a coordinate order read backwards fails here
  expect_gt(b1$table$estimate[[1L]], 0.7)
  expect_lt(b1$table$estimate[[2L]], 0.7)

  # and the numbers are the ones the chart implies, computed apart from the
  # summary
  eta <- unlist(fit@hyper$mu[["a | g"]])
  sig <- parameters7::param_value(parameters7::dr_prod(2L), eta)
  expect_equal(b1$table$estimate[[1L]], sqrt(sig[1, 1]), tolerance = 1e-8)
  expect_equal(b1$table$estimate[[2L]], sqrt(sig[2, 2]), tolerance = 1e-8)
  expect_equal(b1$table$estimate[[3L]],
               sig[1, 2] / sqrt(sig[1, 1] * sig[2, 2]), tolerance = 1e-8)

  # the members carry none of it, and say where it is instead of reporting
  # that there is nothing to report
  for (p in c("mu", "sigma")) {
    bl <- Filter(function(z) identical(z$kind, "random"), s@tables[[p]])
    expect_length(bl, 1L)
    expect_identical(nrow(bl[[1L]]$table), 0L)
    expect_match(bl[[1L]]$note, "covariance block 'a | g', reported above",
                 fixed = TRUE)
  }

  # the class's degrees of freedom are the members' added up
  expect_true(any(grepl(sprintf("edf %.2f", sum(fit@edf$edf[
    fit@edf$term == "random(~1 | a | g)"])), b1$bits)))
})


test_that("a class reaching a subformula is reported and both members say so", {
  skip_on_cran()
  set.seed(910)
  m <- 12L
  ni <- 14L
  dd <- data.frame(id = factor(rep(seq_len(m), each = ni)),
                   x = rep(seq(0.2, 4, length.out = ni), m))
  a_i <- 5 + stats::rnorm(m, 0, 1)
  g_i <- stats::rnorm(m, 0, 0.5)
  dd$y <- a_i[as.integer(dd$id)] * exp(-0.8 * dd$x) +
    stats::rnorm(m * ni, 0, exp(-1 + g_i[as.integer(dd$id)]))
  fit <- statmod(y ~ 0 + nl(~ a * exp(-r * x), a ~ 1 + random(~ 1 | u | id)) |
                   sigma ~ random(~ 1 | u | id),
                 distributions7::gaussian1_distrib(), dd,
                 outer_criterion = reml())
  s <- summary(fit)
  expect_length(s@classes, 1L)
  # the effect is on `a` and not on the mean, which is a different quantity
  # in the same equation
  expect_identical(s@classes[[1L]]$coords$label,
                   c("mu:a:(Intercept)", "sigma:(Intercept)"))

  out <- utils::capture.output(print(s))
  expect_identical(sum(grepl("reported above", out, fixed = TRUE)), 2L)
  # the nl term keeps its own reading: the compartment for `a`, the rows of
  # its other parameters, and no grouping indicator printed one per line
  expect_true(any(grepl("^  a  ~ ", out)))
  expect_false(any(grepl("nl.a.random.1", out, fixed = TRUE)))
  # AND A POINTER NEVER OUTLIVES THE SECTION IT POINTS AT
  expect_true(any(grepl("shared covariance blocks", out, fixed = TRUE)))
})


test_that("a covariance of one term stays in its block, renamed", {
  skip_on_cran()
  set.seed(808)
  m <- 25L
  ni <- 10L
  g <- factor(rep(seq_len(m), each = ni))
  x <- stats::rnorm(m * ni)
  b <- stats::rnorm(m, 0, 1.1)
  dd <- data.frame(g = g, x = x)
  dd$y <- 1 + 0.5 * x + b[as.integer(g)] + stats::rnorm(m * ni, 0, 0.6)
  fit <- statmod(y ~ x + random(~ 1 + x | g),
                 distributions7::gaussian1_distrib(), dd,
                 outer_criterion = reml())
  s <- summary(fit)
  # NOTHING IS SHARED, so nothing is lifted out: the section exists for a
  # block that belongs to no equation, and this one belongs to `mu`
  expect_length(s@classes, 0L)
  out <- utils::capture.output(print(s))
  expect_false(any(grepl("shared covariance blocks", out, fixed = TRUE)))

  bl <- Filter(function(z) identical(z$kind, "random"), s@tables$mu)[[1L]]
  expect_identical(bl$table$name[1:3],
                   c("sd[(Intercept)]", "sd[x]", "cor[(Intercept), x]"))
})
