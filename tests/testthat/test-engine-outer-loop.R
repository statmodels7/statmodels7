# The outer search in compiled code (engine_outer_search(): the search, the
# point and the penalty, phase 4) against the R route of outer_fit(): the
# same fit in every property, and the penalty built in C++ identical to the
# one R builds, signed zeros included.

loop_data <- function(n = 300, groups = 15, seed = 5) {
  set.seed(seed)
  g <- factor(rep_len(seq_len(groups), n))
  x <- runif(n)
  z <- runif(n, -1, 1)
  b <- matrix(rnorm(2 * groups, sd = 0.5), groups)
  data.frame(g = g, x = x, z = z,
             y = sin(5 * x) + b[g, 1] + b[g, 2] * z + rnorm(n, sd = 0.3))
}

# a fit on the R route of outer_fit() and on the compiled one
loop_fits <- function(...) {
  fl <- engine_flags
  old <- mget(c("outer_search", "outer_point", "outer_pen"), envir = fl)
  on.exit(for (k in names(old)) assign(k, old[[k]], envir = fl))
  for (k in names(old)) assign(k, FALSE, envir = fl)
  a <- statmod(..., engine = "auto")
  for (k in names(old)) assign(k, TRUE, envir = fl)
  b <- statmod(..., engine = "auto")
  list(r = a, c = b)
}

# identical values; a closure is compared without its environment, which
# belongs to each run
same_value <- function(a, b) {
  identical(a, b) || isTRUE(all.equal(a, b, tolerance = 0, check.environment = FALSE))
}

# TEMPORARY diagnostic: the size of a difference between the two routes
gap_report <- function(a, b) {
  if (same_value(a, b)) return("")
  one <- function(x, y) {
    x <- unlist(x); y <- unlist(y)
    if (!is.numeric(x) || !is.numeric(y)) return("non-numeric")
    if (length(x) != length(y)) return(sprintf("lengths %d vs %d", length(x), length(y)))
    sprintf("n %d, max abs %.3g, scale gap %.3g", length(x), max(abs(x - y), na.rm = TRUE), scale_gap(x, y))
  }
  if (is.list(a) && is.list(b) && !is.null(names(a)))
    paste(vapply(intersect(names(a), names(b)), function(k)
      paste0(k, ": ", one(a[[k]], b[[k]])), ""), collapse = "; ")
  else one(a, b)
}

expect_same_fit <- function(f) {
  for (p in setdiff(S7::prop_names(f$r), c("elapsed", "methods"))) {
    a <- S7::prop(f$c, p); b <- S7::prop(f$r, p)
    expect_true(same_value(a, b), info = paste(p, gap_report(a, b)))
  }
  # the optimizer's line search holds a closure over each run's frame: its
  # class and its rule are compared, the rest of the record exactly
  mr <- f$r@methods
  mc <- f$c@methods
  for (k in setdiff(names(mr), "search")) expect_true(same_value(mc[[k]], mr[[k]]), info = k)
  if (!is.null(mr$search)) {
    expect_identical(class(mc$search), class(mr$search))
    expect_identical(mc$search@criterion@label, mr$search@criterion@label)
  }
}

test_that("the compiled search reproduces the R route: smooth, random slope, two equations", {
  d <- loop_data()
  g1 <- distributions7::gaussian1_distrib()
  expect_same_fit(loop_fits(y ~ s(x), distrib = g1, data = d))
  expect_same_fit(loop_fits(y ~ s(x) + random(~ z | g), distrib = g1, data = d))
  expect_same_fit(loop_fits(y ~ s(x) | sigma ~ z, distrib = g1, data = d))
  expect_same_fit(loop_fits(y ~ s(x) + random(~ 1 | g), distrib = g1, data = d,
                            outer_criterion = ml()))
})

test_that("the compiled search reproduces the R route: t priors and a kinked block", {
  d <- loop_data()
  g1 <- distributions7::gaussian1_distrib()
  tp <- distributions7::fixed(distributions7::student_t1_distrib(), mu = 0)
  expect_same_fit(loop_fits(y ~ z + random(~ 1 | g, distrib = tp), distrib = g1, data = d))
  mv2 <- distributions7::fixed(distributions7::mvstudent_t1_distrib(2), mu1 = 0, mu2 = 0)
  expect_same_fit(loop_fits(y ~ z + random(~ z | g, distrib = mv2), distrib = g1, data = d))
  expect_same_fit(loop_fits(y ~ lasso(~ x + z) + random(~ 1 | g), distrib = g1, data = d))
})

test_that("the penalty built in C++ is the one R builds", {
  d <- loop_data()
  g1 <- distributions7::gaussian1_distrib()
  for (form in list(y ~ s(x) + random(~ z | g),
                    y ~ z + random(~ 1 | u | g) + random(~ 0 + z | u | g),
                    y ~ te(x, z))) {
    f <- statmod(form, distrib = g1, data = d)
    sp <- f@spec
    dz <- statmod_design(sp)
    idx <- outer_hyper_index(sp, statmod_blocks(sp, dz))
    beta <- unlist(f@coefficients[sp@distrib@params], use.names = FALSE)
    ps <- engine_pen_spec(sp, dz, idx, f@hyper, beta, "qr")
    expect_false(is.null(ps))
    v <- hyper_to_eta(f@hyper, idx)
    r <- engine_pen_build(ps, v, beta)
    expect_true(r$ok)
    cf <- statmod_objective(sp, f@hyper, dz, FALSE, NULL)$split(beta)
    em <- engine_model(sp, dz, cf, f@hyper, "qr")
    expect_identical(r$pen$S, em$pen$S, num.eq = FALSE)
    expect_identical(r$pen$C, em$pen$C, num.eq = FALSE)
    expect_identical(r$pen$units, em$pen$units, num.eq = FALSE)
  }
})
