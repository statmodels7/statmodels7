# The outer search in compiled code (engine_outer_search(): the search, the
# point and the penalty, phase 4) against the R route of outer_fit(): the
# same fit in every property (numbers to rounding), and the penalty built in
# C++ identical to the one R builds, signed zeros included.

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

# the same fit up to rounding: numeric leaves agree to 1e-12 on the scale
# of the R route's values, floor 1 (on Linux the two routes have differed
# by one ulp in the objective, run to run, and a block's change is a
# difference of two objectives near 1e-14); every other leaf, the structure
# and the attributes exactly
same_fit_value <- function(a, b, tol = 1e-12) {
  if (same_value(a, b)) return(TRUE)
  if (is.numeric(a) && is.numeric(b)) return(scale_gap(a, b, floor = 1) <= tol)
  if (is.list(a) && is.list(b) && !is.object(a) && !is.object(b)) {
    if (length(a) != length(b) || !identical(names(a), names(b))) return(FALSE)
    return(all(vapply(seq_along(a), function(k) same_fit_value(a[[k]], b[[k]], tol), TRUE)))
  }
  if (is.data.frame(a) && is.data.frame(b)) {
    if (!identical(dim(a), dim(b)) || !identical(names(a), names(b))) return(FALSE)
    return(all(vapply(seq_along(a), function(k) same_fit_value(a[[k]], b[[k]], tol), TRUE)))
  }
  FALSE
}

# On arm64 macOS R's own BLAS rounds with fused multiply-adds and the
# compiled code does not, so a search along a weakly identified direction
# can drift by 1e-9 and reach the rounding floor on one route only: there
# the fits marked `at_resolution` are compared at the search's resolution,
# and the stopping rule of one route may carry the floor rule besides
on_arm_mac <- identical(Sys.info()[["sysname"]], "Darwin")

expect_same_fit <- function(f, at_resolution = FALSE) {
  tol <- if (at_resolution && on_arm_mac) 1e-6 else 1e-12
  for (p in setdiff(S7::prop_names(f$r), c("elapsed", "methods"))) {
    expect_true(same_fit_value(S7::prop(f$c, p), S7::prop(f$r, p), tol), info = p)
  }
  # the optimizer's line search holds a closure over each run's frame: its
  # class and its rule are compared, the rest of the record exactly
  mr <- f$r@methods
  mc <- f$c@methods
  for (k in setdiff(names(mr), "search")) expect_true(same_fit_value(mc[[k]], mr[[k]], tol), info = k)
  if (!is.null(mr$search)) {
    expect_identical(class(mc$search), class(mr$search))
    lc <- mc$search@criterion@label
    lr <- mr$search@criterion@label
    if (at_resolution && on_arm_mac) {
      expect_true(startsWith(lc, lr) || startsWith(lr, lc), info = paste(lc, "|", lr))
    } else {
      expect_identical(lc, lr)
    }
  }
}

test_that("the compiled search reproduces the R route: smooth, random slope, two equations", {
  d <- loop_data()
  g1 <- distributions7::gaussian1_distrib()
  expect_same_fit(loop_fits(y ~ s(x), distrib = g1, data = d))
  expect_same_fit(loop_fits(y ~ s(x) + random(~ z | g), distrib = g1, data = d))
  expect_same_fit(loop_fits(y ~ s(x) | sigma ~ z, distrib = g1, data = d))
  expect_same_fit(loop_fits(y ~ s(x) + random(~ 1 | g), distrib = g1, data = d,
                            outer_criterion = ml()), at_resolution = TRUE)
})

test_that("the compiled search reproduces the R route: t priors and a kinked block", {
  d <- loop_data()
  g1 <- distributions7::gaussian1_distrib()
  tp <- distributions7::fixed(distributions7::student_t1_distrib(), mu = 0)
  expect_same_fit(loop_fits(y ~ z + random(~ 1 | g, distrib = tp), distrib = g1, data = d))
  mv2 <- distributions7::fixed(distributions7::mvstudent_t1_distrib(2), mu1 = 0, mu2 = 0)
  expect_same_fit(loop_fits(y ~ z + random(~ z | g, distrib = mv2), distrib = g1, data = d),
                  at_resolution = TRUE)
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
