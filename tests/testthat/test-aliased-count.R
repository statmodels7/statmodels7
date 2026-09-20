# The count of the active coordinates is their RANK only where the design
# identifies every column. Where it does not, the count overstates the degrees
# of freedom by the deficiency, and a coefficient the model does not identify
# is reported rather than left to be found in a summary nobody printed.
#
# The three things under test are separate and are tested separately: the
# certificate, the rank it falls back on, and what a reader is told.

collinear_data <- function(seed = 11L, nn = 600L) {
  set.seed(seed)
  g <- factor(sample(1:12, nn, TRUE))
  x <- stats::runif(nn, 0, 2)
  amp <- 2 + c(0.9, -0.7, 0.8, -0.6, rep(0, 8))[as.integer(g)]
  data.frame(y = amp * exp(-1.2 * x) + stats::rnorm(nn, 0, 0.15),
             x = x, g = g)
}

lasso_data <- function(seed = 5L, n = 400L, p = 12L) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  b <- c(2, -1.5, 1, rep(0, p - 3))
  d <- data.frame(y = drop(X %*% b) + stats::rnorm(n))
  d$X <- X
  d
}

test_that("the certificate answers for the design and not for one active set", {
  d <- lasso_data()
  spec <- statmod_spec(y ~ lasso(X), distributions7::gaussian1_distrib(), d)
  des <- statmod_design(spec)
  expect_true(design_count_exact(des))

  dc <- collinear_data()
  # a parameter's own intercept beside a complete set of indicators
  sp2 <- statmod_spec(y ~ 0 + nl(~ a * exp(-r * x), a ~ 1 + lasso(~ g)),
                      distributions7::gaussian1_distrib(), dc)
  de2 <- statmod_design(sp2)
  expect_false(design_count_exact(de2))
})

test_that("a block that moves is never certified, whatever its rank", {
  # the certificate cannot be memoized where a term recomputes its own block,
  # so it declines there and active_rank() answers instead
  dc <- collinear_data()
  sp <- statmod_spec(y ~ 0 + nl(~ a * exp(-r * x), a ~ 1 + lasso(~ g)),
                     distributions7::gaussian1_distrib(), dc)
  de <- statmod_design(sp)
  expect_true(length(attr(de, "refresh")) > 0L)
  expect_false(design_count_exact(de))
})

test_that("active_rank counts the independent active columns", {
  d <- lasso_data()
  spec <- statmod_spec(y ~ lasso(X), distributions7::gaussian1_distrib(), d)
  des <- statmod_design(spec)
  k <- des$mu$npar
  # a full-rank design: the rank of any subset is its size
  expect_identical(active_rank(des, seq_len(k)), as.numeric(k))
  expect_identical(active_rank(des, seq_len(4L)), 4)
  # and a stacked index reaching the second equation adds that equation's own
  tot <- k + des$sigma$npar
  expect_identical(active_rank(des, seq_len(tot)), as.numeric(tot))
})

test_that("block_column_rank agrees with the pivot the fit uses", {
  set.seed(1)
  X <- matrix(stats::rnorm(200 * 4), 200, 4)
  expect_identical(block_column_rank(X), qr(X)$rank)
  X[, 4L] <- X[, 1L] + X[, 2L]
  expect_identical(block_column_rank(X), 3L)
  # sparse and dense answer alike on the same matrix
  S <- Matrix::Matrix(X, sparse = TRUE)
  expect_identical(block_column_rank(S), 3L)
})

test_that("the count and the exact rank fit the same model", {
  # THE POINT OF THE CERTIFICATE IS EXACTNESS AND NOT A DIFFERENT FIT. On the
  # collinear model the deficiency is 1 or 2 at path points the selection does
  # not stop at, so the two routes are expected to agree; what would not agree
  # is the trace this replaced, which does not exist there.
  dc <- collinear_data()
  fm <- y ~ 0 + nl(~ a * exp(-r * x), a ~ 1 + lasso(~ g))
  f <- suppressWarnings(statmod(fm, distributions7::gaussian1_distrib(), dc))
  expect_true(is.finite(as.numeric(stats::logLik(f))))
  # the four groups that carry an amplitude survive and the eight that do not
  # are at zero. `coef()` and NOT the stored vector: @coefficients is unnamed
  # and carries the aliased coordinate where the pivot left it, where coef()
  # names the rows and reports that one as missing
  cf <- stats::coef(f)$mu
  keep <- grep("^nl\\.a\\.lasso", names(cf))
  expect_identical(sum(cf[keep] != 0, na.rm = TRUE), 4L)
})

test_that("the criterion reads the rank and not the count where they differ", {
  # THE DISCRIMINATING ASSERTION. Count and rank agree at every point the
  # selection stops at, so no comparison of two FITS can tell the branch from
  # the bare count; what separates them is tau at a point where the active set
  # is dependent, which is where the count overstates by the deficiency.
  dc <- collinear_data()
  fm <- y ~ 0 + nl(~ a * exp(-r * x), a ~ 1 + lasso(~ g))
  f <- suppressWarnings(statmod(fm, distributions7::gaussian1_distrib(), dc))
  spec <- f@spec
  des <- statmod_design_at(spec, f@coefficients, statmod_design(spec))
  act <- rep(TRUE, sum(vapply(des, function(z) z$npar, integer(1))))
  n_act <- sum(act)
  rk <- active_rank(des, which(act))
  # the premise: this model IS dependent at this point, or the test says
  # nothing at all
  expect_lt(rk, n_act)
  r <- statmod_pe(spec, des, f@coefficients, f@hyper, bic(), active = act)
  skip_if(is.null(r), "the penalized information has no Cholesky factor here")
  expect_identical(r$edf, rk)
  expect_false(isTRUE(all.equal(r$edf, n_act)))
})

test_that("a coefficient the model does not identify is reported", {
  dc <- collinear_data()
  fm <- y ~ 0 + nl(~ a * exp(-r * x), a ~ 1 + lasso(~ g))
  expect_warning(f <- statmod(fm, distributions7::gaussian1_distrib(), dc),
                 "not identified at the fitted point")
  expect_identical(length(f@aliased), 1L)
  # and the estimate is missing rather than zero, which is base R's convention.
  # @aliased names a coordinate of the whole model, so the equation's prefix
  # comes off before it indexes that equation's own vector
  bare <- sub("^[^:]+:", "", f@aliased)
  expect_true(is.na(stats::coef(f)$mu[[bare]]))
})

test_that("nothing is reported where every column is identified", {
  # WHAT IS ASSERTED IS THAT NO ALIASING WARNING IS RAISED, not that the fit is
  # silent: a path that reaches the end of its grid says so, which is another
  # warning and another subject
  aliasing_warnings <- function(expr) {
    out <- character(0)
    withCallingHandlers(force(expr), warning = function(cond) {
      if (grepl("not identified at the fitted point", conditionMessage(cond))) {
        out <<- c(out, conditionMessage(cond))
      }
      invokeRestart("muffleWarning")
    })
    out
  }
  d <- lasso_data()
  f <- NULL
  expect_identical(
    aliasing_warnings(f <- statmod(y ~ lasso(X),
                                   distributions7::gaussian1_distrib(), d)),
    character(0))
  expect_identical(length(f@aliased), 0L)
  # the control that keeps the warning from being vacuous: a random effect is
  # deficient in the DESIGN and identified by its own penalty, and a rank test
  # read there rather than on K would name it
  set.seed(7)
  d$g <- factor(sample(1:20, nrow(d), TRUE))
  d$y <- d$y + stats::rnorm(20)[d$g]
  f2 <- NULL
  expect_identical(
    aliasing_warnings(f2 <- statmod(y ~ lasso(X) + random(~ 1 | g),
                                    distributions7::gaussian1_distrib(), d)),
    character(0))
  expect_identical(length(f2@aliased), 0L)
})

test_that("warn_aliased says nothing about an empty set", {
  expect_silent(warn_aliased(character(0)))
  expect_warning(warn_aliased("mu:x1"), "1 coefficient is not identified")
  expect_warning(warn_aliased(c("mu:x1", "mu:x2")),
                 "2 coefficients are not identified")
  # a long list is not printed to the end
  w <- tryCatch(warn_aliased(paste0("mu:x", 1:9)),
                warning = function(cond) conditionMessage(cond))
  expect_true(grepl("and 5 more", w, fixed = TRUE))
})

test_that("a summary prints where a selected column is not identified", {
  # AN ALIASED COEFFICIENT IS NOT A SELECTED ONE. Its estimate is NA, so
  # `!= 0` is NA, and a logical NA indexed a row of nothing but NAs: the
  # label column carried one, the block's width came out NA and the printer
  # died inside formatC(). Reachable from any fit whose selected columns are
  # dependent.
  dc <- collinear_data()
  fm <- y ~ 0 + nl(~ a * exp(-r * x), a ~ 1 + lasso(~ g))
  f <- suppressWarnings(statmod(fm, distributions7::gaussian1_distrib(), dc))
  txt <- utils::capture.output(print(summary(f)))
  expect_true(any(grepl("lasso.g1", txt, fixed = TRUE)))
  # the aliased one is not among the rows reported as selected
  expect_false(any(grepl("lasso.g12", txt, fixed = TRUE)))
})
