# The depth-one exact memo on the structural state: a hit returns the
# stored object itself, a design carrying refreshable terms bypasses it,
# and a fit that reads the same point twice reads one computation.

test_that("structural_memo hits on an identical key and misses otherwise", {
  d <- structure(list(), class = "list")
  sst <- new.env(parent = emptyenv())
  attr(d, "structure") <- sst
  calls <- 0L
  f <- function() { calls <<- calls + 1L; list(v = calls) }
  a <- structural_memo(d, "x", list(k = 1), f)
  b <- structural_memo(d, "x", list(k = 1), f)
  expect_identical(calls, 1L)
  # the stored object itself, not a copy of equal value
  expect_identical(a, b)
  c2 <- structural_memo(d, "x", list(k = 2), f)
  expect_identical(calls, 2L)
  expect_identical(c2$v, 2L)
  # slots do not collide
  structural_memo(d, "y", list(k = 1), f)
  expect_identical(calls, 3L)
})

test_that("a design with no structural state or with refreshable terms bypasses", {
  d <- structure(list(), class = "list")
  calls <- 0L
  f <- function() { calls <<- calls + 1L; calls }
  structural_memo(d, "x", list(k = 1), f)
  structural_memo(d, "x", list(k = 1), f)
  expect_identical(calls, 2L)

  sst <- new.env(parent = emptyenv())
  attr(d, "structure") <- sst
  attr(d, "refresh") <- list(something = TRUE)
  structural_memo(d, "x", list(k = 1), f)
  structural_memo(d, "x", list(k = 1), f)
  expect_identical(calls, 4L)
})

test_that("a panel fit through the memo matches its own criterion identities", {
  # small end-to-end control: the memo must not move a structural fit at all
  set.seed(7)
  groups <- 6L; per <- 30L
  id <- factor(rep(seq_len(groups), each = per))
  n <- length(id)
  x <- rnorm(n)
  omega_i <- rnorm(groups, 0.2, 0.3)
  y <- numeric(n)
  for (g in seq_len(groups)) {
    rows <- which(as.integer(id) == g)
    f <- omega_i[g] / (1 - 0.5)
    s <- 0
    for (k in seq_along(rows)) {
      f <- omega_i[g] + 0.2 * s + 0.5 * f
      eta <- 0.5 * x[rows[k]] + f
      y[rows[k]] <- eta + rnorm(1)
      s <- y[rows[k]] - eta
    }
  }
  d <- data.frame(id = id, t = rep(seq_len(per), groups), x = x, y = y)
  fit <- statmod(y ~ x + gas(p = 1, q = 1, omega ~ random(~1 | id),
                             by = id, time = t),
                 distributions7::gaussian1_distrib(), d,
                 outer_criterion = reml(hessian = "observed"))
  expect_true(is.finite(fit@criterion))
  # vcov reads the full information at the fitted point, which the memo
  # serves: it must be symmetric and finite
  vc <- vcov(fit)
  expect_true(all(is.finite(diag(vc))))
  expect_identical(vc, t(vc))
})
