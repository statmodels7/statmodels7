## A sharp break-point's objective has a kink in the break-point at the
## observations, so at its minimum every Gauss-Newton step is rejected
## whatever the Levenberg damping. iwls_fit() therefore does not escalate the
## damping on a rejected step there: measured on
## seg(x, psi ~ random(~1 | id)), escalating took the outer search from 6
## evaluations and 5.4 s of processor time to 28 and 25.1 s at an unchanged
## verdict, the one-sided slopes at the reported point reading +6.00 and
## -1.63 (piano_stabilita.txt 28d-28e).

test_that("a sharp break-point is recognized and a smoothed one is not", {
  set.seed(2)
  d <- data.frame(x = runif(80, 0, 10), y = rnorm(80))
  sp <- function(f) statmod_spec(f, distributions7::gaussian1_distrib(), d,
                                 NULL, NULL)
  expect_true(has_sharp_breakpoint(sp(y ~ seg(x))))
  expect_true(has_sharp_breakpoint(sp(y ~ jump(x))))
  expect_true(has_sharp_breakpoint(sp(y ~ jseg(x))))
  expect_false(has_sharp_breakpoint(
    sp(y ~ seg(x, smoothed = numericals7::smooth_probit()))))
  expect_false(has_sharp_breakpoint(sp(y ~ s(x, bspline_smooth(k = 8)))))
  expect_false(has_sharp_breakpoint(sp(y ~ x)))
})

test_that("with the damping off a rejected step ends the run at once", {
  ## the line search replaced by one that rejects at every length, which is
  ## what a kink at the minimum does; the calls are counted
  calls <- 0L
  local_mocked_bindings(iwls_line_search = function(...) {
    calls <<- calls + 1L
    list(ok = FALSE, sol = list(dropped = integer(0)), cand = NULL,
         vnew = NA_real_, step_used = 2^-30, errors = 0L)
  })
  obj <- list(fn = function(b) sum((b - 1)^2), gr = function(b) 2 * (b - 1))
  pieces_at <- function(b) list()
  on <- iwls_fit(obj, c(0, 0), iwls(), 10L, pieces_at)
  n_on <- calls
  calls <- 0L
  off <- iwls_fit(obj, c(0, 0), iwls(), 10L, pieces_at, damp_on_reject = FALSE)
  ## eight escalations and the attempt they retry, against one attempt
  expect_equal(n_on, 9L)
  expect_equal(calls, 1L)
  expect_match(on$note, "damping included", fixed = TRUE)
  expect_match(off$note, "no acceptable step", fixed = TRUE)
  expect_false(grepl("damping included", off$note, fixed = TRUE))
})
