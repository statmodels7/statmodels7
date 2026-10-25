# A break-point that jump() holds at the minimum of its profile has a zero row
# and column in the joint information. Beside a score-driven filter, the
# outer gradient's joint movement and the structural table both inverted
# that matrix whole, so the certificate read "the outer gradient is not
# finite" and the filter's loading and persistence lost their standard
# errors (pig2 with jump(tt) + gas(1, 1) | alpha ~ s(tt) on 800 counts; the
# same on the gaussian below).

sim_jump_filter <- function(seed = 1, n = 200) {
  set.seed(seed)
  d <- data.frame(t = seq_len(n))
  d$tt <- d$t / n
  f <- 0
  y <- numeric(n)
  for (i in seq_len(n)) {
    m <- 3 + 0.5 * (i > 0.6 * n) + f
    y[i] <- stats::rnorm(1, m, 0.3)
    f <- 0.1 * (y[i] - m) + 0.8 * f
  }
  d$y <- y
  d
}

test_that("a held break-point beside a filter keeps the outer gradient", {
  skip_on_cran()
  d <- sim_jump_filter()
  fit <- statmod(y ~ jump(tt) + gas(p = 1, q = 1, time = t) | sigma ~ s(tt),
                 distrib = distributions7::gaussian1_distrib(), data = d)
  # the smooth on sigma goes to its edge here (lambda about 2e8), which the
  # certificate reads as a boundary once the gradient exists
  ce <- statmod_certificate(fit)
  expect_true(ce$state %in% c("converged", "boundary"))
  expect_false(identical(ce$reason,
                         "the outer gradient is not finite at the reported point"))
  st <- statmod_structural_table(fit)
  expect_true(all(is.finite(st$se[!st$held])))
})
