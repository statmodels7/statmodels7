# A Jeffreys block beside a filter and a smoothed seg() in an equation that
# the filter does not drive. The penalty's derivatives are exact; the outer
# criterion still differences them, because its exact derivatives have not
# been checked against the criterion for this case (a fit of the panel used
# for that check took over half an hour). The derivatives are read near the
# starting point of the coefficients, with no fit; the negative control shows
# that the block's movement matters there.

test_that("the penalty beside a filter and a smoothed seg() in theta is exact", {
  set.seed(9)
  S <- 4L
  Tn <- 30L
  n <- S * Tn
  d <- data.frame(t = rep(1:Tn, S), id = factor(rep(1:S, each = Tn)),
                  xs = stats::runif(n, 0, 3), z = stats::rnorm(n))
  e <- numeric(n)
  for (i in 2:n) {
    e[i] <- if (d$t[i] == 1) 0 else 0.7 * e[i - 1] + stats::rnorm(1, sd = 0.5)
  }
  lev <- stats::rnorm(S, sd = 0.5)
  d$ys <- stats::rnbinom(n, mu = exp(1.5 + 0.4 * lev[d$id] + 0.5 * e),
                         size = exp(0.3 + 0.4 * d$z +
                                      1.5 * pmax(d$xs - 1.5, 0)))
  sp <- statmod_spec(ys ~ gas(p = 1, q = 1, omega ~ 1 + random(~ 1 | id),
                              time = t, by = id) |
                       theta ~ seg(xs, smoothed =
                                     numericals7::smooth_probit()) + jfp(~ z),
                     distributions7::negbin2_distrib(), d)
  des <- statmod_design(sp)
  expect_false(jfp_filter_outer_ok(sp, des, 1L))
  expect_false(jfp_filter_outer_ok(sp, des, 2L))
  hy <- statmod_hyper_start(sp, des)
  obj <- statmod_objective(sp, hy, des)
  jp <- statmod_joint_pieces(sp, des, obj, hy)
  nb <- jp$nb
  b0 <- statmod_start(sp, des, obj, NULL)
  u0 <- c(b0, jp$zeta())
  cf_at <- function(uu) {
    jp$setz(uu[nb + seq_along(jp$free)])
    obj$split(uu[seq_len(nb)])
  }
  set.seed(2)
  u <- u0 + 0.05 * stats::rnorm(length(u0))
  v <- stats::rnorm(length(u))
  w <- stats::rnorm(length(u))
  rel <- function(a, b) max(abs(a - b)) / max(abs(b))
  cd <- function(f, x, v, h) (f(x + h * v) - f(x - h * v)) / (2 * h)
  rich <- function(f, x, v, h) (4 * cd(f, x, v, h / 2) - cd(f, x, v, h)) / 3
  pen <- function(x, what) jfp_penalty(sp, des, cf_at(x), what, joint = TRUE)
  D3 <- function(x) jfp_penalty_move(sp, des, cf_at(x), v)
  errs <- function() {
    H <- pen(u, "hessian")
    c(rel(as.numeric(H %*% v), rich(function(x) pen(x, "gradient"), u, v,
                                     1e-4)),
      rel(D3(u), rich(function(x) pen(x, "hessian"), u, v, 1e-4)),
      rel(jfp_penalty_move(sp, des, cf_at(u), v, w = w),
          rich(D3, u, w, 1e-4)))
  }
  expect_lt(max(errs()), 1e-8)
  # negative control: without the moving block's parts
  ctl <- testthat::with_mocked_bindings(errs(),
    jfp_static_rows = function(st, dirs) lapply(st$params, function(p)
      matrix(0, st$spec@n_obs, st$joint$mf)),
    jfp_static_curv = function(st, g, dirs = list())
      matrix(0, st$joint$mf, st$joint$mf))
  expect_gt(min(ctl), 1e-6)
})
