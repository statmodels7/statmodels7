# A sharp break-point developed over groups can have its penalized mode with
# one group's break-point exactly on an observation. The scoring step there
# crosses the kink and is rejected at every length, and the inner fit used to
# stop at that point with every other coordinate wherever it happened to be.

kink_panel <- function() {
  set.seed(31)
  m <- 12
  ni <- 30
  id <- factor(rep(seq_len(m), each = ni))
  x <- as.numeric(replicate(m, sort(runif(ni, 0, 10))))
  psi <- rnorm(m, 5, 0.4)
  mu <- 1 + 0.5 * x + 1.8 * pmax(x - psi[as.integer(id)], 0)
  data.frame(y = mu + rnorm(m * ni, 0, 0.4), x = x, id = id)
}

kink_setup <- function(sig) {
  d <- kink_panel()
  f <- y ~ seg(x, psi ~ random(~1 | id))
  spec <- statmod_spec(f, gaussian1_distrib(), d, NULL, NULL)
  inner <- iwls_resolve(iwls(maxit = 10000), spec@distrib)
  cfg <- inner_settings(inner)
  hy <- statmod_hyper_start(spec, statmod_design(spec))
  key <- names(hy$mu)[1L]
  hy$mu[[key]][["sigma"]] <- sig
  run <- function(start) {
    design <- statmod_design(spec)
    blocks <- statmod_blocks(spec, design)
    obj <- statmod_objective(spec, hy, design, cfg$expected, cfg$approx)
    if (is.null(start)) start <- statmod_start(spec, design, obj, NULL)
    r <- statmod_alternate(spec, design, blocks, hy, inner, start,
                           cfg$expected, cfg$approx, cfg$maxit, cfg$tol,
                           verbosity(0), hold_refresh = TRUE)
    kk <- kink_positions(spec, design, r$obj$split(r$par),
                         r$obj$split(seq_along(r$par)))
    list(par = r$par, value = r$value, kinks = kk, gr = r$obj$gr(r$par))
  }
  list(run = run)
}

test_that("the inner fit reaches the mode from a point stopped at a kink", {
  skip_on_cran()
  s <- kink_setup(0.4483794)
  m <- s$run(NULL)
  # the premise: at this smoothing the mode puts a break-point on an
  # observation. Where it does not on some platform there is nothing to test.
  skip_if(!length(m$kinks), "the mode has no break-point on an observation here")
  expect_true(all(m$kinks %in% 4:16))
  # the other coordinates are at their optimum there
  expect_lt(max(abs(m$gr[-m$kinks])), 1e-2)

  # move every OTHER group's deviation off its optimum, keeping the kinked
  # break-point where it is: a step on all the coordinates crosses the kink
  start <- m$par
  dev <- setdiff(5:16, m$kinks)
  start[dev] <- start[dev] + 0.05
  r <- s$run(start)
  expect_lt(r$value - m$value, 1e-6)
  expect_lt(max(abs(r$gr[-r$kinks])), 1e-2)
})

test_that("seg with a developed break-point reads one criterion across routes", {
  skip_on_cran()
  d <- kink_panel()
  f <- y ~ seg(x, psi ~ random(~1 | id))
  a <- statmod(f, gaussian1_distrib(), data = d,
               inner_optimizer = iwls(maxit = 10000))
  b <- statmod(f, gaussian1_distrib(), data = d,
               inner_optimizer = iwls(maxit = 10000),
               outer_optimizer = lbfgs())
  expect_lt(abs(a@criterion - b@criterion), 1e-3)
})
