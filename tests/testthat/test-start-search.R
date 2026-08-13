# A global search for the starting values, run once before the fit.
#
# The first two tests are the ones that decide whether the feature is safe:
# a convex fit must be UNCHANGED by it, and the coordinates it searches must
# be the ones where the problem is not convex. What it buys where the default
# already works is nothing, which is the measured answer and not an omission.

set.seed(12)
n <- 200
dd <- data.frame(x = runif(n), z = runif(n))
dd$y <- 1 + 2 * dd$x + stats::rnorm(n, sd = 0.4)

test_that("a convex fit is not moved by searching first", {
  # If it were, the fit had not converged: on a convex problem the scoring
  # step reaches the same optimum from anywhere, so this is a statement about
  # the fit as much as about the search.
  ref <- statmod(y ~ x + z, distributions7::gaussian1_distrib(), dd)
  set.seed(1)
  got <- statmod(y ~ x + z, distributions7::gaussian1_distrib(), dd,
                 start = start_search(optimizers7::sa(maxit = 15)))
  expect_true(got@converged)
  expect_equal(got@coefficients, ref@coefficients, tolerance = 1e-8)
  expect_equal(got@objective, ref@objective, tolerance = 1e-10)
})

test_that("it searches the non-convex coordinates and leaves the rest", {
  set.seed(3)
  d1 <- data.frame(x = sort(stats::runif(300, 0, 10)),
                   g = factor(rep(1:10, 30)))
  u <- stats::rnorm(10, sd = 0.4)
  d1$y <- 1 + 0.4 * d1$x + 3 * pmax(d1$x - 6, 0) + u[as.integer(d1$g)] +
    stats::rnorm(300, sd = 0.4)
  sp <- statmod_spec(y ~ seg(x) + random(~ 1 | g),
                     distributions7::gaussian1_distrib(), d1)
  de <- statmod_design(sp)
  idx <- statmodels7:::search_coords(sp, de)
  nms <- unlist(lapply(sp@distrib@params, function(p) de[[p]]$coef_names),
                use.names = FALSE)
  picked <- nms[idx]
  # the break-point term, which recomputes its own block, and the intercepts
  expect_true(all(grepl("^seg", picked) | grepl("Intercept", picked)))
  # and NOT the ten random-effect coefficients, whose block is convex: a
  # search over them would spend the budget where it buys nothing
  expect_false(any(grepl("^random", picked)))
  expect_length(picked, 5L)
})

test_that("the search reaches a starting point a bad one cannot", {
  # Against start_origin(), which on a response centred far from zero is a
  # genuinely bad place to begin: this is the mechanism working, not a claim
  # that the DEFAULT needs it.
  set.seed(21)
  ds <- data.frame(x = stats::runif(150))
  ds$y <- 500 + 40 * ds$x + stats::rnorm(150, sd = 5)
  ref <- statmod(y ~ x, distributions7::gaussian1_distrib(), ds)
  set.seed(2)
  z <- statmod(y ~ x, distributions7::gaussian1_distrib(), ds,
               start = start_origin(), inner_optimizer = iwls(maxit = 8L))
  set.seed(2)
  s <- statmod(y ~ x, distributions7::gaussian1_distrib(), ds,
               start = start_search(optimizers7::sa(maxit = 25)),
               inner_optimizer = iwls(maxit = 8L))
  # the searched start is nearer the answer than the zero one, on a budget
  # too short to recover from a bad beginning
  expect_lt(abs(s@objective - ref@objective), abs(z@objective - ref@objective))
})

test_that("over= names the coordinates and refuses what is not one", {
  set.seed(1)
  a <- statmod(y ~ x + z, distributions7::gaussian1_distrib(), dd,
               start = start_search(optimizers7::sa(maxit = 10),
                                    over = c("x", "z")))
  expect_true(a@converged)
  expect_error(statmod(y ~ x, distributions7::gaussian1_distrib(), dd,
                       start = start_search(over = "nope")),
               "not a coefficient")
})

test_that("a filter's own parameters are reachable by the search", {
  # They do not live in the coefficient vector -- the objective reads them
  # from the design's structural state -- so this is the one place a strategy
  # returning coefficients still has to reach past them.
  set.seed(5)
  nt <- 300
  f <- numeric(nt); s <- numeric(nt)
  fl <- 0.5 / (1 - 0.6)
  for (t in seq_len(nt)) {
    fl <- 0.5 + 0.3 * (if (t > 1) s[t - 1] else 0) +
      0.6 * (if (t > 1) f[t - 1] else fl)
    f[t] <- fl
    s[t] <- stats::rnorm(1, sd = 0.5)
  }
  d2 <- data.frame(t = seq_len(nt), y = f + s)
  fml <- y ~ gas(p = 1, q = 1, time = t) - 1
  ref <- statmod(fml, distributions7::gaussian1_distrib(), d2,
                 outer_criterion = NULL)
  set.seed(6)
  got <- statmod(fml, distributions7::gaussian1_distrib(), d2,
                 outer_criterion = NULL,
                 start = start_search(optimizers7::sa(maxit = 15)))
  expect_true(got@converged)
  # the same optimum: the fit is robust here, and what this asserts is that
  # searching over the term's parameters did not damage it
  expect_equal(got@objective, ref@objective, tolerance = 1e-6)
  expect_equal(unlist(got@structural), unlist(ref@structural),
               tolerance = 1e-4)
})

test_that("the constructor validates its arguments", {
  expect_error(start_search(optimizer = "sa"), "must be an optimizers7")
  expect_error(start_search(over = 1), "character vector")
  expect_error(start_search(hyper = NA), "TRUE or FALSE")
  expect_output(print(start_search()), "search")
})
