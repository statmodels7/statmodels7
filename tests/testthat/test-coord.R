# The compiled coordinate descent, and the routes it has to agree with.

set.seed(11)
nc <- 200L
pc <- 20L
xc <- matrix(stats::rnorm(nc * pc), nc, pc)
dc <- data.frame(y = as.numeric(xc %*% c(2, -1.5, 1, rep(0, pc - 3))) +
                   stats::rnorm(nc))
dc$x <- xc

# the kernel written out in R, which the compiled one is held to
coord_descent_r <- function(X, z, w, beta, cut, slope, icept, maxit, tol) {
  p <- ncol(X)
  v <- as.numeric(crossprod(w, X^2))
  r <- z - as.numeric(X %*% beta)
  apply1 <- function(u, j) {
    a <- abs(u)
    k <- which(a <= cut[j, ])[1L]
    if (is.na(k)) k <- ncol(cut)
    sign(u) * (slope[j, k] * a + icept[j, k])
  }
  for (it in seq_len(maxit)) {
    delta <- 0
    for (j in seq_len(p)) {
      if (v[j] <= 0) next
      u <- beta[j] + sum(w * X[, j] * r) / v[j]
      nb <- apply1(u, j)
      if (nb != beta[j]) {
        r <- r - X[, j] * (nb - beta[j])
        delta <- max(delta, abs(nb - beta[j]))
        beta[j] <- nb
      }
    }
    if (delta < tol) break
  }
  list(beta = beta, sweeps = it)
}

setup_block <- function(f, nm, th, data = dc) {
  spec <- statmod_spec(f, distributions7::gaussian1_distrib(), data)
  design <- statmod_design(spec)
  blocks <- statmod_blocks(spec, design)
  hy <- statmod_hyper_merge(spec, statmod_hyper_start(spec),
                            stats::setNames(list(stats::setNames(list(th), nm)),
                                            "mu"))
  obj <- statmod_objective(spec, hy, design)
  list(spec = spec, design = design, block = blocks$sparse[[1L]], hyper = hy,
       obj = obj, beta = statmod_start(spec, design, obj, NULL))
}

test_that("the compiled kernel is the R twin", {
  # the arithmetic is the same in the same order, but a compiler may contract
  # a multiply-add, so the two are compared at a tolerance and not by identity
  b <- setup_block(y ~ lasso(x), "lasso(x)", c(lambda = 8))
  cols <- b$design$mu$blocks[["lasso(x)"]]
  X <- b$design$mu$X[, cols, drop = FALSE]
  w <- rep(1.3, nc)
  z <- dc$y
  v <- as.numeric(crossprod(w, X^2))
  tab <- penalties7::penalty_prox_spec(b$block$penalty, list(lambda = 8), 1 / v)
  b0 <- numeric(ncol(X))
  all_cols <- seq_len(ncol(X)) - 1L
  r <- coord_descent_r(X, z, w, b0, tab$cut, tab$slope, tab$icept, 500L, 1e-12)
  # and the two ways of holding the gradient are the same arithmetic seen
  # from two sides, so they have to agree with each other as well
  for (cov in c(FALSE, TRUE)) {
    a <- coord_descent(X, z, w, b0, tab$cut, tab$slope, tab$icept, all_cols,
                       500L, 1e-12, cov)
    expect_equal(a$beta, r$beta, tolerance = 1e-10,
                 label = if (cov) "covariance" else "naive")
  }
  a <- coord_descent(X, z, w, b0, tab$cut, tab$slope, tab$icept, all_cols,
                     500L, 1e-12, FALSE)
  # the sweep counts are NOT comparable and asserting they were is what this
  # line used to do: the compiled kernel cycles over the active set between
  # full sweeps, so it takes more passes and each one is cheaper. What the
  # twin establishes is the arithmetic, not the schedule.
  expect_gt(a$sweeps, 0L)
  expect_gt(r$sweeps, 0L)
})

test_that("coordinate descent and the proximal route reach the same point", {
  # they share the objective and nothing else: one reads the block's columns
  # and the running residual, the other the whole model through fn and gr
  cases <- list(
    list(f = y ~ lasso(x), nm = "lasso(x)", th = c(lambda = 8)),
    list(f = y ~ enet(x), nm = "enet(x)", th = c(lambda = 8, alpha = 0.6)),
    list(f = y ~ scad(x), nm = "scad(x)", th = c(lambda = 8, a = 3.7)),
    list(f = y ~ mcp(x), nm = "mcp(x)", th = c(lambda = 8, gamma = 3)))
  for (cs in cases) {
    b <- setup_block(cs$f, cs$nm, cs$th)
    cd <- sparse_fit(b$obj, b$beta, b$block, b$hyper, spec = b$spec,
                     design = b$design)
    px <- sparse_fit(b$obj, b$beta, b$block, b$hyper)
    expect_false(is.null(cd), label = cs$nm)
    expect_equal(cd$par, px$par, tolerance = 1e-7, label = cs$nm)
    expect_equal(cd$value, px$value, tolerance = 1e-9, label = cs$nm)
  }
})

test_that("a response that is not gaussian rebuilds the working quadratic", {
  # the quadratic is exact for a gaussian with an identity link and one pass
  # is the answer; elsewhere the weights move and the loop repeats
  set.seed(4)
  dp <- data.frame(cnt = stats::rpois(300, exp(0.5)))
  dp$x <- matrix(stats::rnorm(300 * 10), 300, 10)
  spec <- statmod_spec(cnt ~ lasso(x), distributions7::poisson_distrib(), dp)
  design <- statmod_design(spec)
  blocks <- statmod_blocks(spec, design)
  hy <- statmod_hyper_merge(spec, statmod_hyper_start(spec),
                            list(mu = list("lasso(x)" = c(lambda = 3))))
  obj <- statmod_objective(spec, hy, design)
  b0 <- statmod_start(spec, design, obj, NULL)
  cd <- sparse_fit(obj, b0, blocks$sparse[[1L]], hy, spec = spec,
                   design = design)
  px <- sparse_fit(obj, b0, blocks$sparse[[1L]], hy)
  # Looser on the coefficients than the gaussian case above, and the reason is
  # in the subject of the test: here the working quadratic is rebuilt, so the
  # two routes' difference compounds over the passes rather than being bounded
  # by one solve. macOS disagreed at 3e-6 relative on the smallest
  # coefficients where the other four platforms were inside 1e-7. The
  # objective, which both routes minimize, is held to 1e-9, and a rebuilt
  # quadratic that was actually wrong moves the coefficients by percents.
  expect_equal(cd$par, px$par, tolerance = 1e-5)
  expect_equal(cd$value, px$value, tolerance = 1e-9)
})

test_that("a penalty with no table keeps the proximal route", {
  # the route is chosen by what the penalty can describe, not by its name
  heavy <- penalties7::heavy_penalty(n_coef = 5L)
  expect_null(penalties7::penalty_prox_spec(heavy, list(sigma = 1, nu = 4),
                                            rep(0.5, 5)))
  b <- setup_block(y ~ lasso(x), "lasso(x)", c(lambda = 8))
  # the block's own penalty does have one, so the fast route applies
  expect_false(is.null(coord_fit(b$obj, b$beta, b$block, b$hyper, b$spec,
                                 b$design, TRUE, "bartlett")))
})

test_that("the lasso agrees with glmnet on the objective they share", {
  # the lasso is homogeneous of degree one in lambda, so a rescaling of the
  # loss is absorbed exactly and the two objectives are the same problem.
  # Measured: 0.00e+00, 1.4e-14 and 1.1e-13 at three sizes. That does NOT
  # hold for the elastic net, SCAD or MCP, whose shape parameters make the
  # penalty non-homogeneous, so only the timing is comparable there.
  skip_if_not_installed("glmnet")
  n <- 200L
  set.seed(7)
  X <- scale(matrix(stats::rnorm(n * 20L), n, 20L), TRUE, FALSE)
  y <- as.numeric(X %*% c(2, -1.5, 1, rep(0, 17))) + stats::rnorm(n)
  y <- y - mean(y)
  dd <- data.frame(y = y)
  dd$x <- X
  lam <- 0.05 * n
  b <- setup_block(y ~ lasso(x) - 1 | sigma ~ 1, "lasso(x)",
                   c(lambda = lam), dd)
  cd <- sparse_fit(b$obj, b$beta, b$block, b$hyper, spec = b$spec,
                   design = b$design)
  cf <- b$obj$split(cd$par)
  s <- exp(cf$sigma[1L])
  g <- glmnet::glmnet(X, y, alpha = 1, lambda = lam * s^2 / n,
                      standardize = FALSE, intercept = FALSE)
  ours <- function(bb) sum((y - X %*% bb)^2) / (2 * s^2) +
    penalties7::penalty_value(b$block$penalty, bb, list(lambda = lam))
  expect_lt(abs(ours(cf$mu) - ours(as.numeric(g$beta))), 1e-8)
  expect_identical(abs(cf$mu) > 1e-8, abs(as.numeric(g$beta)) > 1e-8)
})

test_that("a whole fit goes through the compiled route and still converges", {
  fit <- statmod(y ~ lasso(x), distributions7::gaussian1_distrib(), dc,
                 hyper = list(mu = list("lasso(x)" = c(lambda = 8))))
  expect_true(fit@converged)
  expect_lt(sum(abs(fit@coefficients$mu[-1L]) > 1e-8), pc)
  # and the selection still answers the smoothing parameter
  hi <- statmod(y ~ lasso(x), distributions7::gaussian1_distrib(), dc,
                hyper = list(mu = list("lasso(x)" = c(lambda = 200))))
  expect_lt(sum(abs(hi@coefficients$mu[-1L]) > 1e-8),
            sum(abs(fit@coefficients$mu[-1L]) > 1e-8))
})


test_that("screening is checked, so a rule that discards too much is exact", {
  # The sequential strong rule assumes the gradient moves no faster than the
  # threshold, which is not a theorem: it can discard a coordinate that
  # belongs in the fit. What makes the answer exact is reading the gradient
  # over EVERY column at the point reached and putting back whatever exceeds
  # the kink. Screening the block down to one column and letting the check
  # rebuild it is that mechanism run at its worst case.
  b <- setup_block(y ~ lasso(x), "lasso(x)", c(lambda = 8))
  cols <- b$design$mu$blocks[["lasso(x)"]]
  X <- b$design$mu$X[, cols, drop = FALSE]
  w <- rep(1, nc)
  z <- dc$y
  v <- as.numeric(crossprod(w, X^2))
  th <- list(lambda = 8)
  s <- kink_scale(b$block$penalty, th)

  answer <- NULL
  for (keep0 in list(seq_len(ncol(X)), 1L, c(1L, 5L))) {
    keep <- keep0
    b0 <- numeric(ncol(X))
    repeat {
      tab <- penalties7::penalty_prox_spec(b$block$penalty, th, 1 / v[keep])
      out <- coord_descent(X, z, w, b0, tab$cut, tab$slope, tab$icept,
                           as.integer(keep - 1L), 500L, 1e-12, FALSE)
      back <- setdiff(which(abs(out$grad) > s * (1 + 1e-10)), keep)
      if (!length(back)) break
      keep <- sort(c(keep, back))
    }
    if (is.null(answer)) answer <- out$beta else
      expect_equal(out$beta, answer, tolerance = 1e-9)
  }
  # the check is what does it: without putting anything back, one column
  # gives a different answer
  tab <- penalties7::penalty_prox_spec(b$block$penalty, th, 1 / v[1L])
  bad <- coord_descent(X, z, w, numeric(ncol(X)), tab$cut, tab$slope,
                       tab$icept, 0L, 500L, 1e-12, FALSE)
  expect_gt(max(abs(bad$beta - answer)), 1e-3)
})

test_that("the screening rule keeps every coordinate that is alive", {
  b <- setup_block(y ~ lasso(x), "lasso(x)", c(lambda = 8))
  cols <- b$design$mu$blocks[["lasso(x)"]]
  X <- b$design$mu$X[, cols, drop = FALSE]
  w <- rep(1, nc)
  s <- kink_scale(b$block$penalty, list(lambda = 8))
  # a coordinate away from zero is never screened out, whatever its gradient
  beta <- numeric(ncol(X))
  beta[7L] <- 0.4
  expect_true(7L %in% coord_screen(X, w, dc$y, beta, s, 1e6))
  # with no previous point the rule falls back to its global form
  expect_gt(length(coord_screen(X, w, dc$y, numeric(ncol(X)), s, NULL)), 0L)
  # and a kink so wide that nothing survives still leaves one coordinate to
  # visit rather than an empty sweep
  expect_identical(length(coord_screen(X, w, dc$y, numeric(ncol(X)), 1e8,
                                       1e8)), 1L)
})


test_that("the path carries the previous kink onto the blocks", {
  # the previous point travels on the blocks rather than through the argument
  # list of every layer between the path and the descent
  b <- setup_block(y ~ lasso(x), "lasso(x)", c(lambda = 8))
  expect_null(b$block$prev_kink)
  bk <- blocks_at_kink(list(sparse = list(b$block)), b$hyper)
  expect_equal(bk$sparse[[1L]]$prev_kink, 8, tolerance = 1e-9)
  # a kink of 1e-300 is still a kink, and the first version of this asserted
  # otherwise; what leaves nothing to screen against is a penalty with no kink
  hy0 <- b$hyper
  hy0$mu[["lasso(x)"]][["lambda"]] <- 1e-300
  expect_equal(blocks_at_kink(list(sparse = list(b$block)),
                              hy0)$sparse[[1L]]$prev_kink, 1e-300)
  spec2 <- statmod_spec(y ~ ridge(x), distributions7::gaussian1_distrib(), dc)
  pen <- modelterms7::term_penalty(spec2@terms$mu[["ridge(x)"]])
  expect_equal(kink_scale(pen, list(sigma = 1)), 0)
})

test_that("screening along a path does not change where the path lands", {
  # the rule discards, the check puts back, and the answer is the one the
  # unscreened path reaches. That is the property worth pinning: the speed is
  # a measurement and this is a fact.
  set.seed(21)
  n <- 150L
  p <- 60L
  X <- scale(matrix(stats::rnorm(n * p), n, p), TRUE, FALSE)
  y <- as.numeric(X %*% c(rep(2, 4), rep(0, p - 4))) + stats::rnorm(n)
  dd <- data.frame(y = y - mean(y))
  dd$x <- X
  spec <- statmod_spec(y ~ lasso(x) - 1 | sigma ~ 1,
                       distributions7::gaussian1_distrib(), dd)
  design <- statmod_design(spec)
  blocks <- statmod_blocks(spec, design)
  hy <- statmod_hyper_start(spec)
  obj <- statmod_objective(spec, hy, design)
  b0 <- statmod_start(spec, design, obj, NULL)
  top <- path_null_score(obj, b0, blocks$sparse[[1L]], hy)
  vals <- path_values(blocks$sparse[[1L]]$penalty, hy$mu[[1L]], "lambda", top,
                      12L, 1e-2)

  walk <- function(screen) {
    warm <- b0
    last <- NULL
    out <- list()
    for (v in vals) {
      h <- hyper_set(hy, list(parameter = "mu", term = "lasso(x)",
                              name = "lambda"), v)
      bk <- if (screen && !is.null(last)) blocks_at_kink(blocks, last) else
        blocks
      r <- statmod_alternate(spec, design, bk, h, iwls(), warm, TRUE,
                             "bartlett", 100L, 1e-8, verbosity(0))
      warm <- r$par
      last <- h
      out[[length(out) + 1L]] <- r$obj$split(r$par)$mu
    }
    out
  }
  a <- walk(FALSE)
  b <- walk(TRUE)
  for (j in seq_along(a)) {
    expect_equal(b[[j]], a[[j]], tolerance = 1e-7,
                 label = sprintf("point %d of the path", j))
  }
})
