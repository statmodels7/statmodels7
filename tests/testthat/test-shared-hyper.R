# Data both files below use: two smooths whose wanted smoothing is far apart,
# so that sharing is a real restriction and not a formality.
shared_data <- function(n = 300, seed = 21) {
  set.seed(seed)
  d <- data.frame(x = runif(n, 0, 10), z = runif(n, 0, 10))
  d$y <- sin(3 * d$x) + 0.2 * d$z + rnorm(n, sd = 0.4)
  d
}

test_that("a label collapses the outer index to one row over two members", {
  d <- shared_data()
  sp <- statmod_spec(y ~ s(x, k = 8, id = "L") + s(z, k = 8, id = "L"),
                     gaussian1_distrib(), d)
  bl <- statmod_blocks(sp, statmod_design(sp))
  idx <- outer_hyper_index(sp, bl)
  expect_identical(nrow(idx), 1L)
  mem <- index_members(idx)
  expect_identical(nrow(mem), 2L)
  expect_identical(mem$row, c(1L, 1L))
  expect_identical(mem$name, c("lambda", "lambda"))

  # and without the label there are two rows, each its own member: the member
  # table is then the index, which is what lets every loop over members be
  # the loop over rows that was there before sharing existed
  sp2 <- statmod_spec(y ~ s(x, k = 8) + s(z, k = 8), gaussian1_distrib(), d)
  bl2 <- statmod_blocks(sp2, statmod_design(sp2))
  idx2 <- outer_hyper_index(sp2, bl2)
  expect_identical(nrow(idx2), 2L)
  mem2 <- index_members(idx2)
  expect_identical(mem2$row, 1:2)
  expect_identical(mem2[, c("parameter", "term", "name")],
                   idx2[, c("parameter", "term", "name")])
})

test_that("the estimated value is written under every member's own key", {
  d <- shared_data()
  fit <- statmod(y ~ s(x, k = 8, id = "L") + s(z, k = 8, id = "L"),
                 gaussian1_distrib(), d, outer_criterion = reml())
  h <- hyper(fit)
  expect_identical(nrow(h), 2L)
  # THE REPRESENTATION CHECK: sharing does not merge the penalties, so each
  # keeps its own entry in the store and the two carry the same number
  expect_identical(h$estimate[[1L]], h$estimate[[2L]])
  expect_false(any(h$held))

  # the negative control: unshared, the two are genuinely different
  free <- statmod(y ~ s(x, k = 8) + s(z, k = 8), gaussian1_distrib(), d,
                  outer_criterion = reml())
  hf <- hyper(free)
  expect_gt(abs(diff(log(hf$estimate))), 1)
})

test_that("sharing is a restriction, so the criterion cannot rise", {
  d <- shared_data()
  a <- statmod(y ~ s(x, k = 8) + s(z, k = 8), gaussian1_distrib(), d,
               outer_criterion = reml())
  b <- statmod(y ~ s(x, k = 8, id = "L") + s(z, k = 8, id = "L"),
               gaussian1_distrib(), d, outer_criterion = reml())
  expect_gt(a@criterion, b@criterion)
})

test_that("the effective degrees of freedom stay per term", {
  # one hyperparameter does not make one count: the blocks differ, so the
  # traces differ, and a count taken per hyperparameter would report the two
  # as equal
  d <- shared_data()
  fit <- statmod(y ~ s(x, k = 8, id = "L") + s(z, k = 8, id = "L"),
                 gaussian1_distrib(), d, outer_criterion = reml())
  e <- fit@edf[grepl("^s\\(", fit@edf$term), ]
  expect_identical(nrow(e), 2L)
  expect_false(isTRUE(all.equal(e$edf[[1L]], e$edf[[2L]])))
})

test_that("the group's gradient is the sum of its members'", {
  # the two references share no arithmetic with each other: the members are
  # read from a DIFFERENT index, built from a formula carrying no label
  d <- shared_data()
  meth <- reml()
  lam <- 0.35

  mode_at <- function(v) {
    f <- bquote(y ~ s(x, k = 8, lambda = .(v)) + s(z, k = 8, lambda = .(v)))
    statmod(stats::as.formula(f), gaussian1_distrib(), d)@coefficients
  }
  prep <- function(f) {
    sp <- statmod_spec(f, gaussian1_distrib(), d)
    de <- statmod_design(sp)
    bl <- statmod_blocks(sp, de)
    list(spec = sp, design = de, blocks = bl,
         idx = outer_hyper_index(sp, bl))
  }
  set_all <- function(pp, v) {
    h <- statmod_hyper_start(pp$spec, pp$design)
    for (p in names(h)) {
      for (k in names(h[[p]])) {
        if ("lambda" %in% names(h[[p]][[k]])) h[[p]][[k]][["lambda"]] <- v
      }
    }
    h
  }
  grad_at <- function(pp, v) {
    h <- set_all(pp, v)
    statmod_marginal_grad(pp$spec, pp$design, mode_at(v), h, meth, pp$idx,
                          integrated_basis(pp$spec, pp$design, meth@kind))
  }

  sh <- prep(y ~ s(x, k = 8, id = "L") + s(z, k = 8, id = "L"))
  fr <- prep(y ~ s(x, k = 8) + s(z, k = 8))
  g_sh <- grad_at(sh, lam)
  g_fr <- grad_at(fr, lam)
  expect_length(g_sh, 1L)
  expect_length(g_fr, 2L)
  expect_equal(g_sh, sum(g_fr), tolerance = 1e-10)
  # and the members are large and of OPPOSITE sign, so the agreement is not
  # two small numbers happening to match
  expect_gt(max(abs(g_fr)), 1)
  expect_lt(prod(g_fr), 0)
})

test_that("a shared axis is swept once by a path", {
  set.seed(4)
  n <- 300
  X <- matrix(rnorm(n * 6), n, 6)
  colnames(X) <- c("w1", "w2", "w3", "v1", "v2", "v3")
  d <- as.data.frame(X)
  d$y <- 1.5 * X[, 1] - 1.2 * X[, 4] + rnorm(n, sd = 0.6)

  sp <- statmod_spec(y ~ lasso(~ w1 + w2 + w3, id = "K") +
                       lasso(~ v1 + v2 + v3, id = "K"),
                     gaussian1_distrib(), d)
  de <- statmod_design(sp)
  bl <- statmod_blocks(sp, de)
  rows <- path_rows(sp, bl, statmod_hyper_start(sp, de), bic())
  expect_identical(nrow(rows), 1L)
  expect_identical(nrow(rows$members[[1L]]), 2L)

  fit <- statmod(y ~ lasso(~ w1 + w2 + w3, id = "K") +
                   lasso(~ v1 + v2 + v3, id = "K"),
                 gaussian1_distrib(), d, outer_criterion = bic())
  h <- hyper(fit)
  expect_identical(nrow(h), 2L)
  expect_identical(h$estimate[[1L]], h$estimate[[2L]])
})

test_that("a label may not span the two hyperparameter machines", {
  set.seed(4)
  n <- 200
  d <- data.frame(w1 = rnorm(n), w2 = rnorm(n), v1 = rnorm(n), v2 = rnorm(n))
  d$y <- 1.5 * d$w1 - 1.2 * d$v1 + rnorm(n, sd = 0.6)
  expect_error(
    statmod(y ~ lasso(~ w1 + w2, id = "K") + ridge(~ v1 + v2, id = "K"),
            gaussian1_distrib(), d, outer_criterion = bic()),
    "smooth penalty and on a kinked one")
})

test_that("a label may not tie hyperparameters whose links differ", {
  # one free value has to land in both domains, and (0, Inf) is not (0, 1)
  set.seed(4)
  n <- 200
  d <- data.frame(w1 = rnorm(n), w2 = rnorm(n))
  d$y <- 1.5 * d$w1 + rnorm(n, sd = 0.6)
  expect_error(
    index_group(
      data.frame(parameter = c("mu", "mu"), term = c("a", "b"),
                 name = c("lambda", "alpha"), stringsAsFactors = FALSE),
      list(linkfunctions7::log_link(), linkfunctions7::logit_link()),
      c("A", "A")),
    "One value cannot lie in both domains")
})

test_that("a fit with no label is untouched by the machinery", {
  # the negative control for the whole change: with nothing shared, every
  # quantity is what it was, and the exact Hessian is still available
  d <- shared_data()
  fit <- statmod(y ~ s(x, k = 8) + s(z, k = 8), gaussian1_distrib(), d,
                 outer_criterion = reml())
  sp <- fit@spec
  de <- statmod_design(sp)
  bl <- statmod_blocks(sp, de)
  idx <- outer_hyper_index(sp, bl)
  expect_true(outer_gradient_ok(sp, de, idx, reml(), 1L))
  expect_true(outer_gradient_ok(sp, de, idx, reml(), 2L))
})

# The Hessian at a shared row, and the identity that judges it. Holding the
# members of a group at one value makes the shared model the free one read
# along the diagonal eta -> (eta, eta), so its criterion is a composition and
# its Hessian is the free one conjugated by the sharing map J. That identity
# shares no arithmetic with the assembly, and reading ONE member of a group --
# which is what a table keyed by the pair of hyperparameter names could hold --
# fails it by a quarter of the entry rather than by a rounding.
shared_hess <- function(formula, data, eta) {
  spec <- statmod_spec(formula, gaussian1_distrib(), data)
  des <- statmod_design(spec)
  idx <- outer_hyper_index(spec, statmod_blocks(spec, des))
  hy <- eta_to_hyper(eta, idx, statmod_hyper_start(spec, des))
  r <- fit_at_hyper(formula, gaussian1_distrib(), data, hy)
  as.matrix(statmod_marginal_hess(r$spec, r$design, r$coefficients, hy,
                                  reml(), idx, NULL))
}

shared_hess_data <- function(n = 300, seed = 21) {
  set.seed(seed)
  d <- data.frame(x = runif(n, 0, 10), z = runif(n, 0, 10),
                  w = runif(n, 0, 10))
  d$y <- sin(3 * d$x) + 0.2 * d$z + 0.5 * cos(d$w) + rnorm(n, sd = 0.4)
  d
}

test_that("a shared group has an exact outer Hessian at both orders", {
  d <- shared_data()
  sp <- statmod_spec(y ~ s(x, k = 8, id = "L") + s(z, k = 8, id = "L"),
                     gaussian1_distrib(), d)
  de <- statmod_design(sp)
  idx <- outer_hyper_index(sp, statmod_blocks(sp, de))
  expect_true(outer_gradient_ok(sp, de, idx, reml(), 1L))
  expect_true(outer_gradient_ok(sp, de, idx, reml(), 2L))
})

test_that("the shared Hessian is the unshared one read along the diagonal", {
  d <- shared_hess_data()
  eta <- c(0.5, 0.5)
  a <- shared_hess(y ~ s(x, k = 8, id = "L") + s(z, k = 8, id = "L") +
                     s(w, k = 8), d, eta)
  b <- shared_hess(y ~ s(x, k = 8) + s(z, k = 8) + s(w, k = 8), d,
                   c(eta[[1L]], eta[[1L]], eta[[2L]]))
  J <- matrix(c(1, 1, 0, 0, 0, 1), nrow = 3L)
  expect_equal(a, t(J) %*% b %*% J, tolerance = 1e-10)

  # and the identity is not vacuous: the entry a per-term table could have
  # held is one member own curvature, which is a different number
  expect_gt(abs(a[1L, 1L] - b[1L, 1L]), 0.2)
})

test_that("two hyperparameters of ONE unit may sit in two rows", {
  # the shape a table keyed by the pair of NAMES could not represent at all:
  # a tensor smooth whose first margin is shared with another term and whose
  # second is free, so the within-unit cross derivative spans two rows
  d <- shared_hess_data(n = 400)
  eta <- c(0.4, -0.2)
  fs <- y ~ te(x, z, k = c(5, 5), id = c(lambda1 = "L")) +
    s(w, k = 8, id = "L")
  sp <- statmod_spec(fs, gaussian1_distrib(), d)
  mem <- index_members(outer_hyper_index(sp, statmod_blocks(sp,
                                          statmod_design(sp))))
  expect_identical(mem$row, c(1L, 2L, 1L))

  a <- shared_hess(fs, d, eta)
  b <- shared_hess(y ~ te(x, z, k = c(5, 5)) + s(w, k = 8), d,
                   c(eta[[1L]], eta[[2L]], eta[[1L]]))
  J <- matrix(c(1, 0, 1, 0, 1, 0), nrow = 3L)
  expect_equal(a, t(J) %*% b %*% J, tolerance = 1e-10)
  # the cross entry is the unit own, and it is not zero
  expect_gt(abs(a[1L, 2L]), 0.1)
})

test_that("a shared hyperparameter carries a standard error", {
  d <- shared_data()
  fit <- statmod(y ~ s(x, k = 8, id = "L") + s(z, k = 8, id = "L"),
                 gaussian1_distrib(), d, outer_criterion = reml())
  V <- statmod_hyper_vcov(fit@spec, statmod_design(fit@spec),
                          fit@coefficients, fit@hyper, fit@methods$outer)
  expect_false(is.null(V))
  expect_identical(dim(V), c(1L, 1L))
  expect_gt(V[[1L, 1L]], 0)
  # and the summary reports it rather than an empty cell. hyper() carries no
  # se column at all, so an assertion made there passes on NULL whatever the
  # fit did; the number a reader sees is in the summary table.
  bl <- Filter(function(z) identical(z$kind, "smooth"),
               summary(fit)@tables$mu)[[1L]]
  r <- bl$table[bl$table$role == "estimated", ]
  expect_identical(nrow(r), 1L)
  expect_true(is.finite(r$se))
  expect_lt(r$lower, r$estimate)
  expect_gt(r$upper, r$estimate)
})

test_that("a value held on one member is held for the whole group", {
  d <- shared_data()
  f <- y ~ s(x, k = 8, lambda = 2, id = "L") + s(z, k = 8, id = "L")
  sp <- statmod_spec(f, gaussian1_distrib(), d)
  de <- statmod_design(sp)

  # the store starts with the held value under BOTH keys. $ is not used on
  # it: an entry is an atomic NAMED VECTOR, so $ is an error there and a
  # partial match elsewhere
  h0 <- statmod_hyper_start(sp, de)
  expect_identical(unname(h0$mu[[1L]][["lambda"]]), 2)
  expect_identical(unname(h0$mu[[2L]][["lambda"]]), 2)

  # and there is nothing left for either machine to estimate
  bl <- statmod_blocks(sp, de)
  expect_identical(nrow(outer_hyper_index(sp, bl)), 0L)

  fit <- statmod(f, gaussian1_distrib(), d, outer_criterion = reml())
  h <- hyper(fit)
  expect_identical(h$estimate, c(2, 2))
  expect_true(all(h$held))
  expect_identical(h$id, c("L", "L"))
})

test_that("two members held at different values are a contradiction", {
  d <- shared_data()
  expect_error(
    statmod(y ~ s(x, k = 8, lambda = 2, id = "L") +
              s(z, k = 8, lambda = 5, id = "L"),
            gaussian1_distrib(), d, outer_criterion = reml()),
    "is held at 2 in")
  # the same value twice says the same thing twice and is fine
  expect_no_error(
    statmod(y ~ s(x, k = 8, lambda = 2, id = "L") +
              s(z, k = 8, lambda = 2, id = "L"),
            gaussian1_distrib(), d, outer_criterion = reml()))
})

test_that("a held label reaches a path too", {
  set.seed(4)
  n <- 250
  X <- matrix(rnorm(n * 4), n, 4)
  colnames(X) <- c("w1", "w2", "v1", "v2")
  d <- as.data.frame(X)
  d$y <- 1.5 * X[, 1] - 1.2 * X[, 3] + rnorm(n, sd = 0.6)
  fit <- statmod(y ~ lasso(~ w1 + w2, lambda = 20, id = "K") +
                   lasso(~ v1 + v2, id = "K"),
                 gaussian1_distrib(), d, outer_criterion = bic())
  h <- hyper(fit)
  expect_identical(h$estimate, c(20, 20))
  expect_true(all(h$held))
})
