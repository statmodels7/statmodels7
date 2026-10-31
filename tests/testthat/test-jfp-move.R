# The third and fourth derivatives of the Jeffreys penalty where a predictor
# is not linear in the coefficients and where the block's rows move, against
# an independent twin that enumerates every assignment of the directions to
# the three factors of the block's matrix, every set partition in the
# predictors (Faa di Bruno) and every set partition of the log-determinant's
# derivative, and against a numerical Jacobian of the order below.

jfp_twin_partitions <- function(x) {
  if (!length(x)) return(list(list()))
  first <- x[1L]
  out <- list()
  for (p in jfp_twin_partitions(x[-1L])) {
    out[[length(out) + 1L]] <- c(list(first), p)
    for (i in seq_along(p)) {
      q <- p
      q[[i]] <- c(first, q[[i]])
      out[[length(out) + 1L]] <- q
    }
  }
  out
}

jfp_twin_perms <- function(x) {
  if (length(x) <= 1L) return(list(x))
  out <- list()
  for (i in seq_along(x)) {
    for (p in jfp_twin_perms(x[-i])) out[[length(out) + 1L]] <- c(x[i], p)
  }
  out
}

# the derivative of log|I_tilde| along the slots `slots` (names into `dirs`),
# for the one block of the state
jfp_twin_dlogdet <- function(st, dirs, slots) {
  bs <- st$blocks[[1L]]
  params <- st$params
  np <- length(params)
  a <- jfp_live(bs)
  n <- st$spec@n_obs
  w <- st$spec@weights
  units <- refresh_units(st$spec, st$design, st$coef, params, st$npar,
                         st$offs)
  memo <- new.env()
  # the derivative of every predictor along a list of directions
  eta_d <- function(ds) {
    key <- paste("eta", paste(names(ds), collapse = ","))
    if (!is.null(memo[[key]])) return(memo[[key]])
    r <- length(ds)
    out <- lapply(seq_len(np), function(cc) {
      if (r == 1L) {
        return(as.numeric(as.matrix(st$design[[params[cc]]]$X) %*%
                            ds[[1L]][st$offs[cc] + seq_len(st$npar[cc])]))
      }
      acc <- numeric(n)
      for (un in units) {
        if (un$a != cc) next
        dd <- lapply(ds, function(z) z[un$ra])
        D <- switch(r - 1L,
          modelterms7::term_block_deriv(un$built, coef = un$bt, v = dd[[1]]),
          modelterms7::term_block_deriv2(un$built, coef = un$bt, v = dd[[1]],
                                         u = dd[[2]]),
          modelterms7::term_block_deriv3(un$built, coef = un$bt, v = dd[[1]],
                                         u = dd[[2]], s = dd[[3]]),
          modelterms7::term_block_deriv4(un$built, coef = un$bt, v = dd[[1]],
                                         u = dd[[2]], s = dd[[3]], r = dd[[4]]))
        acc <- acc + as.numeric(as.matrix(D) %*% dd[[r]])
      }
      acc
    })
    memo[[key]] <- out
    out
  }
  Jarr <- function(cs) {
    k <- length(cs)
    if (k == 0L) return(-rep_len(st$E[[hess_key(params, a, a)]], n))
    D <- jfp_darray(st, k)
    key <- switch(k,
      distributions7::dexpected_key(params, a, a, cs[1]),
      distributions7::d2expected_key(params, a, a, cs[1], cs[2]),
      distributions7::d3expected_key(params, a, a, cs[1], cs[2], cs[3]),
      distributions7::d4expected_key(params, a, a, cs[1], cs[2], cs[3], cs[4]))
    -rep_len(D[[key]], n)
  }
  # J along a set of slots: Faa di Bruno over the set partitions
  J_d <- function(S) {
    if (!length(S)) return(Jarr(integer(0)))
    out <- numeric(n)
    for (pp in jfp_twin_partitions(S)) {
      etas <- lapply(pp, function(b) eta_d(dirs[b]))
      combos <- as.matrix(expand.grid(rep(list(seq_len(np)), length(pp))))
      for (j in seq_len(nrow(combos))) {
        cs <- combos[j, ]
        f <- Jarr(cs)
        for (b in seq_along(pp)) f <- f * etas[[b]][[cs[b]]]
        out <- out + f
      }
    }
    out
  }
  V_d <- function(S) {
    if (!length(S)) return(bs$rows[[a]])
    if (!identical(bs$block$kind, "refresh")) return(NULL)
    dd <- lapply(dirs[S], function(z) z[bs$tidx])
    D <- switch(length(S),
      modelterms7::term_block_deriv(bs$tm, coef = bs$bt, v = dd[[1]]),
      modelterms7::term_block_deriv2(bs$tm, coef = bs$bt, v = dd[[1]],
                                     u = dd[[2]]),
      modelterms7::term_block_deriv3(bs$tm, coef = bs$bt, v = dd[[1]],
                                     u = dd[[2]], s = dd[[3]]),
      modelterms7::term_block_deriv4(bs$tm, coef = bs$bt, v = dd[[1]],
                                     u = dd[[2]], s = dd[[3]], r = dd[[4]]))
    as.matrix(D)[, bs$pos, drop = FALSE] %*% bs$block$C
  }
  # the block's matrix along a set of slots: every assignment to (J, L, R)
  I_d <- function(S) {
    p <- ncol(bs$rows[[a]])
    out <- matrix(0, p, p)
    for (code in seq_len(3^length(S)) - 1L) {
      g <- integer(length(S))
      x <- code
      for (i in seq_along(S)) {
        g[i] <- x %% 3L
        x <- x %/% 3L
      }
      L <- V_d(S[g == 1L])
      R <- V_d(S[g == 2L])
      if (is.null(L) || is.null(R)) next
      out <- out + crossprod(L * (w * J_d(S[g == 0L])), R)
    }
    out
  }
  Ii <- solve(I_d(character(0)))
  tot <- 0
  for (pp in jfp_twin_partitions(slots)) {
    A <- lapply(pp, function(b) Ii %*% I_d(b))
    k <- length(pp)
    ords <- if (k == 1L) list(1L) else
      lapply(jfp_twin_perms(seq_len(k)[-1L]), function(o) c(1L, o))
    for (o in ords) {
      P <- diag(nrow(Ii))
      for (i in o) P <- P %*% A[[i]]
      tot <- tot + (-1)^(k - 1L) * sum(diag(P))
    }
  }
  tot
}

jfp_twin_check <- function(fm, distrib, d) {
  fit <- suppressWarnings(statmod(fm, distrib = distrib, data = d,
                                  outer_criterion = NULL))
  sp <- fit@spec
  des <- statmod_design(sp)
  params <- sp@distrib@params
  cf <- fit@coefficients[params]
  npar <- lengths(cf)
  m <- sum(npar)
  offs <- cumsum(npar) - npar
  split <- function(u) stats::setNames(lapply(seq_along(params), function(a)
    u[offs[a] + seq_len(npar[a])]), params)
  set.seed(1)
  u0 <- unlist(cf, use.names = FALSE) + stats::rnorm(m, sd = 0.03)
  st_at <- function(u) jfp_state(sp, des, split(u), jfp_blocks(sp, des),
                                 emethod = jfp_emethod(sp@distrib))
  v <- stats::rnorm(m)
  w <- stats::rnorm(m)
  st <- st_at(u0)
  D3 <- jfp_deriv3_dir(st, v)
  D4 <- jfp_deriv4_dir(st, v, w)
  T3 <- T4 <- matrix(0, m, m)
  for (k in seq_len(m)) {
    for (l in k:m) {
      dirs <- list(k = replace(numeric(m), k, 1), l = replace(numeric(m), l, 1),
                   v = v, w = w)
      T3[k, l] <- T3[l, k] <- jfp_twin_dlogdet(st, dirs, c("k", "l", "v"))
      T4[k, l] <- T4[l, k] <- jfp_twin_dlogdet(st, dirs,
                                               c("k", "l", "v", "w"))
    }
  }
  a <- st$blocks[[1L]]$block$power
  n3 <- matrix(numDeriv::jacobian(function(t)
    as.vector(jfp_hessian(st_at(u0 + t * v))$hessian), 0), m)
  t3 <- jfp_terms3(st, v)
  t3$dense <- list()
  st0 <- st_at(u0)
  st0$units <- list()
  st0$blocks[[1L]]$block$kind <- "static"
  list(twin3 = max(abs(D3 - a * T3)) / max(abs(T3)),
       twin4 = max(abs(D4 - a * T4)) / max(abs(T4)),
       num3 = max(abs(D3 - n3)) / max(abs(n3)),
       ctrl_dense = max(abs(jfp_terms_matrix(st, t3) - n3)) / max(abs(n3)),
       ctrl_static = max(abs(jfp_deriv4_dir(st0, v, w) - a * T4)) /
         max(abs(T4)))
}

test_that("the moving third and fourth derivatives match an enumerating twin", {
  skip_if_not_installed("numDeriv")
  set.seed(5)
  n <- 50
  d <- data.frame(x = stats::runif(n, 0.2, 3), g = stats::rnorm(n),
                  h = stats::rbinom(n, 1, 0.4), z = stats::rbinom(n, 1, 0.3))
  d$y <- stats::rpois(n, exp(log(3 + 0.5 * d$g) - 0.6 * d$x))
  d$ys <- 1 + 0.5 * d$x + 1.5 * pmax(d$x - 1.5 - 0.3 * d$h, 0) +
    stats::rnorm(n, sd = 0.3)
  d$yn <- stats::rnbinom(n, size = 3, mu = exp(1 - 0.5 * d$x))
  P <- distributions7::poisson_distrib()
  cases <- list(
    # rows that move, in the subformula of nl()
    list(y ~ 0 + nl(~ exp(a - r * x), a ~ 1 + jfp(~ g)), P),
    # rows that move, in the break-point of a smoothed seg()
    list(ys ~ x + seg(x, psi ~ 1 + jfp(~ h),
                      smoothed = numericals7::smooth_probit()),
         distributions7::gaussian1_distrib()),
    # a block that does not move, its information moved by nl() in mu
    list(yn ~ 0 + nl(~ exp(a - r * x)) | theta ~ jfp(~ z),
         distributions7::negbin2_distrib()))
  for (cs in cases) {
    r <- jfp_twin_check(cs[[1]], cs[[2]], d)
    expect_lt(r$twin3, 1e-10)
    expect_lt(r$twin4, 1e-10)
    expect_lt(r$num3, 1e-7)
    expect_gt(r$ctrl_dense, 1e-2)
    expect_gt(r$ctrl_static, 1e-2)
  }
})

test_that("the outer gradient and Hessian are exact with a block that moves", {
  skip_if_not_installed("numDeriv")
  set.seed(6)
  n <- 200
  d <- data.frame(x = stats::runif(n, -2, 2), h = stats::rbinom(n, 1, 0.3),
                  g = factor(rep(1:20, 10)))
  d$y <- stats::rbinom(n, 1, stats::plogis(2 * (d$x - 0.3) +
                                             stats::rnorm(20, sd = 0.5)[d$g]))
  d$y[d$h == 1] <- as.numeric(d$x[d$h == 1] > 0.5)
  B <- distributions7::bernoulli_distrib()
  fm <- y ~ 0 + nl(~ b * (x - e), b ~ 1 + jfp(~ h),
                   start = list(b = 1, e = 0)) + random(~1 | g)
  method <- ml(information = observed_info())
  fit0 <- statmod(fm, distrib = B, data = d, outer_criterion = NULL)
  spec <- fit0@spec
  design <- statmod_design(spec)
  idx <- outer_hyper_index(spec, statmod_blocks(spec, design))
  expect_true(outer_gradient_ok(spec, design, idx, method, 2L))
  basis <- integrated_basis(spec, design, method@kind)
  # the mode on the fit's own specification, polished by Newton on the
  # exact gradient
  at <- function(eta) {
    hy <- eta_to_hyper(eta, idx, fit0@hyper)
    obj <- statmod_objective(spec, hy, design, FALSE)
    par <- obj$stack(fit_at_hyper_spec(spec, design, hy, B))
    for (k in 1:8) {
      gg <- obj$gr(par)
      if (max(abs(gg)) < 1e-10) break
      J <- numDeriv::jacobian(obj$gr, par)
      par <- par - as.numeric(solve((J + t(J)) / 2, gg))
    }
    list(hy = hy, cf = obj$split(par))
  }
  fn <- function(eta) {
    a <- at(eta)
    statmod_marginal(spec, design, a$cf, a$hy, method, basis = basis)$value
  }
  gr <- function(eta) {
    a <- at(eta)
    statmod_marginal_grad(spec, design, a$cf, a$hy, method, idx, basis)
  }
  eta <- hyper_to_eta(fit0@hyper, idx) + 0.3
  a <- at(eta)
  g <- gr(eta)
  gn <- (fn(eta + 1e-4) - fn(eta - 1e-4)) / 2e-4
  expect_lt(abs(g - gn) / abs(gn), 1e-7)
  H <- statmod_marginal_hess(spec, design, a$cf, a$hy, method, idx, basis)
  Hn <- (gr(eta + 1e-4) - gr(eta - 1e-4)) / 2e-4
  expect_lt(abs(H - Hn) / abs(Hn), 1e-7)
})
