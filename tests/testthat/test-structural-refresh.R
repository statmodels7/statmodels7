# A term that recomputes its block (nl(), seg()) beside a penalized filter:
# the exact outer gradient of the structural branch reads the block at the
# mode, the mode's movement by the likelihood's own curvature, and the
# movement of K with the block (structural_refresh_move())

# each fit takes a minute or more, so a case is fitted once per file
structural_refresh_cache <- new.env()

structural_refresh_case <- function(fm, distrib) {
  key <- paste(deparse(fm), collapse = "")
  if (!is.null(structural_refresh_cache[[key]])) {
    return(structural_refresh_cache[[key]])
  }
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
  d$y <- 1.5 * exp(-0.8 * d$xs) + lev[d$id] + e + stats::rnorm(n, sd = 0.5)
  d$yn <- stats::rnbinom(n, mu = exp(0.5 + 1.5 * exp(-0.8 * d$xs) +
                                       0.4 * lev[d$id] + 0.5 * e),
                         size = exp(1 + 0.4 * d$z))
  f <- suppressWarnings(statmod(fm, distrib = distrib, data = d))
  sp <- f@spec
  des <- statmod_design(sp)
  sst <- statmod_structural_state(des)
  key <- names(sst$zeta)[1]
  free <- setdiff(names(sst$zeta[[key]]), sst$held[[key]])
  nb <- sum(vapply(des, function(z) z$npar, integer(1)))
  obj <- statmod_objective(sp, f@hyper, des)
  jp <- statmod_joint_pieces(sp, des, obj, f@hyper)
  u <- c(unlist(f@coefficients[sp@distrib@params], use.names = FALSE),
         as.numeric(sst$zeta[[key]][free]))
  cf_at <- function(uu) {
    jp$setz(uu[nb + seq_along(free)])
    obj$split(uu[seq_len(nb)])
  }
  out <- list(f = f, sp = sp, des = des, nb = nb, jp = jp, u = u,
              cf_at = cf_at)
  assign(key, out, envir = structural_refresh_cache)
  out
}

test_that("the block's movement beside a penalized filter is exact", {
  skip_if_not_installed("numDeriv")
  cases <- list(
    list(y ~ nl(~ b * exp(-r * xs), start = list(b = 1, r = 0.5)) +
           gas(p = 1, q = 1, omega ~ random(~ 1 | id), time = t, by = id),
         distributions7::gaussian1_distrib()),
    list(yn ~ nl(~ b * exp(-r * xs), start = list(b = 1, r = 0.5)) +
           gas(p = 1, q = 1, omega ~ 1 + random(~ 1 | id), time = t,
               by = id) | theta ~ z,
         distributions7::negbin2_distrib()))
  for (cs in cases) {
    k <- structural_refresh_case(cs[[1]], cs[[2]])
    sp <- k$sp
    des <- k$des
    nb <- k$nb
    # away from the mode, where the scores that weight the block's second
    # derivative are not small
    set.seed(2)
    u <- k$u + 0.1 * stats::rnorm(length(k$u))
    cf <- k$cf_at(u)
    K <- statmod_marginal_full(sp, des, cf, k$f@hyper, NULL)
    desr <- statmod_design_at(sp, cf, des)
    jd <- joint_design_rows(sp, desr, cf)
    rm <- structural_refresh_move(sp, desr, cf, jd)
    expect_false(is.null(rm))
    # the mode moves by the penalized objective's own curvature
    Hn <- numDeriv::jacobian(k$jp$gr, u)
    Hn <- (Hn + t(Hn)) / 2
    Km <- K
    Km[seq_len(nb), seq_len(nb)] <- Km[seq_len(nb), seq_len(nb)] + rm$Dm
    expect_lt(max(abs(Km - Hn)) / max(abs(Hn)), 1e-7)
    # the control: K alone is the Gauss-Newton matrix of the block
    expect_gt(max(abs(K - Hn)) / max(abs(Hn)), 1e-4)
    # K moves along v: the trace against M, against Richardson differences
    # of K over the joint vector
    set.seed(1)
    v <- stats::rnorm(length(u))
    M <- crossprod(matrix(stats::rnorm(length(u)^2), length(u))) / length(u)
    Kf <- function(uu) statmod_marginal_full(sp, des, k$cf_at(uu),
                                             k$f@hyper, NULL)
    h <- 1e-4
    dKn <- (4 * (Kf(u + h / 2 * v) - Kf(u - h / 2 * v)) / h -
              (Kf(u + h * v) - Kf(u - h * v)) / (2 * h)) / 3
    cf <- k$cf_at(u)
    st <- structural_grad_parts_impl(sp, desr, cf, jd, M)
    base <- sum(st$u * v) + structural_chain_extra(sp, desr, jd, M, st, v)
    mv <- rm$trace(v, M)
    expect_lt(abs(base + mv - sum(M * dKn)) / abs(sum(M * dKn)), 1e-7)
    # the control: without the block's movement
    expect_gt(abs(base - sum(M * dKn)) / abs(sum(M * dKn)), 1e-5)
  }
})

test_that("the block's second movement beside a penalized filter is exact", {
  # one case: on the negative binomial the filter's own directions dominate
  # every trace against a random M, and a difference cannot separate the
  # block's share from them
  cases <- list(
    list(y ~ nl(~ b * exp(-r * xs), start = list(b = 1, r = 0.5)) +
           gas(p = 1, q = 1, omega ~ random(~ 1 | id), time = t, by = id),
         distributions7::gaussian1_distrib()))
  for (cs in cases) {
    k <- structural_refresh_case(cs[[1]], cs[[2]])
    sp <- k$sp
    des <- k$des
    nb <- k$nb
    set.seed(2)
    u <- k$u + 0.1 * stats::rnorm(length(k$u))
    cf <- k$cf_at(u)
    desr <- statmod_design_at(sp, cf, des)
    jd <- joint_design_rows(sp, desr, cf)
    rm <- structural_refresh_move(sp, desr, cf, jd)
    set.seed(1)
    v1 <- stats::rnorm(length(u))
    v2 <- stats::rnorm(length(u))
    M <- crossprod(matrix(stats::rnorm(length(u)^2), length(u))) / length(u)
    st <- structural_grad_parts_impl(sp, desr, cf, jd, M)
    rs <- structural_refresh_second(sp, desr, cf, jd, M, rm, list(v1, v2))
    expect_false(is.null(rs))
    # the movement of K with the rows, as a matrix
    expect_lt(abs(sum(M * rs$R[[1]]) - rm$trace(v1, M)) /
                abs(rm$trace(v1, M)), 1e-9)
    # the third derivative of the penalized objective along v2 and v1:
    # K + Dm moves with the rows and with the block's second derivative
    Kmode <- function(uu) {
      cfu <- k$cf_at(uu)
      dsu <- statmod_design_at(sp, cfu, des)
      K <- statmod_marginal_full(sp, dsu, cfu, k$f@hyper, NULL)
      r <- structural_refresh_move(sp, dsu, cfu, joint_design_rows(sp, dsu, cfu))
      K[seq_len(nb), seq_len(nb)] <- K[seq_len(nb), seq_len(nb)] + r$Dm
      K
    }
    h <- 1e-3
    rich <- function(f) (4 * (f(h / 2) - f(-h / 2)) / h -
                           (f(h) - f(-h)) / (2 * h)) / 3
    dKm <- rich(function(s) as.numeric(Kmode(u + s * v2) %*% v1))
    k_cfu <- k$cf_at(u)
    an <- as.numeric((structural_dk_matrix(sp, desr, jd, st, v2) +
                        rs$R[[2]]) %*% v1) + rs$mode(1, 2)
    expect_lt(max(abs(an - dKm)) / max(abs(dKm)), 1e-8)
    # the control: without the movement of the block's own term
    expect_gt(max(abs(an - rs$mode(1, 2) - dKm)) / max(abs(dKm)), 1e-4)
    # the second movement of K along v1 and v2, against M
    Kf <- function(uu) statmod_marginal_full(sp, des, k$cf_at(uu),
                                             k$f@hyper, NULL)
    f2 <- function(s) sum(M * (Kf(u + s * (v1 + v2)) - Kf(u + s * (v1 - v2)) -
                                 Kf(u - s * (v1 - v2)) +
                                 Kf(u - s * (v1 + v2)))) / (4 * s^2)
    num <- (4 * f2(h / 2) - f2(h)) / 3
    k_cfu <- k$cf_at(u)
    blk4 <- structural_blocks4(sp, jd, st)
    fixed <- structural_chain_extra2(sp, desr, jd, M, st, blk4, v1, v2)
    an2 <- fixed + rs$trace(1, 2)
    expect_lt(abs(an2 - num) / abs(num), 1e-8)
    # the control: the rows held fixed
    expect_gt(abs(fixed - num) / abs(num), 1e-4)
  }
})

test_that("the outer Hessian is exact beside a recomputed block", {
  k <- structural_refresh_case(
    y ~ nl(~ b * exp(-r * xs), start = list(b = 1, r = 0.5)) +
      gas(p = 1, q = 1, omega ~ random(~ 1 | id), time = t, by = id),
    distributions7::gaussian1_distrib())
  des <- k$des
  idx <- outer_hyper_index(k$sp, statmod_blocks(k$sp, des))
  m <- reml(information = observed_info())
  expect_true(outer_gradient_ok(k$sp, des, idx, m, 1L))
  expect_true(outer_gradient_ok(k$sp, des, idx, m, 2L))
})

test_that("a held break-point beside a penalized filter leaves K", {
  # jump() holds its break-point's slot, whose row of K is zero; the
  # criterion's value leaves it out of the determinant, and so do the
  # gradient, the Hessian and the joint smoother
  kj <- structural_refresh_case(
    y ~ jump(xs) + gas(p = 1, q = 1, omega ~ random(~ 1 | id), time = t,
                       by = id),
    distributions7::gaussian1_distrib())
  sp <- kj$sp
  des <- kj$des
  cf <- kj$cf_at(kj$u)
  hy <- kj$f@hyper
  m <- reml(information = observed_info())
  idj <- outer_hyper_index(sp, statmod_blocks(sp, des))
  expect_true(length(term_held_stack(sp, des)) > 0)
  expect_true(outer_gradient_ok(sp, des, idj, m, 2L))
  g <- statmod_structural_grad(sp, des, cf, hy, m, idj)
  expect_true(!is.null(g) && all(is.finite(g)))
  h <- statmod_structural_hess(sp, des, cf, hy, m, idj)
  expect_true(!is.null(h) && all(is.finite(h)))
  js <- joint_smoother_diag(sp, cf, des, hy)
  expect_false(isTRUE(js$failed))
  expect_equal(unname(js$beta[term_held_stack(sp, des)]), 1)
})

test_that("beside an unpenalized filter the mode moves by the total score", {
  skip_if_not_installed("numDeriv")
  # mode_curvature() weighted by the total derivative in the static
  # predictor (statmod_score_obs()) gives the joint objective's curvature;
  # weighted by the family's score alone it does not
  k <- structural_refresh_case(
    yn ~ nl(~ b * exp(-r * xs), start = list(b = 1, r = 0.5)) +
      gas(p = 1, q = 1, time = t, by = id) | theta ~ random(~ 1 | id),
    distributions7::negbin2_distrib())
  sp <- k$sp
  des <- k$des
  nb <- k$nb
  set.seed(2)
  u <- k$u + 0.1 * stats::rnorm(length(k$u))
  cf <- k$cf_at(u)
  desr <- statmod_design_at(sp, cf, des)
  params <- sp@distrib@params
  npar <- vapply(desr, function(z) z$npar, integer(1))
  offs <- cumsum(npar) - npar
  w <- sp@weights
  gl <- lapply(statmod_score_obs(sp, cf, desr), function(x) x / w)
  J <- statmod_full_information(sp, cf, des)
  S <- as_dense(statmod_penalty_at(sp, cf, k$f@hyper, des, "hessian"))
  ib <- seq_len(nb)
  Jt <- J
  Jt[ib, ib] <- Jt[ib, ib] + S +
    mode_curvature(sp, desr, cf, params, npar, offs, nb, gl = gl)
  Js <- J
  Js[ib, ib] <- Js[ib, ib] + S +
    mode_curvature(sp, desr, cf, params, npar, offs, nb)
  Hn <- numDeriv::jacobian(k$jp$gr, u)
  Hn <- (Hn + t(Hn)) / 2
  # on the block's own entries, where the two weightings differ
  blk <- unlist(lapply(attr(desr, "refresh"), function(r)
    offs[match(r$param, params)] + desr[[r$param]]$blocks[[r$term]]))
  e_t <- max(abs((Jt - Hn)[blk, blk])) / max(abs(Hn[blk, blk]))
  e_s <- max(abs((Js - Hn)[blk, blk])) / max(abs(Hn[blk, blk]))
  expect_lt(max(abs(Jt - Hn)) / max(abs(Hn)), 1e-7)
  expect_lt(e_t, 1e-7)
  expect_gt(e_s, 100 * e_t)
})
