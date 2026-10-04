# The compiled route of the score-driven recursion's second, third and fourth
# orders (modelterms7's gas_curvature_gen.cpp) against the R route it
# mirrors, on both of the term's routes: the scalar one and the subformula
# one, with the autoregressive coefficients varying by observation at q = 3
# and q = 4, where the Levinson-Durbin map's third and fourth derivatives are
# not identically zero.

gen_twin_data <- function() {
  set.seed(11)
  m_grp <- 6
  ti <- 30
  dd <- data.frame(id = factor(rep(seq_len(m_grp), each = ti)),
                   t = rep(seq_len(ti), m_grp))
  n <- nrow(dd)
  dd$x <- stats::rnorm(n)
  dd$z <- stats::runif(n, -1, 1)
  dd$y <- stats::rpois(n, exp(0.8 + 0.2 * dd$x + 0.3 * sin(dd$t / 4)))
  dd
}

# the largest relative gap between the two routes over the three orders,
# at a point away from any optimum and along random directions
gen_twin_gap <- function(form, dd) {
  rel <- function(a, b) {
    a <- unlist(a)
    b <- unlist(b)
    max(abs(a - b)) / max(1, max(abs(b)))
  }
  spec <- statmod_spec(form, poisson_distrib(), dd)
  design <- statmod_design(spec)
  hyper <- statmod_hyper_start(spec, design)
  obj <- statmod_objective(spec, hyper, design, FALSE, "bartlett")
  beta <- statmod_start(spec, design, obj, NULL)
  sst <- statmod_structural_state(design)
  key <- names(sst$zeta)[[1L]]
  free <- setdiff(names(sst$zeta[[key]]), sst$held[[key]])
  sst$zeta[[key]][free] <- sst$zeta[[key]][free] +
    stats::runif(length(free), -0.1, 0.1)
  sst$key <- NULL
  coef <- obj$split(beta + stats::runif(length(beta), -0.05, 0.05))
  jd <- joint_design_rows(spec, design, coef)
  st <- structural_grad_parts_impl(spec, design, coef, jd,
                                   diag(length(jd$keep)))
  f <- jd$f
  vfull <- stats::rnorm(jd$mf)
  wfull <- stats::rnorm(jd$mf)
  call_t <- function(gen, ...) {
    gen(f$tm, f$eta_static, spec@response,
        function(e, i) st$sd_at[i], function(e, i) st$cd_at[i],
        f$psi, st$w * st$s_at, st$seed, ...)
  }
  data_args <- list(score_values = st$sd_at, curvature_values = st$cd_at,
                    blocks_data = st$blocks_data)
  r2 <- call_t(modelterms7::term_curvature, st$blocks(NULL))
  c2 <- do.call(call_t, c(list(modelterms7::term_curvature, st$blocks(NULL)),
                          data_args))
  r3 <- call_t(modelterms7::term_third, st$blocks(vfull), vfull)
  c3 <- do.call(call_t, c(list(modelterms7::term_third, st$blocks(vfull),
                               vfull), data_args))
  D5 <- distributions7::distrib_deriv5(spec@distrib, spec@response,
                                       jd$ev$theta, scale = "link")
  dr4 <- filter_driving(spec, jd$ev$theta, jd$ap, filter_scaling(f$tm),
                        st$gl, st$H, st$D3, st$D4, D5)
  blk4 <- .structural_blocks(jd$params, jd$ap, jd$V, dr4$H, dr4$D3, dr4$D4,
                             jd$n, dr4$D5)
  bd4 <- structural_blocks_data(jd$params, jd$ap, jd$V, dr4$H, dr4$D3, jd$n,
                                D4 = dr4$D4, D5 = dr4$D5)
  dirs <- list(vfull, wfull)
  r4 <- call_t(modelterms7::term_fourth, blk4(dirs), dirs)
  c4 <- call_t(modelterms7::term_fourth, blk4(dirs), dirs,
               score_values = st$sd_at, curvature_values = st$cd_at,
               blocks_data = bd4)
  c(order2 = rel(r2, c2), order3 = rel(r3$curvature, c3$curvature),
    dphi3 = rel(r3$dphi, c3$dphi), order4 = rel(r4$curvature, c4$curvature),
    dphi4 = rel(r4$dphi, c4$dphi), dpsi4 = rel(r4$dpsi, c4$dpsi))
}

test_that("the compiled higher orders reproduce the R recursion", {
  skip_on_cran()
  dd <- gen_twin_data()
  forms <- list(
    y ~ x + gas(p = 1, q = 1, by = id, time = t),
    y ~ x + gas(p = 2, q = 2, by = id, time = t, scaling = 1),
    y ~ x + gas(p = 1, q = 1, by = id, time = t, omega ~ random(~1 | id)),
    y ~ x + gas(p = 1, q = 1, by = id, time = t, pacf1 ~ random(~1 | id),
                scaling = 0.5),
    y ~ x + gas(p = 1, q = 3, by = id, time = t, pacf2 ~ z),
    y ~ x + gas(p = 1, q = 4, by = id, time = t, pacf1 ~ z, pacf3 ~ z))
  for (fm in forms) {
    gaps <- gen_twin_gap(fm, dd)
    expect_lt(max(gaps), 1e-12, label = deparse(fm[[3L]])[1L])
  }
})
