# The coefficients a marginal criterion estimates: reml(marginal =) and
# ml(marginal =).

marginal_data <- function(n = 60, seed = 1) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n), z = stats::runif(n))
  d$y <- 1 + 2 * d$x - d$z + stats::rnorm(n, sd = 1.5)
  d
}

test_that("reml() returns the sigma and the criterion of lm(REML)", {
  d <- marginal_data()
  g <- distributions7::gaussian1_distrib()
  l <- stats::lm(y ~ x + z, d)
  f <- statmod(y ~ x + z, g, d)
  expect_equal(exp(f@coefficients$sigma), summary(l)$sigma, tolerance = 1e-7)
  expect_equal(f@criterion, as.numeric(stats::logLik(l, REML = TRUE)),
               tolerance = 1e-8)
  expect_equal(unname(f@coefficients$mu), unname(stats::coef(l)),
               tolerance = 1e-8)
})

test_that("ml() returns sqrt(rss/n), and marginal = 'none' the joint mode", {
  d <- marginal_data()
  g <- distributions7::gaussian1_distrib()
  l <- stats::lm(y ~ x + z, d)
  ml_sigma <- sqrt(sum(stats::residuals(l)^2) / nrow(d))
  f_ml <- statmod(y ~ x + z, g, d, outer_criterion = ml())
  f_none <- statmod(y ~ x + z, g, d, outer_criterion = reml(marginal = "none"))
  expect_equal(exp(f_ml@coefficients$sigma), ml_sigma, tolerance = 1e-7)
  expect_equal(exp(f_none@coefficients$sigma), ml_sigma, tolerance = 1e-7)
  # "none" with no penalty runs no criterion at all
  expect_true(is.na(f_none@criterion))
})

test_that("with every parameter named, reml() and ml() are one criterion", {
  set.seed(11)
  m <- 20L
  gr <- factor(rep(seq_len(m), each = 8L))
  x <- stats::runif(length(gr))
  eta <- 0.3 + 0.8 * x + stats::rnorm(m, 0, 0.5)[gr]
  dat <- data.frame(g = gr, x = x,
                    y = stats::rgamma(length(gr), shape = 3,
                                      rate = 3 / exp(eta)))
  d <- distributions7::gamma1_distrib()
  fr <- statmod(y ~ x + random(~ 1 | g), d, dat,
                outer_criterion = reml(marginal = "all"))
  fm <- statmod(y ~ x + random(~ 1 | g), d, dat,
                outer_criterion = ml(marginal = "all"))
  expect_equal(fr@criterion, fm@criterion, tolerance = 1e-10)
  expect_equal(unlist(fr@coefficients), unlist(fm@coefficients),
               tolerance = 1e-8)
})

test_that("the variance carries the criterion's curvature in sigma", {
  d <- marginal_data()
  g <- distributions7::gaussian1_distrib()
  l <- stats::lm(y ~ x + z, d)
  f <- statmod(y ~ x + z, g, d)
  V <- vcov(f, readable = FALSE)
  expect_equal(unname(V[1:3, 1:3]), unname(stats::vcov(l)), tolerance = 1e-8)
  # -d2 l_M / d(log sigma)^2 is 2(n - p) at the REML optimum
  expect_equal(V[4, 4], 1 / (2 * (nrow(d) - 3)), tolerance = 1e-8)
  expect_identical(statmod_certificate(f)$state, "converged")
})

test_that("the exact gradient and Hessian over (eta, gamma) match differences", {
  set.seed(2)
  n <- 200
  d <- data.frame(x = stats::runif(n), z = stats::runif(n),
                  g = factor(sample(15, n, TRUE)))
  d$y <- sin(3 * d$x) + stats::rnorm(15, 0, 0.4)[d$g] +
    stats::rnorm(n, sd = exp(-0.5 + 0.8 * d$z))
  G <- distributions7::gaussian1_distrib()
  method <- reml()
  fit <- statmod(y ~ x + random(~ 1 | g) | sigma ~ z, G, d,
                 outer_criterion = method)
  spec <- fit@spec
  design <- statmod_design(spec)
  blocks <- statmod_blocks(spec, design)
  idx <- outer_hyper_index(spec, blocks)
  gam <- marginal_coords(spec, design, method)
  nh <- nrow(idx)
  ng <- length(gam$where)
  expect_identical(ng, 2L)
  basis <- integrated_basis(spec, design, method@kind, gamma = TRUE)
  beta0 <- unlist(fit@coefficients[spec@distrib@params], use.names = FALSE)
  set.seed(5)
  v0 <- c(hyper_to_eta(fit@hyper, idx), beta0[gam$where]) +
    stats::rnorm(nh + ng, 0, 0.15)
  at <- function(v) {
    hy <- eta_to_hyper(v[seq_len(nh)], idx, fit@hyper)
    sp <- statmod_hold(spec, gam, v[nh + seq_len(ng)])
    b <- beta0
    b[gam$where] <- v[nh + seq_len(ng)]
    res <- statmod_alternate(sp, design, blocks, hy,
                             iwls(tol = 1e-12, maxit = 500), b, FALSE, "opg",
                             500, 1e-12, verbosity(0))
    # the mode polished by Newton steps on the free coordinates: a reference
    # refitted to a stall is biased, not noisy
    par <- res$par
    free <- setdiff(seq_along(par), gam$where)
    for (it in 1:4) {
      gr <- res$obj$gr(par)
      he <- as.matrix(res$obj$he(par))
      par[free] <- par[free] - solve(he[free, free], gr[free])
    }
    cf <- res$obj$split(par)
    list(sp = sp, hy = hy, cf = cf,
         val = statmod_marginal(sp, design, cf, hy, method, "opg", basis)$value,
         g = statmod_marginal_grad(sp, design, cf, hy, method, idx, basis,
                                   gam = gam))
  }
  a <- at(v0)
  H <- statmod_marginal_hess(a$sp, design, a$cf, a$hy, method, idx, basis,
                             gam = gam)
  gap <- function(h) {
    fd_g <- vapply(seq_along(v0), function(j) {
      e <- replace(numeric(length(v0)), j, h)
      (at(v0 + e)$val - at(v0 - e)$val) / (2 * h)
    }, numeric(1))
    fd_H <- sapply(seq_along(v0), function(j) {
      e <- replace(numeric(length(v0)), j, h)
      (at(v0 + e)$g - at(v0 - e)$g) / (2 * h)
    })
    fd_H <- (fd_H + t(fd_H)) / 2
    c(g = max(abs(a$g - fd_g)) / max(abs(a$g)),
      H = max(abs(H - fd_H)) / max(abs(H)))
  }
  g1 <- gap(1e-2)
  g2 <- gap(1e-3)
  expect_lt(g2[["g"]], 1e-5)
  expect_lt(g2[["H"]], 1e-5)
  # O(h^2): a missing term would be flat in the step
  expect_gt(g1[["g"]] / g2[["g"]], 30)
  expect_gt(g1[["H"]] / g2[["H"]], 30)
})

test_that("marginal is checked against the family", {
  d <- marginal_data()
  g <- distributions7::gaussian1_distrib()
  expect_error(statmod(y ~ x, g, d, outer_criterion = reml(marginal = "nu")),
               "Its parameters are: mu, sigma")
  expect_error(reml(marginal = c("none", "sigma")), "alone")
  expect_error(reml(marginal = NA_character_), "'marginal' must be")
  d$w <- stats::runif(nrow(d))
  expect_error(statmod(y ~ x | sigma ~ te(z, w), g, d,
                       outer_criterion = reml(marginal = "sigma")),
               "null space")
})

test_that("a parameter the default names keeps the joint mode where not covered", {
  d <- marginal_data()
  d$w <- stats::runif(nrow(d))
  g <- distributions7::gaussian1_distrib()
  spec <- statmod_spec(y ~ x | sigma ~ te(z, w), g, d)
  design <- statmod_design(spec)
  mc <- marginal_coords(spec, design, reml())
  expect_length(mc$where, 0L)
  expect_match(mc$skipped[["sigma"]], "null space")
})

test_that("a block that is a working linearization keeps the joint mode", {
  # A sharp jseg()'s block is not a Jacobian, so a determinant over its
  # columns reads no curvature: measured, two fits reaching the same
  # break-point read -132.9 and -109.8 according to where they began.
  set.seed(2)
  d <- data.frame(x = sort(stats::runif(200, 0, 10)))
  d$y <- 1 + 0.4 * d$x + 1.5 * (d$x > 6) + stats::rnorm(200, sd = 0.3)
  g <- distributions7::gaussian1_distrib()
  spec <- statmod_spec(y ~ jseg(x), g, d)
  design <- statmod_design(spec)
  mc <- marginal_coords(spec, design, reml())
  expect_length(mc$where, 0L)
  expect_match(mc$skipped[["sigma"]], "working linearization")
  expect_error(statmod(y ~ jseg(x, n_boot = 0), g, d,
                       outer_criterion = reml(marginal = "sigma")),
               "working linearization")
})

test_that("an estimated coefficient at its limit is named a boundary", {
  # nu estimated on the criterion runs to its gaussian limit on gaussian
  # data, where the criterion's curvature in it vanishes: the certificate
  # names it by equation and coefficient, and gives it no hyperparameter key
  set.seed(11)
  d <- data.frame(x = stats::runif(300))
  d$y <- 1 + 2 * d$x + stats::rnorm(300, sd = 0.5)
  f <- statmod(y ~ x, distributions7::student_t1_distrib(), d)
  ct <- statmod_certificate(f)
  expect_identical(ct$state, "boundary")
  expect_identical(ct$boundary, "nu/(Intercept)")
  expect_length(ct$boundary_key, 0L)
})
