# The table of a structural term's own parameters in summary() reads the
# variance of the selected convention. The reference is written in the test:
# the delta method over the block of vcov(readable = FALSE) that belongs to
# the term's parameters, and, for the inverse information, the explicit
# inverse of the joint information.

tbl_gas_data <- function(n = 150, omega = 0.3, kappa1 = 0.4, b1 = 0.7,
                         seed = 21) {
  set.seed(seed)
  f <- numeric(n)
  y <- numeric(n)
  s <- 0
  f0 <- omega / (1 - b1)
  for (t in seq_len(n)) {
    f[t] <- omega + kappa1 * (if (t > 1) s else 0) +
      b1 * (if (t > 1) f[t - 1] else f0)
    y[t] <- f[t] + stats::rnorm(1)
    s <- y[t] - f[t]
  }
  data.frame(t = seq_len(n), y = y, x = stats::rnorm(n))
}

# the reported quantities' standard errors by the delta method over the
# term's block of vcov(readable = FALSE), for the term `tm_name` of `mu`
tbl_delta_se <- function(fit, vt, tm_name) {
  V <- vcov(fit, vcov_type = vt, readable = FALSE)
  tm <- fit@spec@terms$mu[[tm_name]]
  design <- statmod_design(fit@spec)
  st <- statmod_structural_state(design)
  nm <- modelterms7::term_params(tm)
  free <- setdiff(nm, st$held[[tm_name]])
  rn <- unname(structural_tail_names(fit@spec, design))
  Vz <- V[rn, rn, drop = FALSE]
  rd <- modelterms7::term_readable(tm, st$zeta[[tm_name]])
  J <- rd$jacobian[, match(free, nm), drop = FALSE]
  stats::setNames(sqrt(diag(J %*% Vz %*% t(J))), rd$name)
}

test_that("the table of a gas() term reads the variance of every convention", {
  dd <- tbl_gas_data()
  fit <- statmod(y ~ x + gas(p = 1, q = 1, time = t) - 1,
                 distributions7::gaussian1_distrib(), dd)
  set.seed(1)
  fit <- add_bootstrap(fit, n = 8, resample = parametric_resample())
  types <- list(
    bayesian = bayesian_vcov(),
    frequentist = frequentist_vcov(),
    hac = frequentist_vcov(sandwich = hac_sandwich()),
    unconditional = unconditional_vcov(),
    bootstrap = bootstrap_vcov())
  tm_name <- attr(statmod_design(fit@spec), "structural")[[1L]]$term
  se <- list()
  for (k in names(types)) {
    tb <- statmod_structural_table(fit, 0.95, types[[k]])
    expect_true(all(is.finite(tb$se)), info = k)
    se[[k]] <- stats::setNames(tb$se, tb$name)
    # against the readable matrix, where it is a delta method too
    if (k != "bootstrap") {
      Vr <- vcov(fit, vcov_type = types[[k]], readable = TRUE)
      ref <- sqrt(diag(Vr))[paste0("mu:gas.", tb$name)]
      expect_equal(unname(tb$se), unname(ref), tolerance = 1e-10, info = k)
    }
    # and against the delta method written here
    expect_equal(se[[k]], tbl_delta_se(fit, types[[k]], tm_name)[tb$name],
                 tolerance = 1e-10, info = k)
  }
  # the conventions differ where they should: the HAC and the bootstrap
  # variances are not the inverse information
  expect_gt(max(abs(se$hac / se$bayesian - 1)), 0.05)
  expect_gt(max(abs(se$bootstrap / se$bayesian - 1)), 0.01)

  # the inverse information is unchanged: the tail of the explicit inverse
  spec <- fit@spec
  design <- statmod_design(spec)
  I <- statmod_full_information(spec, fit@coefficients, design)
  Vi <- solve(I)
  nb <- sum(vapply(design, function(d) d$npar, integer(1)))
  tail_i <- nb + 1:3
  tm <- spec@terms$mu[[tm_name]]
  rd <- modelterms7::term_readable(tm, statmod_structural_state(design)$zeta[[tm_name]])
  ref <- sqrt(diag(rd$jacobian %*% Vi[tail_i, tail_i] %*% t(rd$jacobian)))
  expect_equal(unname(se$bayesian), unname(ref), tolerance = 1e-10)
  expect_identical(
    statmod_structural_table(fit, 0.95, NULL)$se,
    statmod_structural_table(fit, 0.95, bayesian_vcov())$se)

  # the summary prints the table of the convention it was asked for
  sm <- summary(fit, vcov_type = types$hac)
  expect_equal(sm@structural$se, unname(se$hac), tolerance = 1e-10)
  sb <- summary(fit, vcov_type = types$bootstrap)
  expect_equal(sb@structural$se, unname(se$bootstrap), tolerance = 1e-10)
})

test_that("the break-point position of seg() reads the convention in the summary", {
  set.seed(3)
  x <- sort(stats::runif(150, 0, 10))
  y <- 1 + 0.5 * x + 1.5 * pmax(x - 6, 0) + stats::rnorm(150, sd = 0.4)
  fit <- statmod(y ~ seg(x), distributions7::gaussian1_distrib(),
                 data.frame(x = x, y = y))
  set.seed(1)
  fit <- add_bootstrap(fit, n = 8, resample = case_resample())
  for (vt in list(bayesian_vcov(), frequentist_vcov(), bootstrap_vcov())) {
    sm <- summary(fit, vcov_type = vt)
    blk <- Filter(function(b) identical(b$kind, "breakpoint"),
                  sm@tables$mu)[[1L]]$table
    ref <- sqrt(diag(vcov(fit, vcov_type = vt, readable = TRUE)))
    expect_equal(unname(blk$se), unname(ref[paste0("mu:seg.", rownames(blk))]),
                 tolerance = 1e-10, info = vt@name)
  }
})

test_that("a regime() table reads the bootstrap variance", {
  set.seed(3)
  n <- 150
  P <- matrix(c(0.93, 0.07, 0.06, 0.94), 2, 2, byrow = TRUE)
  s <- integer(n)
  s[1L] <- 1L
  for (t in 2:n) s[t] <- sample.int(2L, 1L, prob = P[s[t - 1L], ])
  dd <- data.frame(t = seq_len(n), y = stats::rnorm(n, c(0, 3)[s]))
  fit <- statmod(y ~ regime(2, time = t) - 1,
                 distributions7::gaussian1_distrib(), dd)
  set.seed(8)
  fit <- add_bootstrap(fit, n = 6)
  tb <- statmod_structural_table(fit, 0.95, bootstrap_vcov())
  tb0 <- statmod_structural_table(fit, 0.95, bayesian_vcov())
  expect_identical(tb$name, tb0$name)
  expect_true(any(is.finite(tb$se)))
  expect_false(isTRUE(all.equal(tb$se, tb0$se)))
  # the delta method over the stored replicates' covariance
  V <- vcov(fit, vcov_type = bootstrap_vcov(), readable = FALSE)
  design <- statmod_design(fit@spec)
  tm <- fit@spec@terms$mu[[names(fit@spec@terms$mu)[1L]]]
  nm <- modelterms7::term_params(tm)
  rd <- modelterms7::term_readable(
    tm, statmod_structural_state(design)$zeta[[1L]])
  rn <- rownames(V)[startsWith(rownames(V), "mu:") & !grepl("sigma", rownames(V))]
  rn <- tail(rn, length(nm))
  ref <- sqrt(pmax(diag(rd$jacobian %*% V[rn, rn] %*% t(rd$jacobian)), 0))
  fin <- is.finite(tb$se)
  expect_equal(unname(tb$se[fin]), unname(ref[fin]), tolerance = 1e-8)
})
