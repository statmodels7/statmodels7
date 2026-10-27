G <- distributions7::gaussian1_distrib()
set.seed(11)
n <- 120
pd <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n),
                 id = factor(sample(1:10, n, TRUE)), x = runif(n))
u0 <- rnorm(10)
u1 <- 0.5 * u0 + rnorm(10, sd = 0.5)
pd$y <- 1 + pd$x1 + 0.5 * pd$x2 + u0[pd$id] + u1[pd$id] * pd$x2 +
  sin(2 * pi * pd$x) + rnorm(n, sd = 0.5)

pctx <- function(fit, name, data = NULL, param = "mu") {
  list(param = param, name = name, vars = term_variables(
    fit@spec@terms[[param]][[name]], name), cols = NULL, data = data)
}
pterm <- function(fit, name, param = "mu") fit@spec@terms[[param]][[name]]

rname <- "ridge(~x1 + x2 + x3)"
fr <- statmod(y ~ ridge(~ x1 + x2 + x3), G, pd)
rname <- names(fr@spec@terms$mu)[2]
lname <- "lasso(~x1 + x2 + x3)"
fl <- statmod(y ~ lasso(~ x1 + x2 + x3), G, pd)
lname <- names(fl@spec@terms$mu)[2]

test_that("the unpenalized estimates of a ridge are the least squares ones", {
  d <- plot_data(plot_shrinkage(sort = FALSE), pterm(fr, rname), fr,
                 pctx(fr, rname))$estimates
  ols <- coef(lm(y ~ x1 + x2 + x3, pd))[-1]
  expect_equal(d$unpenalized, unname(ols), tolerance = 1e-8)
  b <- fr@coefficients$mu[-1]
  expect_equal(d$penalized, unname(b), tolerance = 1e-12)
  expect_true(all(abs(d$penalized) < abs(d$unpenalized)))
})

test_that("a lasso's zero coefficient has an unpenalized estimate", {
  d <- plot_data(plot_shrinkage(sort = FALSE), pterm(fl, lname), fl,
                 pctx(fl, lname))$estimates
  ols <- coef(lm(y ~ x1 + x2 + x3, pd))[-1]
  expect_equal(d$unpenalized, unname(ols), tolerance = 1e-8)
  expect_equal(d$penalized, unname(fl@coefficients$mu[-1]), tolerance = 1e-12)
})

test_that("more coefficients than the data identify is an error", {
  dup <- pd
  dup$x1b <- dup$x1
  f <- statmod(y ~ ridge(~ x1 + x1b + x2), G, dup)
  nm <- names(f@spec@terms$mu)[2]
  expect_error(autoplot(f, term = nm, plot_type = plot_shrinkage()),
               "cannot be inverted")
})

test_that("the quadratic path passes through the fit and the ridge solution", {
  lam <- hyper(fr)$estimate
  dp <- plot_data(plot_path(), pterm(fr, rname), fr, pctx(fr, rname))
  p <- dp$paths
  at <- p[abs(p$log_lambda - log(lam)) < 1e-9, ]
  expect_equal(at$estimate, unname(fr@coefficients$mu[-1]), tolerance = 1e-9)
  # closed form with the fitted dispersion, intercept unpenalized
  s2 <- fr@fitted$sigma[1]^2
  X <- cbind(1, as.matrix(pd[c("x1", "x2", "x3")]))
  P <- diag(c(0, 1, 1, 1))
  for (l in c(0.1 * lam, 10 * lam)) {
    ref <- solve(crossprod(X) / s2 + l * P, crossprod(X, pd$y) / s2)[-1]
    q <- plot_data(plot_path(log_lambda = log(l) + c(0, 1e-3)), pterm(fr, rname),
                   fr, pctx(fr, rname))$paths
    expect_equal(q$estimate[q$log_lambda == log(l)], drop(ref), tolerance = 1e-6)
  }
})

test_that("the quadratic path agrees with a refit at a fixed lambda", {
  lam <- hyper(fr)$estimate
  l <- 8 * lam
  f2 <- statmod(y ~ ridge(~ x1 + x2 + x3, lambda = l), G, pd)
  q <- plot_data(plot_path(log_lambda = log(l) + c(0, 1e-3)), pterm(fr, rname),
                 fr, pctx(fr, rname))$paths
  expect_equal(q$estimate[q$log_lambda == log(l)],
               unname(f2@coefficients$mu[-1]), tolerance = 0.02)
})

test_that("the lasso path passes through the fit and agrees with refits", {
  lam <- hyper(fl)$estimate
  dp <- plot_data(plot_path(), pterm(fl, lname), fl, pctx(fl, lname))
  p <- dp$paths
  at <- p[abs(p$log_lambda - log(lam)) < 1e-9, ]
  expect_equal(at$estimate, unname(fl@coefficients$mu[-1]), tolerance = 1e-6)
  expect_equal(dp$count$n[abs(dp$count$log_lambda - log(lam)) < 1e-9],
               sum(fl@coefficients$mu[-1] != 0))
  expect_equal(max(dp$count$n), 3)
  expect_equal(dp$count$n[which.max(dp$count$log_lambda)], 0)
  for (l in lam * c(0.3, 2)) {
    f2 <- statmod(y ~ lasso(~ x1 + x2 + x3, lambda = l), G, pd)
    q <- plot_data(plot_path(log_lambda = log(l) + c(0, 1e-3)),
                   pterm(fl, lname), fl, pctx(fl, lname))$paths
    expect_equal(q$estimate[q$log_lambda == log(l)],
                 unname(f2@coefficients$mu[-1]), tolerance = 0.05)
  }
})

rfit <- statmod(y ~ x1 + random(~ 1 | id), G, pd)
gname <- names(rfit@spec@terms$mu)[2]

test_that("the caterpillar holds the coefficients, their se and the raw means", {
  d <- plot_data(plot_caterpillar(sort = FALSE), pterm(rfit, gname), rfit,
                 pctx(rfit, gname))$effects
  b <- rfit@coefficients$mu[3:12]
  expect_equal(d$estimate, unname(b), tolerance = 1e-12)
  V <- vcov(rfit, vcov_type = bayesian_vcov(), readable = FALSE)
  nm <- paste0("mu:", modelterms7::term_coef_names(pterm(rfit, gname)))
  expect_equal(d$se, sqrt(diag(V[nm, nm])), ignore_attr = TRUE,
               tolerance = 1e-10)
  expect_equal(d$upper - d$lower, 2 * stats::qnorm(0.975) * d$se,
               tolerance = 1e-10)
  gm <- tapply(pd$y - rfit@fitted$mu, pd$id, mean) + b
  expect_equal(d$raw, as.vector(gm), tolerance = 1e-8)
  s <- plot_data(plot_caterpillar(), pterm(rfit, gname), rfit,
                 pctx(rfit, gname))$effects
  expect_equal(levels(s$level), as.character(1:10)[order(b)])
})

test_that("a random effect has no unpenalized estimate beside an intercept", {
  expect_error(autoplot(rfit, term = gname, plot_type = plot_shrinkage()),
               "cannot be inverted")
})

test_that("the prior density integrates to one and has the fitted scale", {
  d <- plot_data(plot_prior(), pterm(rfit, gname), rfit, pctx(rfit, gname))
  expect_equal(sum(d$density$density) * diff(d$density$x[1:2]), 1,
               tolerance = 1e-3)
  s <- hyper(rfit)$estimate[1]
  expect_equal(d$density$density, stats::dnorm(d$density$x, 0, s),
               tolerance = 1e-10)
})

sfit <- statmod(y ~ x1 + random(~ 1 + x2 | id), G, pd)
sname <- names(sfit@spec@terms$mu)[2]

test_that("a random slope gives a panel for each coefficient and an ellipse", {
  d <- plot_data(plot_caterpillar(), pterm(sfit, sname), sfit,
                 pctx(sfit, sname))$effects
  expect_equal(nlevels(d$coef), 2)
  expect_equal(nrow(d), 20)
  pc <- plot_data(plot_random_cov(), pterm(sfit, sname), sfit,
                  pctx(sfit, sname))
  S <- attr(pc, "Sigma")
  q <- rowSums((as.matrix(pc$ellipse) %*% solve(S)) * as.matrix(pc$ellipse))
  expect_equal(q, rep(stats::qchisq(0.95, 2), length(q)), tolerance = 1e-8)
  rd <- penalties7::penalty_readable(
    statmod_penalized(sfit@spec, statmod_design(sfit@spec))[[1]]$penalty,
    as.list(sfit@hyper$mu[[1]]))
  expect_equal(sqrt(diag(S)), unname(rd$value[1:2]), tolerance = 1e-10)
  expect_error(autoplot(rfit, term = gname, plot_type = plot_random_cov()),
               "at least two")
})

test_that("the curves of a random smooth are the group deviations", {
  cf <- statmod(y ~ s(x) + random(~ s(x) | id), G, pd)
  cn <- names(cf@spec@terms$mu)[3]
  d <- plot_data(plot_curves(n_grid = 5), pterm(cf, cn), cf,
                 pctx(cf, cn, pd))$curves
  expect_equal(nrow(d), 50)
  nd <- data.frame(x = d$x, id = d$group)
  full <- predict(cf, what = "mu", newdata = nd)
  zero <- predict(cf, what = "mu", newdata = nd, random = "zero")
  expect_equal(d$fit, unname(full - zero), tolerance = 1e-8)
  pop <- plot_data(plot_curves(n_grid = 5, population = TRUE), pterm(cf, cn),
                   cf, pctx(cf, cn, pd))
  expect_equal(pop$curves$fit, unname(full) - cf@coefficients$mu[[1]],
               tolerance = 1e-8)
  expect_equal(nrow(pop$population), 5)
  expect_error(autoplot(cf, term = cn), "does not keep its data")
})

test_that("the defaults and the smooth shrinkage run", {
  expect_s3_class(default_plot(pterm(fr, rname)), "statmodels7::PlotShrinkage")
  expect_s3_class(default_plot(pterm(rfit, gname)),
                  "statmodels7::PlotCaterpillar")
  sf <- statmod(y ~ x1 + s(x), G, pd)
  d <- plot_data(plot_shrinkage(), pterm(sf, "s(x)"), sf, pctx(sf, "s(x)"))
  expect_true(nrow(d$estimates) > 3)
  pdf(NULL)
  on.exit(dev.off())
  expect_s3_class(print(autoplot(fl, term = lname, plot_type = plot_path())),
                  "ggplot")
  expect_error(plot_path(log_lambda = c(2, 1)), "increasing")
  expect_error(plot_random_cov(which = c(1, 1)), "different")
})

test_that("a t prior has its own density", {
  tf <- statmod(y ~ x1 + random(~ 1 | id, distrib = distributions7::fixed(
    distributions7::student_t1_distrib(), mu = 0)), G, pd)
  nm <- names(tf@spec@terms$mu)[2]
  d <- plot_data(plot_prior(), pterm(tf, nm), tf, pctx(tf, nm))$density
  expect_true(all(is.finite(d$density)))
  expect_equal(sum(d$density) * diff(d$x[1:2]), 1, tolerance = 0.02)
})
