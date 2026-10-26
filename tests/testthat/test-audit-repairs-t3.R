# Repairs found by the documentation audit of tranche 3 (inference and
# prediction). Every block fails on 08b94b63.

g1 <- distributions7::gaussian1_distrib()

t3_lin <- function(n = 40, seed = 2) {
  set.seed(seed)
  d <- data.frame(x = runif(n))
  d$y <- 1 + d$x + rnorm(n)
  d
}

t3_gas <- function() {
  set.seed(4)
  n <- 100
  e <- rnorm(n)
  f <- numeric(n)
  y <- numeric(n)
  for (t in seq_len(n)) {
    if (t > 1) f[t] <- 0.4 * e[t - 1] + 0.7 * f[t - 1]
    y[t] <- 2 + f[t] + e[t]
  }
  d <- data.frame(y = y, t = seq_len(n), x = rnorm(n))
  statmod(y ~ x + gas(p = 1, q = 1, time = t), g1, d)
}

t3_ranef <- function(distrib = NULL) {
  set.seed(6)
  d <- data.frame(g = factor(rep(1:10, each = 6)), x = rnorm(60))
  d$y <- 1 + d$x + rnorm(10)[d$g] + rnorm(60, sd = 0.5)
  f <- if (is.null(distrib)) y ~ x + random(~1 | g) else
    y ~ x + random(~1 | g, distrib = distrib)
  environment(f) <- environment()
  statmod(f, g1, d)
}

t3_smooth_data <- function() {
  set.seed(1)
  d <- data.frame(x = runif(200), z = runif(200))
  d$y <- sin(2 * pi * d$x) + rnorm(200, sd = 0.3)
  d
}

test_that("rstatmod() keeps the supplied columns in every replicate", {
  set.seed(1)
  d0 <- data.frame(z = runif(30))
  s <- rstatmod(y ~ x + z, g1, d0, n_sim = 2,
                covariates = list(x = function(n) rnorm(n)),
                par = list(mu = c(1, 2, 3), sigma = 0))
  expect_identical(s$data[[2]]$z, d0$z)
  expect_false(identical(s$data[[1]]$x, s$data[[2]]$x))
})

test_that("a count with a fractional part is rejected, not truncated", {
  d <- t3_lin(30)
  expect_error(rstatmod(y ~ x, g1, d, n_sim = 2.5,
                        par = list(mu = c(1, 2), sigma = 0)),
               "'n_sim' must be a positive whole number")
  expect_error(rstatmod(y ~ x, g1, d, n = 30.7,
                        par = list(mu = c(1, 2), sigma = 0)),
               "'n' is 30.7")
  expect_error(rstatmod(y ~ 1, g1, n = 10.5, par = list(mu = 1, sigma = 0)),
               "'n' must be a positive whole number")
  fit <- statmod(y ~ x, g1, d)
  expect_error(simulate(fit, nsim = 2.5), "'nsim' must be a positive whole")
})

test_that("a par function that returns a missing value says so", {
  d <- t3_lin(30)
  expect_error(rstatmod(y ~ x, g1, d,
                        par = list(mu = function(k) c(1, NA), sigma = 0)),
               "returned a missing value")
})

test_that("fitted() rejects an argument it would ignore", {
  fit <- statmod(y ~ x, g1, t3_lin())
  expect_error(fitted(fit, newdata = data.frame(x = c(0, 1))),
               "fitted\\(\\) takes 'what' alone")
  expect_length(fitted(fit, "mu"), 40L)
})

test_that("predict() rejects arguments that no route reads", {
  fit <- statmod(y ~ x, g1, t3_lin())
  nd <- data.frame(x = c(0, 1))
  expect_error(predict(fit, what = "mu", newdata = nd, bogus = 3),
               "takes no argument 'bogus'")
  expect_error(predict(fit, what = "mu", newdata = nd, vcov_type = unconditional_vcov()),
               "read only with")
  expect_error(predict(fit, what = "mu", se = TRUE, readable = FALSE),
               "takes no argument 'readable'")
  a <- predict(fit, what = "mu", newdata = nd, se = TRUE, vcov_type = frequentist_vcov())
  expect_true(all(is.finite(a$se)))
})

test_that("predict() validates the level", {
  fit <- statmod(y ~ x, g1, t3_lin())
  expect_error(predict(fit, what = "mu", se = TRUE, level = 95),
               "strictly between 0 and 1")
  expect_error(predict(fit, what = "mu", se = TRUE, level = -1),
               "strictly between 0 and 1")
})

test_that("the forecast warning names only an equation that a filter reaches", {
  fg <- t3_gas()
  fut <- data.frame(t = 101:103, x = 0)
  expect_no_warning(predict(fg, what = "sigma", newdata = fut, se = TRUE))
  expect_warning(predict(fg, what = "mu", newdata = fut, se = TRUE),
                 "past the series")
})

test_that("a continued filter keeps its static predictor apart", {
  fg <- t3_gas()
  fut <- data.frame(t = 101:103, x = 0)
  ns <- asNamespace("statmodels7")
  ss <- ns$spec_at(fg, fut, need_response = FALSE)
  ep <- ns$statmod_eta_continued(fg, ss, statmod_design(ss), deriv = FALSE)
  expect_false(isTRUE(all.equal(ep$eta$mu, ep$eta_static$mu)))
  b <- fg@coefficients$mu
  expect_equal(ep$eta_static$mu, rep(unname(b[[1L]]), 3))
})

test_that("an empty newdata is named in the message", {
  fit <- statmod(y ~ x, g1, t3_lin())
  expect_error(predict(fit, what = "mu", newdata = t3_lin()[0, ]),
               "'newdata' has no rows")
})

test_that("a parameter with no coefficients has a standard error of zero", {
  fit <- statmod(y ~ x | sigma ~ 0, g1, t3_lin())
  p <- predict(fit, what = "sigma", se = TRUE)
  expect_equal(p$se, rep(0, 40))
  expect_equal(p$lower, p$fit)
})

test_that("an unknown 'what' is named the same way under every random mode", {
  fr <- t3_ranef()
  nd <- data.frame(x = 0, g = factor("new"))
  expect_error(predict(fr, what = "foo", newdata = nd, random = "marginal"),
               "is neither a parameter of this distribution nor a")
})

test_that("a group interval under a Student t prior answers a missing row", {
  ft <- t3_ranef(distributions7::fixed(distributions7::student_t1_distrib(),
                                       mu = 0))
  nd <- data.frame(x = c(0.5, NA), g = factor("new"))
  p <- predict(ft, what = "mu", newdata = nd, random = "zero",
               interval = "group", se = TRUE)
  expect_true(is.finite(p$lower[1]) && is.finite(p$upper[1]))
  expect_true(is.na(p$lower[2]) && is.na(p$upper[2]))
})

test_that("a Laplace prior is not averaged as a Gaussian", {
  lp <- distributions7::fixed(distributions7::laplace_distrib(), mu = 0)
  fl <- suppressWarnings(t3_ranef(lp))
  nd <- data.frame(x = 0, g = factor("new"))
  expect_error(predict(fl, what = "mu", newdata = nd, random = "marginal"),
               "is not Gaussian")
})

test_that("the note on failed restricted refits has one negation", {
  src <- paste(deparse(statmodels7:::restricted_stat_rows), collapse = "\n")
  expect_false(grepl("s did not\"", src, fixed = TRUE))
})

test_that("a smooth under aic() is described as minimized, not as a path", {
  f <- statmod(y ~ s(x, bspline_smooth(k = 10)), g1, t3_smooth_data(),
               outer_criterion = aic())
  nt <- summary(f)@notes
  expect_false(any(grepl("path over", nt)))
  expect_true(any(grepl("minimizing that\n  criterion", nt, fixed = TRUE)))
})

test_that("the second member of a shared hyperparameter has its error", {
  f <- statmod(y ~ s(x, bspline_smooth(k = 8), id = "L") +
                 s(z, bspline_smooth(k = 8), id = "L"), g1, t3_smooth_data())
  s <- summary(f)
  bl <- Filter(function(b) any(b$table$role == "estimated"), s@tables$mu)
  se <- vapply(bl, function(b) b$table$se[b$table$role == "estimated"][[1L]],
               numeric(1))
  expect_length(se, 2L)
  expect_true(all(is.finite(se)))
  expect_equal(se[[1]], se[[2]])
  expect_false(any(grepl("did not reach a maximum", s@notes)))
})

test_that("a hyperparameter at an edge gets the boundary note alone", {
  set.seed(9)
  dn <- data.frame(x = runif(150))
  dn$y <- rnorm(150)
  f <- statmod(y ~ s(x, bspline_smooth(k = 8)), g1, dn)
  expect_identical(statmod_certificate(f)$state, "boundary")
  nt <- summary(f)@notes
  expect_true(any(grepl("run to an edge", nt)))
  expect_false(any(grepl("did not reach a maximum", nt)))
})

test_that("a restricted test reaches the rows of a break-point term", {
  set.seed(3)
  d <- data.frame(x = runif(150))
  d$y <- 1 + 2 * pmax(d$x - 0.5, 0) + rnorm(150, sd = 0.2)
  f <- statmod(y ~ seg(x), g1, d)
  pick <- function(s) Filter(function(b) identical(b$term, "seg(x)"),
                            s@tables$mu)[[1L]]$table
  tw <- pick(summary(f))
  tl <- pick(summary(f, test = lr_test()))
  r <- match("gamma1", tw$name)
  expect_false(isTRUE(all.equal(tw$statistic[r], tl$statistic[r])))
  expect_true(is.na(tl$statistic[match("psi1", tl$name)]))
})

test_that("a small edf correction is printed with its digits", {
  src <- paste(deparse(statmodels7:::summary.StatmodFit), collapse = "\n")
  expect_false(grepl("carry %.3f", src, fixed = TRUE))
})

test_that("the certificate names the coefficients as the reason", {
  src <- paste(deparse(statmodels7:::certificate_core), collapse = "\n")
  expect_true(grepl("differenced over\\n  hyperparameters alone", src,
                    fixed = TRUE))
})

test_that("format_duration() chooses the unit after rounding", {
  expect_identical(statmodels7:::format_duration(c(999.6e-6, 0.9999, 59.99,
                                                   3599.9, 1.5)),
                   c("1 ms", "1 s", "1 min", "1 h", "1.5 s"))
})

test_that("loglik() rejects an unnamed par", {
  fit <- statmod(y ~ x, g1, t3_lin())
  expect_error(loglik(fit, par = c(0, 0)), "must be a named list")
  expect_error(loglik(fit, par = list(c(0, 0), 0)), "must be a named list")
  expect_equal(loglik(fit, par = list(mu = c(0, 0), sigma = 0)),
               sum(dnorm(t3_lin()$y, 0, 1, log = TRUE)))
})

test_that("the marginal log-likelihood counts what it estimates", {
  d <- t3_smooth_data()
  fm <- statmod(y ~ s(x, bspline_smooth(k = 10)), g1, d, outer_criterion = ml())
  # intercept, the smooth's straight line, sigma's intercept and lambda
  expect_identical(attr(logLik(fm, loglik_type = "marginal"), "df"), 4)
  fh <- statmod(y ~ s(x, bspline_smooth(k = 10), hyper = c(lambda = 3)) +
                  s(z, bspline_smooth(k = 6)), g1, d, outer_criterion = ml())
  # mu's intercept, the two straight lines, sigma's intercept and the one
  # lambda that is estimated
  expect_identical(attr(logLik(fh, loglik_type = "marginal"), "df"), 5)
  # a random intercept counts as lme4 counts it: two fixed effects, the
  # residual scale and one variance (lme4::lmer(y ~ x + (1 | g), REML =
  # FALSE) gives df 4 and the same log-likelihood, -190.6884, on these data)
  set.seed(2)
  dr <- data.frame(g = factor(rep(1:20, each = 6)), x = rnorm(120))
  dr$y <- 1 + dr$x + rnorm(20)[dr$g] + rnorm(120)
  fr <- statmod(y ~ x + random(~1 | g), g1, dr, outer_criterion = ml())
  expect_identical(attr(logLik(fr, loglik_type = "marginal"), "df"), 4)
  expect_equal(as.numeric(logLik(fr, loglik_type = "marginal")), -190.6884,
               tolerance = 1e-6)
})

test_that("the inner reading leaves out the direction a kink holds", {
  set.seed(3)
  d <- data.frame(x = runif(150))
  d$y <- 1 + 2 * pmax(d$x - 0.5, 0) + rnorm(150, sd = 0.2)
  f <- statmod(y ~ seg(x), g1, d)
  ns <- asNamespace("statmodels7")
  de <- statmod_design(f@spec)
  obj <- ns$statmod_objective(f@spec, f@hyper, de, FALSE, distributions7::expect_bartlett())
  b <- obj$stack(f@coefficients)
  kc <- ns$kink_constraints(f@spec, de, f@coefficients,
                            obj$split(seq_along(b)), length(b))
  expect_identical(nrow(kc), 1L)
  ct <- statmod_certificate(f)
  expect_lt(ct$inner[["gradient"]], 1e-3)
  expect_lt(ct$mode_error, 1e-6)
})

test_that("a group interval past the link's domain is NA with a warning", {
  set.seed(1)
  g <- factor(rep(1:30, each = 5))
  x <- rnorm(150)
  b <- rnorm(30)
  d <- data.frame(y = rpois(150, (1.2 + 0.2 * x + b[g])^2), x = x, g = g)
  pd <- distributions7::poisson_distrib(link_mu = linkfunctions7::sqrt_link())
  f <- statmod(y ~ x + random(~1 | g), pd, d)
  nd <- data.frame(x = c(-3, 0), g = factor("new"))
  pl <- suppressWarnings(predict(f, what = "link:mu", newdata = nd,
                                 random = "zero", interval = "group"))
  expect_true(pl$lower[[1]] < 0)
  expect_warning(p <- predict(f, what = "mu", newdata = nd, random = "zero",
                              interval = "group"),
                 "leaves the domain")
  expect_true(is.na(p$lower[[1]]) && is.na(p$upper[[1]]))
})

test_that("a fit that converged is not reported as not converging", {
  f <- statmod(y ~ s(x, bspline_smooth(k = 10)), g1, t3_smooth_data(),
               outer_criterion = aic())
  expect_true(f@converged)
  nt <- summary(f)@notes
  expect_false(any(grepl("did not converge", nt)))
  expect_false(any(grepl("still moving", nt)))
  expect_true(any(grepl("not certified as a maximum of an outer", nt,
                        fixed = TRUE)))
})

test_that("a hyperparameter that no criterion moved is marked start", {
  f <- statmod(y ~ s(x, bspline_smooth(k = 10)), g1, t3_smooth_data(),
               outer_criterion = NULL)
  out <- paste(capture.output(print(summary(f))), collapse = "\n")
  expect_match(out, "lambda [start]", fixed = TRUE)
  expect_true(any(grepl("marked start", summary(f)@notes)))
  h <- statmod(y ~ s(x, bspline_smooth(k = 10), hyper = c(lambda = 2)), g1,
               t3_smooth_data(), outer_criterion = NULL)
  expect_match(paste(capture.output(print(summary(h))), collapse = "\n"),
               "lambda [fixed]", fixed = TRUE)
})

test_that("summary() shows the weight total beside the rows", {
  d <- t3_lin()
  f <- statmod(y ~ x, g1, d, weights = rep(2, 40))
  out <- capture.output(print(summary(f)))
  expect_true(any(grepl("prior weights summing to 80, not 40", out,
                        fixed = TRUE)))
  out1 <- capture.output(print(summary(statmod(y ~ x, g1, d))))
  expect_false(any(grepl("prior weights", out1)))
})

t3_delta_ref <- function(fit, p = "mu") {
  # the delta method on the predictor that statmod_eta() assembles, which
  # reads term_posterior() and term_levels() and carries no derivative,
  # differenced over the coefficients and the term's free parameters
  ns <- asNamespace("statmodels7")
  spec <- fit@spec
  de <- statmod_design(spec)
  st <- ns$statmod_structural_state(de)
  term <- names(st$zeta)[1]
  z0 <- st$zeta[[term]]
  free <- setdiff(names(z0), st$held[[term]])
  cf <- fit@coefficients
  b0 <- unlist(cf, use.names = FALSE)
  grp <- factor(rep(names(cf), lengths(cf)), levels = names(cf))
  f <- function(u) {
    cc <- stats::setNames(split(u[seq_along(b0)], grp), names(cf))
    z <- z0
    z[free] <- u[length(b0) + seq_along(free)]
    st$zeta[[term]] <- z
    on.exit(st$zeta[[term]] <- z0)
    as.numeric(ns$statmod_eta(spec, de, cc)$eta[[p]])
  }
  J <- numDeriv::jacobian(f, c(b0, z0[free]))
  V <- vcov(fit, readable = FALSE)
  ck <- unlist(lapply(names(cf), function(q) if (de[[q]]$npar)
    paste(q, de[[q]]$coef_names, sep = ":")))
  key <- c(ck, setdiff(rownames(V), ck))
  sqrt(pmax(rowSums((J %*% V[key, key]) * J), 0))
}

test_that("a prediction beside a regime term carries the posterior's motion", {
  skip_if_not_installed("numDeriv")
  set.seed(11)
  d <- data.frame(t = 1:200)
  st <- c(rep(1, 90), rep(2, 20), rep(1, 90))
  d$y <- c(-1, 2)[st] + rnorm(200, sd = 0.5)
  f <- statmod(y ~ regime(k = 2, time = t), g1, d)
  p <- predict(f, what = "mu", se = TRUE)
  ref <- t3_delta_ref(f)
  expect_lt(max(abs(p$se - ref) / ref), 1e-7)
  # the rows of the second state are read off its own level
  expect_gt(p$se[100], 2.5 * p$se[1])
})

test_that("a prediction beside a marginal break-point term does too", {
  skip_if_not_installed("numDeriv")
  set.seed(7)
  id <- rep(1:8, each = 10)
  x <- as.numeric(replicate(8, sort(runif(10, 0, 10))))
  pp <- 5 + rnorm(8, 0, 0.5)
  d <- data.frame(id = factor(id), x = x)
  d$y <- 1 + 2 * (x >= pp[id]) + rnorm(80, 0, 0.4)
  f <- statmod(y ~ jump(x, psi ~ random(~1 | id), marginal = TRUE), g1, d)
  p <- predict(f, what = "mu", se = TRUE)
  ref <- t3_delta_ref(f)
  expect_lt(max(abs(p$se - ref) / ref), 1e-7)
})

test_that("a filter that explodes is restarted whatever its flag says", {
  # on ubuntu's R 4.5.3 the first start of this fit stopped at -3299.1 with
  # the optimizer's flag TRUE and a Newton decrement of 3e51, and the second
  # start was never tried; the flag is forced TRUE here to reach that case on
  # every platform
  sb <- data.frame(month = seq_len(nrow(datasets::Seatbelts)),
                   share = datasets::Seatbelts[, "front"] /
                     (datasets::Seatbelts[, "front"] +
                        datasets::Seatbelts[, "rear"]),
                   law = datasets::Seatbelts[, "law"])
  ns <- asNamespace("statmodels7")
  orig <- ns$statmod_fit_joint
  forced <- function(...) {
    r <- orig(...)
    r$converged <- TRUE
    r
  }
  utils::assignInNamespace("statmod_fit_joint", forced, "statmodels7")
  on.exit(utils::assignInNamespace("statmod_fit_joint", orig, "statmodels7"))
  fit <- statmod(share ~ law + gas(p = 1, q = 1, time = month),
                 distributions7::beta1_distrib(), sb)
  expect_equal(as.numeric(logLik(fit)), 450.4818421, tolerance = 1e-6)
})
