# Estimating the hyperparameters by a marginal criterion.

set.seed(21)
n <- 300
ds <- data.frame(x = runif(n, -2, 2))
ds$y <- sin(1.4 * ds$x) + stats::rnorm(n, sd = 0.3)

test_that("the criterion is the Laplace formula, assembled independently", {
  skip_if_not_installed("numDeriv")
  # the observed information is asked for so that numDeriv's Hessian of the
  # penalized objective is the same matrix the criterion uses
  fit <- statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), ds,
                 outer_criterion = reml(hessian = "observed"))
  spec <- fit@spec
  design <- statmod_design(spec)
  obj <- statmod_objective(spec, fit@hyper, design, expected = FALSE)
  b <- unlist(fit@coefficients[spec@distrib@params], use.names = FALSE)

  ll <- as.numeric(stats::logLik(fit))
  rho <- statmod_penalty_at(spec, fit@coefficients, fit@hyper, design, "value")
  M <- numDeriv::hessian(obj$fn, b)
  by_hand <- ll - rho + length(b) / 2 * log(2 * pi) -
    determinant(M, logarithm = TRUE)$modulus[[1L]] / 2

  m <- statmod_marginal(spec, design, fit@coefficients, fit@hyper,
                        reml(hessian = "observed"))
  expect_equal(m$value, as.numeric(by_hand), tolerance = 1e-5)
  expect_equal(m$value, fit@criterion, tolerance = 1e-8)
  expect_identical(m$q, length(b))
})

test_that("REML is Wood's criterion, reached by the other route", {
  # ell - rho + (p/2) log 2pi - (1/2) log|H+S| must equal
  # ell - (lambda/2) b'Pb + (1/2) log|lambda P|_+ - (1/2) log|H+S|
  #      + (M_p/2) log 2pi,   M_p = p - rank
  # the two differing only by how the penalty's normalizing constant is
  # written. Keeping that constant is what makes the first form correct with
  # nothing added by hand, so the identity is worth asserting rather than
  # assuming. The sign of the last term is the whole content of the check:
  # written the other way round the two forms are apart by exactly
  # (p - r) log 2pi, which here is the intercept, the smooth's linear
  # component and the scale's intercept.
  fit <- statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), ds,
                 outer_criterion = reml())
  spec <- fit@spec
  design <- statmod_design(spec)
  nm <- names(spec@terms$mu)[vapply(spec@terms$mu, function(t)
    !is.null(modelterms7::term_penalty(t)), logical(1))]
  pen <- modelterms7::term_penalty(spec@terms$mu[[nm]])
  th <- as.list(fit@hyper$mu[[nm]])
  bt <- fit@coefficients$mu[design$mu$blocks[[nm]]]

  P <- penalties7::penalty_matrix(pen, th)
  quad <- as.numeric(crossprod(bt, P %*% bt)) / 2
  lpd <- penalties7::penalty_logpdet(pen, th)$value
  r <- penalties7::penalty_rank(pen)
  p_tot <- length(unlist(fit@coefficients[spec@distrib@params]))

  # the EXPECTED information on both sides: the criterion is asked for the
  # one the reference below builds, rather than the default, which is the
  # observed one since it buys the exact outer gradient
  m <- statmod_marginal(spec, design, fit@coefficients, fit@hyper,
                        reml("expected"))
  H <- statmod_information_at(spec, fit@coefficients, design, TRUE)
  S <- statmod_penalty_at(spec, fit@coefficients, fit@hyper, design, "hessian")
  logdet <- determinant(H + S, logarithm = TRUE)$modulus[[1L]]
  wood <- as.numeric(stats::logLik(fit)) - quad + lpd / 2 - logdet / 2 +
    (p_tot - r) / 2 * log(2 * pi)

  expect_equal(m$value, wood, tolerance = 1e-8)
  # and the count of unpenalized directions is what that term carries
  expect_equal(p_tot - r, 3)
})

test_that("the reported hyperparameter is where the criterion is best", {
  fit <- statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), ds,
                 outer_criterion = reml())
  nm <- names(fit@hyper$mu)[1L]
  lam <- fit@hyper$mu[[nm]][["lambda"]]
  at <- function(v) {
    h <- fit@hyper
    h$mu[[nm]][["lambda"]] <- v
    statmod(y ~ s(x, k = 10, lambda = v),
            distributions7::gaussian1_distrib(), ds)
  }
  crit <- function(f) {
    statmod_marginal(f@spec, statmod_design(f@spec), f@coefficients,
                     f@hyper, reml())$value
  }
  best <- fit@criterion
  expect_gt(best, crit(at(lam * 3)))
  expect_gt(best, crit(at(lam / 3)))
  expect_gt(best, crit(at(lam * 1.4)))
  expect_gt(best, crit(at(lam / 1.4)))
})

test_that("a smoothing parameter estimated by REML lands on a sane edf", {
  fit <- statmod(y ~ s(x, k = 12), distributions7::gaussian1_distrib(), ds,
                 outer_criterion = reml())
  e <- fit@edf$edf[fit@edf$term != "linpar" & fit@edf$parameter == "mu"]
  # a sine over four periods of the covariate is neither a straight line nor
  # an interpolation of 300 points
  expect_gt(e, 3)
  expect_lt(e, 11)
  # and the fit is close to the truth it was drawn from
  expect_lt(sqrt(mean((predict(fit, "mu") - sin(1.4 * ds$x))^2)), 0.06)
  expect_true(is.finite(fit@criterion))
  expect_gt(nrow(fit@history$outer), 3)
})

test_that("ML puts less variance on a random effect than REML", {
  # the classic downward bias: profiling the fixed effects rather than
  # integrating them leaves the variance component too small
  set.seed(22)
  m <- 25
  g <- factor(rep(paste0("g", seq_len(m)), each = 10))
  u <- stats::rnorm(m, sd = 0.8)
  dv <- data.frame(g = g, x = runif(250))
  dv$y <- 1 + 2 * dv$x + u[as.integer(g)] + stats::rnorm(250, sd = 0.5)

  r <- statmod(y ~ x + random(~ 1 | g),
               distributions7::gaussian1_distrib(), dv, outer_criterion = reml())
  l <- statmod(y ~ x + random(~ 1 | g),
               distributions7::gaussian1_distrib(), dv, outer_criterion = ml())
  nm <- "random(~1 | g)"
  # the effects' penalty carries the PRECISION, so the variance component
  # a reader wants is 1/sqrt(lambda) 
  sd_r <- 1 / sqrt(r@hyper$mu[[nm]][[1L]]) 
  sd_l <- 1 / sqrt(l@hyper$mu[[nm]][[1L]])

  expect_lt(sd_l, sd_r)
  # both near the 0.8 they were drawn from
  expect_gt(sd_r, 0.5)
  expect_lt(sd_r, 1.2)
  # and the fixed effects are recovered
  expect_equal(r@coefficients$mu[2L], 2, tolerance = 0.25)
})

test_that("a kinked penalty keeps the hyperparameter it was given", {
  # a Laplace approximation asks for a second derivative, and a lasso has none
  # at zero, so its lambda is not something a marginal criterion can estimate
  set.seed(23)
  dl <- ds
  for (j in 1:3) dl[[paste0("n", j)]] <- stats::rnorm(n)
  fit <- statmod(y ~ s(x, k = 10) + lasso(~ n1 + n2 + n3, lambda = 40),
                 distributions7::gaussian1_distrib(), dl,
                 outer_criterion = reml())
  ln <- grep("^lasso", names(fit@hyper$mu), value = TRUE)
  expect_equal(unname(fit@hyper$mu[[ln]][["lambda"]]), 40)
  # while the smooth's did move
  expect_false(isTRUE(all.equal(
    unname(fit@hyper$mu[["s(x, k = 10)"]][["lambda"]]), 1)))
})

test_that("the criterion applies to a smooth penalty and to nothing else", {
  # It comes into play if and only if the model carries one, which is a
  # property of the MODEL and not of how the argument was written: a model
  # with no penalty fits and estimates nothing, whether or not the default
  # was typed out.
  a <- statmod(y ~ x, distributions7::gaussian1_distrib(), ds)
  b <- statmod(y ~ x, distributions7::gaussian1_distrib(), ds,
               outer_criterion = reml())
  expect_equal(a@coefficients, b@coefficients)
  expect_length(unlist(a@hyper), 0L)
  expect_true(is.na(a@criterion))

  # a kinked penalty is not a marginal criterion's business: it is chosen by
  # `sparse_criterion`, a path over its own values, which is bic() by default
  k <- statmod(y ~ lasso(~ x), distributions7::gaussian1_distrib(), ds)
  expect_false(isTRUE(all.equal(unname(unlist(k@hyper)), 1)))
  expect_true(is.finite(k@criterion))
  # and holding it is asking for that
  h <- statmod(y ~ lasso(~ x), distributions7::gaussian1_distrib(), ds,
               sparse_criterion = NULL)
  expect_equal(unname(unlist(h@hyper)), 1)

  # and where there IS a smooth penalty it acts, without being asked
  m <- statmod(y ~ s(x, k = 8), distributions7::gaussian1_distrib(), ds)
  expect_false(isTRUE(all.equal(
    unname(m@hyper$mu[["s(x, k = 8)"]][["lambda"]]), 1)))
  expect_false(is.na(m@criterion))

  expect_error(statmod(y ~ s(x, k = 8), distributions7::gaussian1_distrib(),
                       ds, outer_criterion = "reml"),
               "reml(), ml(), aic(), bic(), cv() or NULL", fixed = TRUE)
})

test_that("ml refuses a penalty whose null space it cannot read", {
  # an anisotropic tensor smooth carries an additive penalty, which reports a
  # rank but exposes no null basis; guessing which directions are profiled
  # would give a criterion that is arithmetic without a meaning
  set.seed(24)
  dt <- data.frame(x1 = runif(200, -1, 1), x2 = runif(200, -1, 1))
  dt$y <- dt$x1^2 + dt$x2 + stats::rnorm(200, sd = 0.3)
  expect_error(statmod(y ~ te(x1, x2, k = 4),
                       distributions7::gaussian1_distrib(), dt,
                       outer_criterion = ml()),
               "cannot read the null space")
  # while reml integrates everything and needs no such basis
  fit <- statmod(y ~ te(x1, x2, k = 4),
                 distributions7::gaussian1_distrib(), dt,
                 outer_criterion = reml())
  expect_true(is.finite(fit@criterion))
  # one smoothing parameter per margin, both moved off the probe value
  th <- fit@hyper$mu[["te(x1, x2, k = 4)"]]
  expect_gte(length(th), 2L)
})

test_that("the summary marks an estimated hyperparameter as estimated", {
  fit <- statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), ds,
                 outer_criterion = reml())
  s <- summary(fit)
  kinds <- vapply(s@tables$mu, `[[`, character(1), "kind")
  b <- s@tables$mu[[which(kinds == "smooth")]]
  expect_identical(b$table$role[b$table$name == "lambda"], "estimated")
  # and it now carries a standard error and an interval, so what says the
  # criterion estimated it is the note rather than a mark in the cell
  r <- b$table[b$table$name == "lambda", , drop = FALSE]
  expect_identical(r$source, "reml")
  expect_true(is.finite(r$se))
  expect_true(any(grepl("REML", s@notes)))
  expect_output(print(fit), "REML")
})

test_that("verbose names the outer loop among its switches", {
  expect_true(verbosity(1)$outer)
  expect_false(verbosity(0)$outer)
  expect_true(verbosity(c(outer = TRUE))$outer)
  expect_false(verbosity(c(outer = TRUE))$blocks)
  expect_error(verbosity(c(wrong = TRUE)), "'outer'", fixed = TRUE)
})


test_that("an unavailable point is a barrier the search steps back from", {
  # It must be FINITE. An infinite one is useless to a method that
  # differences its own gradient -- the probe lands in the unavailable
  # region, the difference is non-finite and the search stops -- which is
  # what ended a panel fit at its 73rd evaluation with the criterion still
  # improving.
  set.seed(41)
  n <- 300
  dv <- data.frame(x = stats::runif(n, -2, 2))
  dv$y <- sin(1.5 * dv$x) + stats::rnorm(n, sd = 0.3)
  spec <- statmod_spec(y ~ s(x, k = 10), distributions7::gaussian1_distrib(),
                       dv)
  design <- statmod_design(spec)
  blocks <- statmod_blocks(spec, design)
  hyper <- statmod_hyper_start(spec)
  obj <- statmod_objective(spec, hyper, design)
  beta <- statmodels7:::statmod_start(spec, design, obj, NULL)

  res <- statmodels7:::outer_fit(spec, design, blocks, hyper, iwls(),
                                 reml(), NULL, beta, "bartlett", 100L, 1e-6,
                                 statmodels7:::verbosity(0))
  expect_true(is.finite(res$criterion))
  # the search reports which optimizer ran, defaults included
  expect_true(S7::S7_inherits(res$optimizer, optimizers7::optimizer))
})

test_that("the outer search uses the optimizer it was given", {
  set.seed(42)
  n <- 300
  dv <- data.frame(x = stats::runif(n, -2, 2))
  dv$y <- sin(1.5 * dv$x) + stats::rnorm(n, sd = 0.3)
  fml <- y ~ s(x, k = 8)
  for (o in list(optimizers7::nelder_mead(), optimizers7::bfgs(),
                 optimizers7::lbfgs())) {
    f <- statmod(fml, distributions7::gaussian1_distrib(), dv,
                 outer_optimizer = o)
    expect_identical(f@methods$search@name, o@name)
  }
  # and with none given it records the one it chose, which is the question a
  # reader of a fit asks
  f0 <- statmod(fml, distributions7::gaussian1_distrib(), dv)
  expect_true(S7::S7_inherits(f0@methods$search, optimizers7::optimizer))
})

test_that("a trace names the term without repeating its specification", {
  long <- paste0("gas(p = 1, q = 1, time = t, by = ~ridge(~id), links = ",
                 "list(omega = identity_link()))::omega::ridge(~id)")
  short <- statmodels7:::short_keys(long)
  expect_lt(nchar(short), nchar(long) / 2)
  # what distinguishes one entry of a term from another is kept whole
  expect_true(endsWith(short, "::omega::ridge(~id)"))
  expect_true(startsWith(short, "gas(p = 1, ...)"))
  # the first argument survives, so two smooths stay apart
  expect_identical(statmodels7:::short_keys(c("s(x, k = 20)", "s(z, k = 8)")),
                   c("s(x, ...)", "s(z, ...)"))
  # and where shortening would collide, nothing is shortened
  expect_identical(statmodels7:::short_keys(c("s(x, k = 20)", "s(x, k = 8)")),
                   c("s(x, k = 20)", "s(x, k = 8)"))
  # a call with one argument is already short
  expect_identical(statmodels7:::short_keys("random(~1 | g)"), "random(~1 | g)")
})

test_that("a hyperparameter chosen by a path is marked estimated, not held", {
  # THE DEFECT: the summary asked `methods$outer`, which carries the MARGINAL
  # criterion alone, so a lambda a path had chosen was reported as held at a
  # value the caller had given -- of a number the caller never saw.
  set.seed(41)
  n2 <- 300
  ds <- data.frame(x = runif(n2))
  for (j in 1:6) ds[[paste0("n", j)]] <- stats::rnorm(n2)
  ds$y <- 1 + 2 * ds$x + 1.5 * ds$n1 + stats::rnorm(n2, sd = 0.4)
  fit <- statmod(y ~ x + lasso(~ n1 + n2 + n3 + n4 + n5 + n6),
                 distributions7::gaussian1_distrib(), ds,
                 sparse_criterion = bic())
  s <- summary(fit)
  kinds <- vapply(s@tables$mu, `[[`, character(1), "kind")
  b <- s@tables$mu[[which(kinds == "selection")]]
  r <- b$table[b$table$name == "lambda", , drop = FALSE]
  expect_identical(r$role, "estimated")
  expect_identical(r$source, "bic")
  # and it is not the value it came in with
  expect_false(isTRUE(all.equal(r$estimate, 1)))
  expect_output(print(s), "(bic)", fixed = TRUE)
  expect_false(any(grepl("held at the value", s@notes)))
})

test_that("the hyperparameters head their block", {
  set.seed(42)
  n2 <- 200
  ds <- data.frame(x = runif(n2))
  ds$y <- sin(4 * ds$x) + stats::rnorm(n2, sd = 0.3)
  s <- summary(statmod(y ~ s(x, k = 8), distributions7::gaussian1_distrib(),
                       ds, outer_criterion = reml()))
  kinds <- vapply(s@tables$mu, `[[`, character(1), "kind")
  b <- s@tables$mu[[which(kinds == "smooth")]]
  expect_identical(b$table$name[[1L]], "lambda")
})

test_that("a hyperparameter estimated by REML carries a standard error", {
  set.seed(43)
  n2 <- 200
  ds <- data.frame(x = runif(n2))
  ds$y <- sin(6 * ds$x) + stats::rnorm(n2, sd = 0.3)
  fit <- statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), ds,
                 outer_criterion = reml())
  sp <- fit@spec
  de <- statmod_design(sp)
  V <- statmod_hyper_vcov(sp, de, fit@coefficients, fit@hyper,
                          fit@methods$outer)
  expect_false(is.null(V))
  se_eta <- sqrt(diag(as.matrix(V)))
  expect_true(all(is.finite(se_eta) & se_eta > 0))

  # the reference shares no arithmetic with the exact Hessian: the criterion
  # itself, the mode refitted at each probe, differenced twice on the free
  # scale. The step is 0.05 because a second difference amplifies the
  # criterion's own rounding by h^-2 -- at 1e-3 the same comparison reads
  # 0.40 against 0.59 and measures the arithmetic rather than the curvature.
  blk <- statmod_blocks(sp, de)
  idx <- outer_hyper_index(sp, blk)
  eta0 <- hyper_to_eta(fit@hyper, idx)
  basis <- integrated_basis(sp, de, "reml")
  crit <- function(e) {
    hy <- eta_to_hyper(e, idx, fit@hyper)
    r <- statmod_alternate(sp, de, blk, hy, iwls(),
                           unlist(fit@coefficients, use.names = FALSE),
                           TRUE, "bartlett", 200L, 1e-6, verbosity(0))
    statmod_marginal(sp, de, r$obj$split(r$par), hy, fit@methods$outer,
                     "bartlett", basis)$value
  }
  h <- 0.05
  H <- (crit(eta0 - h) - 2 * crit(eta0) + crit(eta0 + h)) / h^2
  expect_equal(se_eta[[1L]], sqrt(-1 / H), tolerance = 1e-3)

  # and the summary reports it beside the estimate, with the interval built
  # on that scale and mapped back, so a positive quantity keeps a positive
  # lower end and no test is printed
  s <- summary(fit)
  kinds <- vapply(s@tables$mu, `[[`, character(1), "kind")
  b <- s@tables$mu[[which(kinds == "smooth")]]
  r <- b$table[b$table$name == "lambda", , drop = FALSE]
  expect_true(is.finite(r$se) && r$se > 0)
  expect_gt(r$lower, 0)
  expect_lt(r$lower, r$estimate)
  expect_gt(r$upper, r$estimate)
  expect_true(is.na(r$statistic))
  expect_identical(r$source, "reml")
})

test_that("a path that cannot score a single point says so", {
  # the matrix a penalized term is built from lives in the calling
  # environment rather than in `data`, so no fold can be rebuilt. Every
  # deviance came back NA, the path chose nothing, and the fit returned the
  # DEFAULT hyperparameter reporting success.
  set.seed(44)
  n2 <- 120
  Z <- matrix(stats::rnorm(n2 * 5), n2, 5)
  colnames(Z) <- paste0("z", 1:5)
  ds <- data.frame(Z = Z, y = as.numeric(Z %*% c(2, 0, 0, 1, 0)) +
                     stats::rnorm(n2, sd = 0.4))
  expect_error(
    statmod(y ~ 1 + lasso(Z), distributions7::gaussian1_distrib(), ds,
            sparse_criterion = cv()),
    "could not score a single point")

  # the same model with the matrix as a column of `data` fits and chooses
  ds2 <- data.frame(y = ds$y)
  ds2$Z <- Z
  fit <- statmod(y ~ 1 + lasso(Z), distributions7::gaussian1_distrib(), ds2,
                 sparse_criterion = cv())
  expect_false(isTRUE(all.equal(
    unname(fit@hyper$mu[["lasso(Z)"]][["lambda"]]), 1)))
})
