# statmod(): the front end, the alternation and what a fit exposes.

set.seed(1)
n <- 120
dd <- data.frame(x = runif(n), z = runif(n),
                 g = factor(rep(letters[1:4], length.out = n)))
dd$y <- 1 + 2 * dd$x + rnorm(n, sd = 0.6)
dd$R <- matrix(rnorm(n * 3), n, 3)
dd$L <- matrix(rnorm(n * 4), n, 4)

test_that("a gaussian fit is least squares", {
  fit <- statmod(y ~ x + g, distributions7::gaussian1_distrib(), dd)
  ref <- stats::lm(y ~ x + g, dd)
  expect_true(fit@converged)
  expect_equal(fit@coefficients$mu, unname(stats::coef(ref)), tolerance = 1e-10)
  expect_equal(fit@loglik, as.numeric(stats::logLik(ref)), tolerance = 1e-10)
})

test_that("loglik() and logLik() agree to the last digit", {
  # the cheapest check that the callable route and the fitting route are the
  # same model: one rebuilds the design from the specification, the other
  # reports what the fit reached
  fit <- statmod(y ~ x + g, distributions7::gaussian1_distrib(), dd)
  expect_identical(loglik(fit), as.numeric(stats::logLik(fit)))
})

test_that("the model is callable at other parameters and other data", {
  fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
  at_zero <- loglik(fit, par = list(mu = c(0, 0), sigma = 0))
  expect_lt(at_zero, loglik(fit))
  # on a subset, the log-likelihood is a sum over fewer observations
  half <- loglik(fit, data = dd[1:60, ])
  expect_gt(half, loglik(fit))
  expect_true(is.finite(half))
  # the gradient at the optimum is small, to the tolerance the run promised
  expect_lt(max(abs(unlist(gradient(fit)))), 1e-3)
  H <- hessian(fit)
  expect_identical(dim(H), c(3L, 3L))
})

test_that("what loglik() refuses", {
  fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
  expect_error(loglik(fit, par = list(wrong = 1)), "not a parameter")
  expect_error(loglik(fit, par = list(mu = 1)), "length 1")
})

test_that("every parameter can be modelled", {
  fit <- statmod(y ~ x | sigma ~ z, distributions7::gaussian1_distrib(), dd)
  expect_true(fit@converged)
  expect_length(fit@coefficients$mu, 2L)
  expect_length(fit@coefficients$sigma, 2L)
  # modelling the scale cannot make the likelihood worse than holding it fixed
  flat <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
  expect_gte(fit@loglik, flat@loglik - 1e-8)
})

test_that("a differentiable penalty joins the smooth block", {
  fit <- statmod(y ~ x + ridge(R), distributions7::gaussian1_distrib(), dd)
  expect_true(fit@converged)
  # one block, so one sweep and no alternation
  expect_identical(max(fit@history$blocks$sweep), 1L)
  expect_identical(unique(fit@history$blocks$block), "smooth")
  # and the penalized coefficients are smaller than the unpenalized fit's
  unpen <- unname(stats::coef(stats::lm(y ~ x + R, dd)))[3:5]
  design <- statmod_design(fit@spec)
  cols <- design$mu$blocks[["ridge(R)"]]
  expect_lt(sum(fit@coefficients$mu[cols]^2), sum(unpen^2))
})

test_that("a penalty with a kink is fitted apart, and the fit alternates", {
  fit <- statmod(y ~ x + lasso(L), distributions7::gaussian1_distrib(), dd)
  expect_true(fit@converged)
  h <- fit@history$blocks
  expect_true(all(c("smooth", "mu/lasso(L)") %in% h$block))
  expect_gt(max(h$sweep), 1L)
  # the alternation descends: every block leaves the objective no higher
  expect_true(all(h$change > -1e-8))
  # and the objective at the end is the one the fit reports
  expect_equal(min(h$objective), fit@objective, tolerance = 1e-10)
})

test_that("the block split reads the kink and not the term's name", {
  spec <- statmod_spec(y ~ x + ridge(R) + lasso(L),
                       distributions7::gaussian1_distrib(), dd)
  design <- statmod_design(spec)
  b <- statmod_blocks(spec, design)
  expect_length(b$sparse, 1L)
  expect_identical(b$sparse[[1L]]$term, "lasso(L)")
  # the ridge stays in the smooth block with the unpenalized coefficients
  expect_true(all(design$mu$blocks[["ridge(R)"]] %in% b$smooth))
  expect_false(any(design$mu$blocks[["lasso(L)"]] %in% b$smooth))
})

test_that("an optimizer can replace the scoring step", {
  fit <- statmod(y ~ x + g, distributions7::gaussian1_distrib(), dd,
                 inner_optimizer = optimizers7::bfgs())
  ref <- unname(stats::coef(stats::lm(y ~ x + g, dd)))
  expect_equal(fit@coefficients$mu, ref, tolerance = 1e-5)
})

test_that("the elapsed time is reported in the unit it deserves", {
  expect_identical(format_duration(0), "0 s")
  expect_identical(format_duration(3.4e-4), "340 us")
  expect_identical(format_duration(0.25), "250 ms")
  expect_identical(format_duration(1.5), "1.5 s")
  expect_identical(format_duration(90), "1.5 min")
  expect_identical(format_duration(7200), "2 h")
  expect_identical(format_duration(200000), "2.31 d")
  expect_identical(format_duration(NA_real_), "NA")
})

test_that("verbosity names the loops", {
  expect_identical(statmodels7:::verbosity(0),
                   list(outer = FALSE, blocks = FALSE, inner = FALSE,
                        optimizer = FALSE))
  expect_identical(statmodels7:::verbosity(2),
                   list(outer = TRUE, blocks = TRUE, inner = TRUE,
                        optimizer = FALSE))
  expect_identical(statmodels7:::verbosity(c(blocks = TRUE, inner = FALSE)),
                   list(outer = FALSE, blocks = TRUE, inner = FALSE,
                        optimizer = FALSE))
  expect_error(statmodels7:::verbosity(c(wrong = TRUE)), "'blocks'")
  expect_error(statmodels7:::verbosity("loud"), "number from 0 to 3")
  # and it actually prints
  expect_output(statmod(y ~ x, distributions7::gaussian1_distrib(), dd,
                        verbose = 1), "sweep 1")
})

test_that("print says what the fit is", {
  fit <- statmod(y ~ x | sigma ~ z, distributions7::gaussian1_distrib(), dd)
  out <- utils::capture.output(print(fit))
  expect_true(any(grepl("gaussian1", out)))
  expect_true(any(grepl("mu .*~ x", out)))
  expect_true(any(grepl("sigma .*~ z", out)))
  expect_true(any(grepl("log-likelihood", out)))
  expect_true(any(grepl("converged", out)))
})

test_that("prior weights change the fit the way they should", {
  w <- rep(c(1, 3), length.out = n)
  fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd, weights = w)
  ref <- stats::lm(y ~ x, dd, weights = w)
  expect_equal(fit@coefficients$mu, unname(stats::coef(ref)), tolerance = 1e-8)
  # and print says the weights do not sum to n, which is the honest reading
  expect_output(print(fit), "prior weights summing to")
})

test_that("start is validated and used", {
  expect_error(statmod(y ~ x, distributions7::gaussian1_distrib(), dd,
                       start = list(wrong = 1)), "not a parameter")
  expect_error(statmod(y ~ x, distributions7::gaussian1_distrib(), dd,
                       start = list(mu = 1)), "length 1")
  fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd,
                 start = list(mu = c(0, 0), sigma = 0))
  expect_true(fit@converged)
})

test_that("coef names the coefficients", {
  fit <- statmod(y ~ x + g, distributions7::gaussian1_distrib(), dd)
  cf <- stats::coef(fit)
  expect_named(cf, c("mu", "sigma"))
  expect_identical(names(cf$mu)[1L], "(Intercept)")
})

test_that("hyper sets a hyperparameter by key or by label", {
  # a term is keyed by its call as written, which is what keeps two lassos on
  # different covariates apart, and named by its label, which is what anybody
  # would type
  f <- y ~ x + lasso(~ z)
  fit_key <- statmod(f, distributions7::gaussian1_distrib(), dd,
                     hyper = list(mu = list(`lasso(~z)` = c(lambda = 50))))
  fit_lab <- statmod(f, distributions7::gaussian1_distrib(), dd,
                     hyper = list(mu = list(lasso = c(lambda = 50))))
  expect_equal(fit_key@coefficients$mu, fit_lab@coefficients$mu,
               tolerance = 1e-12)
  expect_equal(unname(fit_lab@hyper$mu[[1L]][["lambda"]]), 50)

  # a larger penalty shrinks harder, which is the property a hyperparameter is
  # for and the cheapest evidence that the value reached the penalty at all
  weak <- statmod(f, distributions7::gaussian1_distrib(), dd,
                  hyper = list(mu = list(lasso = c(lambda = 1))))
  strong <- statmod(f, distributions7::gaussian1_distrib(), dd,
                    hyper = list(mu = list(lasso = c(lambda = 500))))
  expect_lt(abs(strong@coefficients$mu[3L]), abs(weak@coefficients$mu[3L]))
})

test_that("hyper is validated and says what is available", {
  f <- y ~ x + lasso(~ z)
  expect_error(statmod(f, distributions7::gaussian1_distrib(), dd,
                       hyper = list(wrong = list(lasso = 1))),
               "not a parameter")
  expect_error(statmod(f, distributions7::gaussian1_distrib(), dd,
                       hyper = list(mu = list(ridge = 1))),
               "names no penalized term")
  expect_error(statmod(f, distributions7::gaussian1_distrib(), dd,
                       hyper = list(mu = list(lasso = c(nu = 1)))),
               "hyperparameters are")
  expect_error(statmod(f, distributions7::gaussian1_distrib(), dd,
                       hyper = list(mu = list(lasso = c(1, 2)))),
               "has length 2")
  # two terms sharing a label is ambiguous, and guessing would set the
  # hyperparameter on the wrong block
  expect_error(statmod(y ~ lasso(~ x) + lasso(~ z),
                       distributions7::gaussian1_distrib(), dd,
                       hyper = list(mu = list(lasso = c(lambda = 5)))),
               "ambiguous")
})


test_that("a smooth term reports an edf between its null space and its rank", {
  # edf() takes the term's unpenalized curvature as its THIRD argument and the
  # hyperparameters as its fourth; a positional call put the second where the
  # first belonged, every smooth term reported NA, and the degrees of freedom
  # then counted the unpenalized terms alone
  set.seed(9)
  n2 <- 300
  ds <- data.frame(x = runif(n2, -2, 2))
  ds$y <- sin(1.4 * ds$x) + stats::rnorm(n2, sd = 0.3)
  fit <- statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), ds)

  e <- fit@edf
  sm <- e[e$term != "linpar" & e$parameter == "mu", , drop = FALSE]
  expect_identical(nrow(sm), 1L)
  expect_false(is.na(sm$edf))
  # the Demmler-Reinsch penalty is rank deficient by exactly one, so the edf
  # runs from the coefficient count down to one and never outside it
  expect_gt(sm$edf, 1 - 1e-8)
  expect_lt(sm$edf, sm$coefficients + 1e-8)

  # a heavier penalty spends less, which is what the number is for
  hard <- statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), ds,
                  hyper = list(mu = list(`s(x)` = c(lambda = 1e6))))
  hs <- hard@edf
  expect_lt(hs$edf[hs$term != "linpar" & hs$parameter == "mu"], sm$edf)

  # and the criteria are built on the total, so they see the smooth
  expect_equal(attr(stats::logLik(fit), "df"), sum(e$edf), tolerance = 1e-12)
  expect_gt(attr(stats::logLik(fit), "df"), 3)
})

test_that("the budget and the tolerance are read off the method", {
  # statmod() carries neither: an argument accepted and ignored is worse than
  # one that errors, and a caller setting iwls(maxit = 20) and a loose
  # maxit = 100 would get one of them with nothing said about the other
  fm <- names(formals(statmod))
  expect_false("maxit" %in% fm)
  expect_false("tol" %in% fm)

  b <- method_budget(iwls(maxit = 20, tol = 1e-4))
  expect_identical(b$maxit, 20L)
  expect_equal(b$tol, 1e-4)

  # an optimizer carries maxit and a criterion, and the tolerance is the
  # largest one that criterion contains: a combined rule stops at whichever
  # part fires first, so asking for more than the loop inside can deliver
  # would be asking the alternation to outlive it
  o <- optimizers7::bfgs(maxit = 33,
                         criterion = optimizers7::crit_grad(1e-5))
  expect_identical(method_budget(o)$maxit, 33L)
  expect_equal(method_budget(o)$tol, 1e-5)

  cmb <- optimizers7::crit_any(optimizers7::crit_grad(1e-9),
                              optimizers7::crit_rel_obj(1e-3))
  expect_equal(criterion_tol(cmb), 1e-3)

  expect_error(statmod(y ~ x, distributions7::gaussian1_distrib(), dd,
                       inner_optimizer = "bfgs"),
               "an optimizers7 optimizer", fixed = TRUE)
})

test_that("the alternation obeys the method's budget", {
  # one sweep is not enough for a lasso to settle beside the smooth block, so
  # a budget of one must come back not converged rather than claiming it
  dl <- dd
  dl$noise <- stats::rnorm(nrow(dl))
  one <- statmod(y ~ x + lasso(~ noise),
                 distributions7::gaussian1_distrib(), dl,
                 inner_optimizer = iwls(maxit = 1L),
                 hyper = list(mu = list(lasso = c(lambda = 5))))
  expect_identical(max(one@history$blocks$sweep), 1L)
  expect_false(one@converged)

  many <- statmod(y ~ x + lasso(~ noise),
                  distributions7::gaussian1_distrib(), dl,
                  hyper = list(mu = list(lasso = c(lambda = 5))))
  expect_true(many@converged)
  expect_gt(max(many@history$blocks$sweep), 1L)
})
