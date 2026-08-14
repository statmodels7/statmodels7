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
  # one block, so one pass and no alternation
  expect_identical(max(fit@history$blocks$pass), 1L)
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
  expect_gt(max(h$pass), 1L)
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
                        verbose = 1), "pass 1")
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

test_that("a term holds its own hyperparameter, and it reaches the fit", {
  # WHICH hyperparameters are held is said where the penalty is named. There
  # is no second place to say it, so there is no key to get right and no
  # ambiguity between two terms sharing a label.
  fit <- statmod(y ~ x + lasso(~ z, lambda = 50),
                 distributions7::gaussian1_distrib(), dd)
  expect_equal(unname(fit@hyper$mu[[1L]][["lambda"]]), 50)

  # a larger penalty shrinks harder, which is the property a hyperparameter is
  # for and the cheapest evidence that the value reached the penalty at all
  weak <- statmod(y ~ x + lasso(~ z, lambda = 1),
                  distributions7::gaussian1_distrib(), dd)
  strong <- statmod(y ~ x + lasso(~ z, lambda = 500),
                    distributions7::gaussian1_distrib(), dd)
  expect_lt(abs(strong@coefficients$mu[3L]), abs(weak@coefficients$mu[3L]))

  # two terms of the same label are two terms with their own values, which the
  # old keying could not express without asking which was which
  two <- statmod(y ~ lasso(~ x, lambda = 1) + lasso(~ z, lambda = 500),
                 distributions7::gaussian1_distrib(), dd)
  expect_length(two@hyper$mu, 2L)
  expect_equal(unname(unlist(lapply(two@hyper$mu, `[[`, "lambda"))),
               c(1, 500))
})

test_that("hyper is gone from statmod() and says where to write it", {
  # two arguments saying the same thing, one of them read by nobody whenever
  # they disagreed
  expect_error(statmod(y ~ x + lasso(~ z),
                       distributions7::gaussian1_distrib(), dd,
                       hyper = list(mu = list(lasso = c(lambda = 5)))),
               "has been removed")
  expect_error(statmod(y ~ x, distributions7::gaussian1_distrib(), dd,
                       nonsense = 1), "unused argument")
  # and the message names the spelling that works
  msg <- tryCatch(statmod(y ~ x + lasso(~ z),
                          distributions7::gaussian1_distrib(), dd,
                          hyper = list(mu = list(lasso = 1))),
                  error = function(e) conditionMessage(e))
  expect_match(msg, "lambda = 3", fixed = TRUE)
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
  hard <- statmod(y ~ s(x, k = 10, lambda = 1e6),
                  distributions7::gaussian1_distrib(), ds)
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
  one <- statmod(y ~ x + lasso(~ noise, lambda = 5),
                 distributions7::gaussian1_distrib(), dl,
                 inner_optimizer = iwls(maxit = 1L))
  expect_identical(max(one@history$blocks$pass), 1L)
  expect_false(one@converged)

  many <- statmod(y ~ x + lasso(~ noise, lambda = 5),
                  distributions7::gaussian1_distrib(), dl)
  expect_true(many@converged)
  expect_gt(max(many@history$blocks$pass), 1L)
})

test_that("the term decides which hyperparameters are estimated", {
  # THE RULE: NULL, the default, means estimated; a value means held. It is
  # the term's answer because the term is where the penalty is named.
  set.seed(61)
  n2 <- 150
  Z <- matrix(stats::rnorm(n2 * 6), n2, 6)
  colnames(Z) <- paste0("z", 1:6)
  de <- data.frame(y = as.numeric(Z %*% c(2, -1.5, rep(0, 4))) +
                     stats::rnorm(n2, sd = 0.5))
  de$Z <- Z
  hy <- function(f) unlist(f@hyper$mu[[1L]])

  # both estimated, neither at the value it started from
  free <- statmod(y ~ enet(Z), distributions7::gaussian1_distrib(), de)
  expect_false(isTRUE(all.equal(unname(hy(free)[["alpha"]]), 0.5)))

  # one held, one estimated
  half <- statmod(y ~ enet(Z, alpha = 0.5),
                  distributions7::gaussian1_distrib(), de)
  expect_equal(unname(hy(half)[["alpha"]]), 0.5)
  expect_false(isTRUE(all.equal(unname(hy(half)[["lambda"]]), 1)))

  # both held, and the fit reports exactly what was written
  both <- statmod(y ~ enet(Z, lambda = 3, alpha = 0.25),
                  distributions7::gaussian1_distrib(), de)
  expect_equal(unname(hy(both)[["lambda"]]), 3)
  expect_equal(unname(hy(both)[["alpha"]]), 0.25)

  # a smooth answers the same way, through a criterion of the other kind
  ds <- data.frame(x = stats::runif(200))
  ds$y <- sin(6 * ds$x) + stats::rnorm(200, sd = 0.3)
  est <- statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), ds)
  hel <- statmod(y ~ s(x, k = 10, lambda = 5),
                 distributions7::gaussian1_distrib(), ds)
  expect_equal(unname(unlist(hel@hyper$mu[[1L]])), 5)
  expect_false(isTRUE(all.equal(unname(unlist(est@hyper$mu[[1L]])), 5)))
  # and holding it is what the summary says
  s <- summary(hel)
  k <- vapply(s@tables$mu, `[[`, character(1), "kind")
  tb <- s@tables$mu[[which(k == "smooth")]]$table
  expect_identical(tb$source[tb$name == "lambda"], "fixed")
})

test_that("a held hyperparameter travels through a subformula", {
  # the value is written in a sub-term of another term's formula, and the
  # entry carries it out to the fit
  set.seed(62)
  n2 <- 200
  dd2 <- data.frame(x = stats::runif(n2, 0, 10),
                    id = factor(rep(1:5, each = n2 / 5)))
  dd2$y <- 1 + 0.5 * dd2$x + stats::rnorm(n2, sd = 0.4)
  fit <- statmod(y ~ nl(~ a * exp(-r * x), a ~ 0 + ridge(~ id, sigma = 0.7),
                        start = list(a = 1, r = 0.2)),
                 distributions7::gaussian1_distrib(), dd2)
  th <- unlist(fit@hyper$mu, use.names = TRUE)
  expect_true(any(abs(th - 0.7) < 1e-12))
})
