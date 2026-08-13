# The subformulas of a term's own parameters, fitted end to end: a
# developed gas level over a panel, nl and seg developments with their
# sub-terms' hyperparameters estimated, and the level's subspace question.

gas_panel <- function(seed = 5, m = 4, Tn = 50, sd = 0.5,
                      om = c(0.1, 0.9, 0.5, 0.3)) {
  set.seed(seed)
  a1 <- 0.3
  b1 <- 0.6
  y <- numeric(m * Tn)
  for (l in seq_len(m)) {
    f <- numeric(Tn)
    s <- numeric(Tn)
    f0 <- om[l] / (1 - b1)
    for (t in seq_len(Tn)) {
      f[t] <- om[l] + a1 * (if (t > 1) s[t - 1] else 0) +
        b1 * (if (t > 1) f[t - 1] else f0)
      y[(l - 1) * Tn + t] <- f[t] + rnorm(1, sd = sd)
      s[t] <- y[(l - 1) * Tn + t] - f[t]
    }
  }
  data.frame(y = y, id = factor(rep(sprintf("s%d", seq_len(m)), each = Tn)),
             t = rep(seq_len(Tn), m))
}

test_that("a developed gas level is fitted end to end", {
  dd <- gas_panel()
  fs <- statmod(y ~ gas(p = 1, q = 1, omega ~ random(~1 | id),
                        by = id, time = t) - 1,
                distributions7::gaussian1_distrib(), dd)
  expect_true(fs@converged)
  # the random intercept's hyperparameter is estimated by the marginal
  # criterion, and the departures answer the simulated levels: group s2
  # carries the largest one, s1 the smallest
  h <- unlist(fs@hyper$mu)
  expect_length(h, 1L)
  zs <- fs@structural[[1L]]$unconstrained
  dev <- zs[grepl(".random", names(zs), fixed = TRUE)]
  expect_identical(unname(which.max(dev)), 2L)
  expect_identical(unname(which.min(dev)), 1L)
  # the loading and the persistence are recovered on the parameter scale;
  # the simulation drives the level with the RAW residual while the model's
  # score is (y - f)/sigma^2, so the loading the model estimates is the
  # simulated one times sigma^2 = 0.25
  expect_equal(unname(exp(zs[["alpha1"]])), 0.3 * 0.25, tolerance = 0.35)
  expect_equal(
    unname(linkfunctions7::linkinv(linkfunctions7::rhobit_link(),
                                   zs[["pacf1"]])), 0.6, tolerance = 0.2)
})

test_that("vcov and summary answer for a penalized structural term", {
  # the labels used to grow past the design when a penalty sat on a
  # structural term's own parameters -- its cols index the term's
  # parameter vector, not the equation's columns -- and vcov() died on
  # duplicate row names
  dd <- gas_panel()
  fit <- statmod(y ~ gas(p = 1, q = 1, omega ~ random(~1 | id),
                         by = id, time = t),
                 distributions7::gaussian1_distrib(), dd)
  v <- vcov(fit)
  expect_true(all(is.finite(diag(v))))
  expect_identical(rownames(v)[1L], "mu:(Intercept)")
  out <- capture.output(print(summary(fit)))
  expect_true(any(grepl("estimated", out)))
})

test_that("the intercept holds a developed level's constant coordinate", {
  dd <- gas_panel()
  fit <- statmod(y ~ gas(p = 1, q = 1, omega ~ random(~1 | id),
                         by = id, time = t),
                 distributions7::gaussian1_distrib(), dd)
  st <- fit@structural[[1L]]
  expect_identical(st$held, "omega.(Intercept)")
  expect_identical(unname(st$unconstrained[["omega.(Intercept)"]]), 0)
  # the intercept carries the stationary level the held coordinate cannot
  b1 <- 0.6
  expect_equal(unname(coef(fit)$mu[[1L]]), mean(c(0.1, 0.9, 0.5, 0.3)) /
                 (1 - b1), tolerance = 0.35)
})

test_that("an unpenalized shared level column is flagged", {
  set.seed(7)
  n <- 80
  dd <- data.frame(y = rnorm(n), z = as.numeric(scale(runif(n))),
                   t = seq_len(n))
  spec <- statmod_spec(y ~ z + gas(p = 1, q = 1, omega ~ z, time = t),
                       distributions7::gaussian1_distrib(), dd)
  expect_warning(statmod_design(spec), "already spans")
  # a penalized development of the same column is identified by its
  # penalty and passes in silence
  spec2 <- statmod_spec(y ~ z + gas(p = 1, q = 1, omega ~ ridge(~z),
                                    time = t),
                        distributions7::gaussian1_distrib(), dd)
  expect_silent(statmod_design(spec2))
})

test_that("nl with a penalized submodel estimates its hyperparameter", {
  set.seed(3)
  n <- 200
  dd <- data.frame(x = runif(n, 0, 3),
                   g = factor(rep(sprintf("g%d", 1:4), length.out = n)))
  amp <- c(g1 = 1.8, g2 = 2.4, g3 = 2.0, g4 = 1.8)
  dd$y <- amp[dd$g] * exp(-1.3 * dd$x) + rnorm(n, sd = 0.1)
  fit <- statmod(y ~ nl(~ a * exp(-r * x), a ~ ridge(~g),
                        start = list(a = 2, r = 1)),
                 distributions7::gaussian1_distrib(), dd)
  expect_true(fit@converged)
  # one hyperparameter, keyed term::parameter::subterm, estimated
  h <- unlist(fit@hyper$mu)
  expect_length(h, 1L)
  cf <- coef(fit)$mu
  expect_equal(unname(cf[["nl.r"]]), 1.3, tolerance = 0.1)
  # the departures answer the amplitudes: g2 is the large one
  dep <- cf[grepl("ridge", names(cf))]
  expect_identical(unname(which.max(dep)), which(grepl("gg2", names(dep))))
})

test_that("seg with per-group break-points is fitted by the layer", {
  set.seed(4)
  n <- 200
  dd <- data.frame(x = runif(n, 0, 10),
                   id = factor(rep(c("a", "b"), length.out = n)))
  truth <- c(a = 4, b = 6)
  dd$y <- 0.5 * dd$x + 2.5 * pmax(dd$x - truth[dd$id], 0) +
    rnorm(n, sd = 0.3)
  fit <- statmod(y ~ seg(x, psi ~ id), distributions7::gaussian1_distrib(),
                 dd)
  expect_true(fit@converged)
  cf <- coef(fit)$mu
  psi_a <- unname(cf[["seg.psi1.(Intercept)"]])
  psi_b <- psi_a + unname(cf[["seg.psi1.idb"]])
  expect_equal(psi_a, 4, tolerance = 0.3)
  expect_equal(psi_b, 6, tolerance = 0.3)
})
