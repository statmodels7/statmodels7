# The extractors of a fitted model: what it reports, on which scale, and the
# one place the readable view is built so that three of them cannot disagree.

fit_simple <- function(seed = 101, n = 200) {
  set.seed(seed)
  d <- data.frame(x = runif(n, -2, 2),
                  g = factor(sample(letters[1:6], n, TRUE)))
  d$y <- 1 + 0.8 * d$x + rnorm(6, 0, 0.5)[as.integer(d$g)] +
    rnorm(n, 0, 0.4)
  statmod(y ~ x + random(~ 1 | g), gaussian1_distrib(), d)
}

fit_break <- function(seed = 71, n = 400, kind = "jseg") {
  set.seed(seed)
  d <- data.frame(x = sort(runif(n, 0, 10)), z = runif(n, -1, 1))
  d$y <- 0.4 * d$x + 2 * pmax(d$x - 6, 0) + 1.5 * (d$x > 6) + 0.5 * d$z +
    rnorm(n, 0, 0.3)
  fm <- switch(kind,
    seg = y ~ z + seg(x, psi = 4),
    jump = y ~ z + jump(x, psi = 4),
    jseg = y ~ z + jseg(x, psi = 4))
  statmod(fm, gaussian1_distrib(), d)
}

fit_gas <- function(seed = 72, n = 500, q = 1L) {
  set.seed(seed)
  y <- numeric(n)
  f <- 1
  for (i in seq_len(n)) {
    y[i] <- f + rnorm(1, 0, 1)
    f <- 0.4 + 0.15 * (y[i] - f) + 0.6 * f
  }
  statmod(y ~ 0 + gas(p = 1, q = q, time = t), gaussian1_distrib(),
          data.frame(y = y, t = seq_len(n)))
}


test_that("hyper reports every hyperparameter, on either scale", {
  fit <- fit_simple()
  h <- hyper(fit)
  expect_true(all(c("parameter", "term", "name", "estimate", "held",
                    "source") %in% names(h)))
  expect_identical(nrow(h), 1L)
  expect_identical(h$name, "sigma")
  expect_identical(h$source, "reml")
  expect_false(h$held)
  # the free scale is the one the outer search ran on, through the
  # hyperparameter's own link
  hl <- hyper(fit, scale = "link")
  expect_equal(hl$estimate, log(h$estimate), tolerance = 1e-12)
  # a model with no penalty answers with no rows, not with NULL
  set.seed(9)
  d <- data.frame(x = runif(50))
  d$y <- d$x + rnorm(50, 0, 0.3)
  h0 <- hyper(statmod(y ~ x, gaussian1_distrib(), d))
  expect_identical(nrow(h0), 0L)
  expect_true(all(c("parameter", "estimate") %in% names(h0)))
})

test_that("a held hyperparameter is reported as held", {
  set.seed(11)
  d <- data.frame(z = runif(200, 0, 1))
  d$y <- sin(3 * d$z) + rnorm(200, 0, 0.3)
  h <- hyper(statmod(y ~ s(z, k = 6, lambda = 10), gaussian1_distrib(), d))
  expect_true(h$held)
  expect_identical(h$source, "fixed")
  expect_equal(h$estimate, 10)
})


test_that("coef reports the quantities, and the coordinates when asked", {
  fit <- fit_break(kind = "jseg")
  a <- coef(fit)$mu
  b <- coef(fit, readable = FALSE)$mu
  # the position is a quantity of the model and the working coordinate it is
  # read off is not: -g/delta against g
  expect_true(any(endsWith(names(a), "psi1")))
  expect_false(any(endsWith(names(a), ".g1")))
  expect_true(any(endsWith(names(b), ".g1")))
  expect_false(any(endsWith(names(b), "psi1")))
  # everything the term does not reparametrize is the same number
  shared <- intersect(names(a), names(b))
  expect_gt(length(shared), 1L)
  expect_equal(unname(a[shared]), unname(b[shared]))
  expect_identical(length(a), length(b))
})

test_that("coef reports a structural term's parameters", {
  fit <- fit_gas()
  a <- coef(fit)$mu
  b <- coef(fit, readable = FALSE)$mu
  # a model whose whole predictor is a filter used to answer numeric(0)
  expect_identical(length(a), 3L)
  expect_true(all(endsWith(names(a), c("omega", "alpha1", "beta1"))))
  # the persistence rides a partial autocorrelation and is REPORTED as the
  # autoregressive coefficient; the loading rides a log
  expect_true(any(endsWith(names(b), "pacf1")))
  expect_equal(unname(a[endsWith(names(a), "alpha1")]),
               exp(unname(b[endsWith(names(b), "alpha1")])))
})

test_that("a developed parameter keeps its coefficients under either reading", {
  set.seed(31)
  m <- 6
  ni <- 40
  id <- factor(rep(seq_len(m), each = ni))
  x <- runif(m * ni, 0, 10)
  psi <- 5 + rnorm(m, 0, 0.5)
  d <- data.frame(y = 0.3 * x + 1.5 * pmax(x - psi[as.integer(id)], 0) +
                    rnorm(m * ni, 0, 0.4), x = x, id = id)
  fit <- statmod(y ~ seg(x, psi ~ random(~ 1 | id), psi = 5),
                 gaussian1_distrib(), d)
  a <- coef(fit)$mu
  b <- coef(fit, readable = FALSE)$mu
  # a development is a vector over covariates with no single value to
  # report, so the term declares nothing for it and nothing is taken away
  expect_identical(length(a), length(b))
  expect_equal(unname(a), unname(b))
})


test_that("vcov and confint report the quantities, and agree with coef", {
  fit <- fit_gas(q = 1L)
  V <- vcov(fit)
  ci <- confint(fit)
  a <- unlist(coef(fit), use.names = TRUE)
  # the three views are built from one map, so their names line up
  expect_identical(rownames(V), rownames(ci))
  expect_equal(ci$estimate, unname(a[match(
    sub("^[^:]+:", "", rownames(ci)), sub("^[^.]+\\.", "", names(a)))]),
    tolerance = 1e-10)
  # the delta method really moves the variance: the loading's is its own
  # scale's, not its coordinate's
  Vr <- vcov(fit, readable = FALSE)
  i <- grep("alpha1$", rownames(V))
  j <- grep("alpha1$", rownames(Vr))
  al <- ci$estimate[grep("alpha1$", rownames(ci))]
  expect_equal(sqrt(V[i, i]), al * sqrt(Vr[j, j]), tolerance = 1e-6)
})

test_that("an interval is built on the scale that keeps its quantity in its set", {
  fit <- fit_gas(q = 1L)
  ci <- confint(fit)
  r <- ci[grep("alpha1$", rownames(ci)), , drop = FALSE]
  # a loading is positive and rides a log, so its interval is asymmetric and
  # its lower end cannot be negative
  expect_gt(r$lower, 0)
  expect_gt(r$upper - r$estimate, r$estimate - r$lower)
})

test_that("vcov selects one equation's submatrix", {
  fit <- fit_simple()
  V <- vcov(fit)
  Vm <- vcov(fit, parameter = "mu")
  expect_lt(nrow(Vm), nrow(V))
  expect_true(all(startsWith(rownames(Vm), "mu:")))
  expect_equal(Vm, V[rownames(Vm), rownames(Vm)])
  expect_error(vcov(fit, parameter = "nope"), "matched nothing")
})

test_that("a block that is a working linearization is held out, not reported at zero", {
  sharp <- fit_break(kind = "jump")
  expect_warning(V <- vcov(sharp), "working linearization")
  nm <- rownames(V)
  own <- grep("^mu:jump\\.", nm)
  rest <- setdiff(seq_along(nm), own)
  # the term's own coefficients would be reported with a standard error of
  # EXACTLY zero -- measured against a bootstrap of 200 resamples, 0.000
  # against 0.063 and 0.540 -- so they are missing instead
  expect_true(all(is.na(diag(V)[own])))
  expect_true(all(is.finite(diag(V)[rest])))
  expect_gt(length(own), 0L)
  # and the continuous construction, whose block IS a Jacobian, keeps
  # everything
  cont <- fit_break(kind = "seg")
  Vc <- expect_no_warning(vcov(cont))
  expect_true(all(is.finite(diag(Vc))))
  expect_true(any(grepl("psi1$", rownames(Vc))))
})

test_that("a summary of a sharp break-point exists and says why it is short", {
  sharp <- fit_break(kind = "jump")
  s <- expect_no_warning(summary(sharp))
  expect_true(any(grepl("working linearization", s@notes)))
  out <- paste(utils::capture.output(print(s)), collapse = "\n")
  expect_match(out, "jump(x, psi = 4)", fixed = TRUE)
})


test_that("the quantile residual is standard normal under a correct model", {
  set.seed(101)
  n <- 2000
  d <- data.frame(x = runif(n, -1, 1))
  # a skewed family, where a Pearson residual would not be normal even here
  d$y <- rgamma(n, shape = 4, rate = 4 / exp(1 + 0.5 * d$x))
  r <- residuals(statmod(y ~ x, gamma1_distrib(), d), seed = 1)
  expect_identical(length(r), as.integer(n))
  expect_lt(abs(mean(r)), 0.1)
  expect_lt(abs(stats::sd(r) - 1), 0.1)
  expect_gt(suppressWarnings(stats::ks.test(r, "pnorm")$p.value), 0.01)
})

test_that("a discrete family is randomized, and reproducibly", {
  set.seed(102)
  n <- 1500
  d <- data.frame(x = runif(n, -1, 1))
  d$y <- rpois(n, exp(1 + 0.8 * d$x))
  fit <- statmod(y ~ x, poisson_distrib(), d)
  r <- residuals(fit, seed = 3)
  expect_gt(suppressWarnings(stats::ks.test(r, "pnorm")$p.value), 0.01)
  # exact again, at the price of being random
  expect_identical(residuals(fit, seed = 3), r)
  expect_false(identical(residuals(fit, seed = 4), r))
  # and the caller's stream is left where it was
  set.seed(555)
  u <- stats::runif(1)
  invisible(residuals(fit, seed = 3))
  set.seed(555)
  expect_identical(stats::runif(1), u)
})

test_that("the quantile residual sees a modelled scale where a response one cannot", {
  set.seed(103)
  n <- 2000
  d <- data.frame(x = runif(n, -1, 1))
  d$y <- rnorm(n, 1 + 0.5 * d$x, exp(-0.5 + 0.8 * d$x))
  fit <- statmod(y ~ x | sigma ~ x, gaussian1_distrib(), d)
  expect_lt(abs(stats::sd(residuals(fit, seed = 1)) - 1), 0.1)
  # the response residual is the numerator alone and carries the varying
  # scale with it
  expect_gt(abs(stats::sd(residuals(fit, type = "response")) - 1), 0.15)
  expect_lt(abs(stats::sd(residuals(fit, type = "pearson")) - 1), 0.1)
})

test_that("a score-driven term's predictor carries its own uncertainty", {
  set.seed(66)
  m <- 300
  y <- numeric(m)
  f <- 1
  for (i in seq_len(m)) {
    y[i] <- f + rnorm(1, 0, 1)
    f <- 0.4 + 0.15 * (y[i] - f) + 0.6 * f
  }
  d <- data.frame(y = y, t = seq_len(m), x = rnorm(m))
  fit <- statmod(y ~ x + gas(p = 1, q = 1, time = t) - 1,
                 gaussian1_distrib(), d)
  spec <- fit@spec
  design <- statmod_design(spec)
  sst <- statmod_structural_state(design)
  tn <- names(sst$zeta)[[1L]]
  free <- setdiff(names(sst$zeta[[tn]]), sst$held[[tn]])
  nb <- design$mu$npar

  # THE REFERENCE differentiates the predictor itself in every estimated
  # parameter, coefficients and the filter's own together, and shares no
  # arithmetic with the delta method: the level is a recursion, so a
  # coefficient of this equation reaches it through the scores as well as
  # through the static part, and leaving that out understates the answer
  eta_of <- function(v) {
    cf <- fit@coefficients
    cf$mu <- v[seq_len(nb)]
    zz <- sst$zeta[[tn]]
    for (k in seq_along(free)) zz[[free[k]]] <- v[nb + k]
    d2 <- statmod_design(spec)
    s2 <- statmod_structural_state(d2)
    s2$zeta[[tn]] <- zz
    s2$key <- NULL
    as.numeric(statmod_eta(spec, d2, cf)$eta$mu)
  }
  v0 <- c(fit@coefficients$mu, unlist(sst$zeta[[tn]])[free])
  J <- numDeriv::jacobian(eta_of, v0)
  key <- c(paste("mu", design$mu$coef_names, sep = ":"),
           structural_tail_names(spec, design))
  V <- as.matrix(vcov(fit, readable = FALSE)[key, key, drop = FALSE])
  ref <- sqrt(pmax(rowSums((J %*% V) * J), 0))
  got <- predict(fit, "link:mu", se = TRUE)$se
  expect_equal(got, ref, tolerance = 1e-6)

  # and the static row alone is NOT the answer, so the check above cannot
  # be satisfied by dropping the propagation
  X <- as.matrix(design$mu$X)
  Vb <- V[seq_len(nb), seq_len(nb), drop = FALSE]
  naive <- sqrt(pmax(rowSums((X %*% Vb) * X), 0))
  expect_gt(max(abs(naive - ref) / ref), 0.05)
})

test_that("a score-driven term is continued past the series", {
  set.seed(66)
  m <- 200
  y <- numeric(m)
  f <- 1
  for (i in seq_len(m)) {
    y[i] <- f + rnorm(1, 0, 1)
    f <- 0.4 + 0.15 * (y[i] - f) + 0.6 * f
  }
  d <- data.frame(y = y, t = seq_len(m), x = rnorm(m))
  fit <- statmod(y ~ x + gas(p = 1, q = 1, time = t) - 1,
                 gaussian1_distrib(), d)
  tn <- names(fit@spec@terms$mu)[[which(vapply(fit@spec@terms$mu,
    function(z) S7::S7_inherits(z, modelterms7::structural_term),
    logical(1)))]]
  nd <- data.frame(y = NA_real_, t = m + 1:5, x = c(1, -1, 0.5, 0, 2))
  got <- predict(fit, "mu", nd)

  # THE REFERENCE rebuilds the term on the whole series and re-runs the
  # FILTER, so it shares the recursion and not the continuation. A gaussian
  # mean's score is (y - eta)/sigma^2, and the continuation's whole content
  # is that past the data that score sits at its conditional mean of zero,
  # which is what a missing response produces here.
  ext <- rbind(d, nd)
  tmx <- modelterms7::term_build(fit@spec@terms$mu[[tn]], ext)
  ost <- statmod_structural_state(statmod_design(fit@spec))
  psi <- structural_psi(tmx, ost$zeta[[tn]])
  s2 <- statmod_respec(fit@spec, ext, need_response = FALSE)
  d2 <- statmod_design_at(s2, fit@coefficients, statmod_design(s2))
  es <- as.numeric(d2$mu$X %*% fit@coefficients$mu)
  sg <- statmod_eta(fit@spec, statmod_design(fit@spec),
                    fit@coefficients)$theta$sigma[[1L]]
  yv <- ext$y
  out <- modelterms7::term_filter(
    tmx, es, yv,
    function(e, i) if (is.na(yv[[i]])) 0 else (yv[[i]] - e) / sg^2,
    function(e, i) if (is.na(yv[[i]])) 0 else -1 / sg^2, psi)
  expect_equal(got, out$eta[m + 1:5])

  # a row inside the observed series is not a continuation of it
  expect_error(predict(fit, "mu", data.frame(y = NA_real_, t = 50, x = 0)),
               "inside the observed series")
  # and a forecast reports no standard error, the uncertainty of the future
  # scores being no delta method
  expect_error(predict(fit, "mu", nd, se = TRUE), "is not reported", fixed = TRUE)
})

test_that("a forecast decays towards the filter's stationary level", {
  set.seed(68)
  m <- 300
  y <- numeric(m)
  f <- 1
  for (i in seq_len(m)) {
    y[i] <- f + rnorm(1, 0, 1)
    f <- 0.4 + 0.15 * (y[i] - f) + 0.6 * f
  }
  d <- data.frame(y = y, t = seq_len(m))
  fit <- statmod(y ~ 0 + gas(p = 1, q = 1, time = t), gaussian1_distrib(), d)
  nd <- data.frame(y = NA_real_, t = m + 1:40)
  fc <- predict(fit, "mu", nd)
  # a property of the model rather than of the code: with the score at its
  # conditional mean the recursion is f <- omega + beta f, whose fixed point
  # is omega/(1 - beta), and the term reports both under their own names
  rd <- modelterms7::term_readable(fit@spec@terms$mu[[1L]],
                                   statmod_structural_state(
                                     statmod_design(fit@spec))$zeta[[1L]])
  om <- rd$value[[which(rd$name == "omega")]]
  bt <- rd$value[[which(rd$name == "beta1")]]
  expect_equal(fc[[40L]], om / (1 - bt), tolerance = 1e-6)
  expect_true(all(abs(diff(fc)) <= abs(diff(fc))[[1L]] + 1e-12))
})

test_that("the ordinary generics answer or say why they cannot", {
  set.seed(70)
  dd <- data.frame(x = runif(60))
  dd$y <- 1 + 2 * dd$x + rnorm(60, sd = 0.4)
  fit <- statmod(y ~ x | sigma ~ x, gaussian1_distrib(), dd)

  expect_identical(nobs(fit), 60L)
  expect_identical(formula(fit), fit@spec@formula)
  expect_identical(family(fit)@distrib_name, "gaussian1")
  expect_equal(weights(fit), rep(1, 60))
  expect_equal(df.residual(fit),
               60 - as.numeric(attr(logLik(fit), "df")))

  # the standard deviation is the RESPONSE's under the fitted law and not
  # whichever parameter is spelled sigma, so it comes from the family's own
  # moment rather than from a name
  expect_equal(sigma(fit), predict(fit, "std_dev"))

  X <- model.matrix(fit, "sigma")
  expect_identical(colnames(X), c("(Intercept)", "x"))
  expect_equal(as.numeric(X %*% coef(fit, readable = FALSE)$sigma),
               predict(fit, "link:sigma"))
  expect_error(model.matrix(fit, "nope"), "must name one of")

  s1 <- simulate(fit, nsim = 2, seed = 7)
  s2 <- simulate(fit, nsim = 2, seed = 7)
  expect_identical(s1, s2)
  expect_identical(dim(s1), c(60L, 2L))
  # the caller's stream is restored, so a seed here does not reach the next
  # draw the caller takes
  set.seed(3)
  a <- rnorm(1)
  set.seed(3)
  invisible(simulate(fit, seed = 99))
  expect_equal(rnorm(1), a)

  expect_error(terms(fit), "one set of terms per distribution parameter")
  expect_error(model.frame(fit), "does not keep the fitting data")
  expect_error(anova(fit), "no")
})
