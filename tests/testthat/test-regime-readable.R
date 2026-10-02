sim_alternating <- function(n, seed) {
  # a chain whose first state is always followed by the second, as the
  # waiting times of Old Faithful are: p11 = 0, so its log-ratio runs to the
  # edge of the chart
  set.seed(seed)
  s <- integer(n)
  s[1] <- 2L
  for (t in 2:n) {
    s[t] <- if (s[t - 1L] == 1L) 2L else if (stats::runif(1) < 0.6) 1L else 2L
  }
  data.frame(t = seq_len(n), y = c(55, 80)[s] + stats::rnorm(n, sd = 6))
}

test_that("a regime fit reports the transition matrix, blank at the edge", {
  dd <- sim_alternating(200, 4)
  fit <- statmod(y ~ regime(k = 2, time = t),
                 distributions7::gaussian1_distrib(), dd)
  z <- fit@structural[[1]]$unconstrained
  skip_if(!(abs(z[["alr1.1"]]) > 8),
          "the first row did not reach the edge on this platform")
  ci <- confint(fit)
  pn <- paste0("mu:regime.p", c("1.1", "1.2", "2.1", "2.2"))
  expect_true(all(pn %in% rownames(ci)))
  P <- as.matrix(parameters7::param_value(parameters7::transition_matrix(2),
                                          z[c("alr1.1", "alr2.1")]))
  expect_equal(ci[pn, "estimate"], as.numeric(t(P)), tolerance = 1e-12)
  # the first row reads the coordinate at the edge, and carries no standard
  # error and no interval; the second row does
  expect_true(all(is.na(ci[pn[1:2], c("se", "lower", "upper")])))
  expect_true(all(is.finite(ci[pn[3:4], "se"])))
  expect_true(all(ci[pn[3:4], "lower"] > 0 & ci[pn[3:4], "upper"] < 1))
  # the second row's standard error is the delta method through the logistic
  # map of its own log-ratio
  V <- vcov(fit, readable = FALSE)
  a <- z[["alr2.1"]]
  p <- stats::plogis(a)
  expect_equal(ci[pn[3], "se"], p * (1 - p) * sqrt(V["mu:regime.alr2.1",
                                                     "mu:regime.alr2.1"]),
               tolerance = 1e-8)
  # the summary agrees, and says why the first row is blank
  out <- paste(utils::capture.output(print(summary(fit), notes = TRUE)),
               collapse = "\n")
  expect_match(out, "p1.1, p1.2 depend on a coordinate", fixed = TRUE)
})

test_that("statmod_latent() gives a regime's smoothed states", {
  dd <- sim_alternating(120, 4)
  fit <- statmod(y ~ regime(k = 2, time = t),
                 distributions7::gaussian1_distrib(), dd)
  lt <- statmod_latent(fit)
  expect_identical(dim(lt), c(120L, 2L))
  expect_equal(rowSums(lt), rep(1, 120), tolerance = 1e-12)
  # a forward-backward pass written out here, at the fitted parameters
  z <- fit@structural[[1]]$unconstrained
  m <- coef(fit)$mu[["(Intercept)"]] + c(0, exp(z[["gap2"]]))
  s <- exp(coef(fit)$sigma[["(Intercept)"]])
  ea <- exp(c(z[["alr1.1"]], z[["alr2.1"]]))
  P <- rbind(c(ea[1], 1) / (1 + ea[1]), c(ea[2], 1) / (1 + ea[2]))
  d0 <- c(P[2, 1], P[1, 2]) / (P[1, 2] + P[2, 1])
  n <- nrow(dd)
  f <- vapply(1:2, function(j) stats::dnorm(dd$y, m[j], s), numeric(n))
  al <- matrix(0, n, 2)
  a <- d0 * f[1, ]
  al[1, ] <- a / sum(a)
  for (t in 2:n) {
    a <- as.vector(al[t - 1, ] %*% P) * f[t, ]
    al[t, ] <- a / sum(a)
  }
  be <- matrix(1, n, 2)
  for (t in (n - 1):1) {
    b <- as.vector(P %*% (f[t + 1, ] * be[t + 1, ]))
    be[t, ] <- b / sum(b)
  }
  g <- al * be
  g <- g / rowSums(g)
  expect_equal(unname(as.matrix(lt)), g, tolerance = 1e-8)
})
