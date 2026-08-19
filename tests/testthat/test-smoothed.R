# The smoothed break-point terms through the layer: a smoothed term is a
# Jacobian block, so it takes the Gauss-Newton embedding of seg() -- no
# working-fit phase -- and the developments of a break-point the sharp
# discontinuous constructions refuse become fittable. The plan's measured
# references (piano_segmented_random.txt, section 5): on seg and, in the
# gaussian, on jseg the smoothed route matches the native one at cor 0.99,
# and the probit correction takes the apparent scale towards the truth.

fitted_term <- function(fit, what) {
  nm <- grep(what, names(fit@spec@terms$mu), fixed = TRUE, value = TRUE)[1L]
  fit@spec@terms$mu[[nm]]
}

test_that("a smoothed seg lands where the native one lands", {
  set.seed(41)
  dd <- data.frame(x = sort(stats::runif(250, 0, 10)))
  dd$y <- 1 + 0.5 * dd$x + 2 * pmax(dd$x - 6, 0) + stats::rnorm(250, sd = 0.3)
  f0 <- statmod(y ~ seg(x), distributions7::gaussian1_distrib(), dd)
  f1 <- statmod(y ~ seg(x, smoothed = penalties7::smooth_probit(h = 0.1)),
                distributions7::gaussian1_distrib(), dd)
  p0 <- as.numeric(modelterms7::seg_psi(fitted_term(f0, "seg(")))
  p1 <- as.numeric(modelterms7::seg_psi(fitted_term(f1, "seg(")))
  expect_lt(abs(p0 - p1), 0.1)
  expect_true(f1@converged)
})

test_that("a smoothed jump is routed as a Jacobian, not through working
          fits", {
  set.seed(42)
  dd <- data.frame(x = sort(stats::runif(300, 0, 10)))
  dd$y <- 1 + 2 * (dd$x >= 6) + stats::rnorm(300, sd = 0.3)
  fit <- statmod(y ~ jump(x, smoothed = penalties7::smooth_probit(h = 0.15)),
                 distributions7::gaussian1_distrib(), dd)
  expect_equal(as.numeric(modelterms7::seg_psi(fitted_term(fit, "jump("))),
               6, tolerance = 0.15)
  # no working-fit phase ran: the trace holds no working block entries
  expect_false(any(grepl("working", fit@history$blocks$block)))
})

test_that("the smoothed jseg fits the random changepoint the sharp one
          refuses, and the probit correction reads through the summary", {
  set.seed(61)
  m <- 12
  ni <- 30
  id <- factor(rep(seq_len(m), each = ni))
  xr <- as.numeric(replicate(m, sort(stats::runif(ni, 0, 10))))
  psi_i <- stats::rnorm(m, 5, 0.5)
  mu <- 1 + 0.5 * xr + 2 * (xr >= psi_i[as.integer(id)]) -
    1.2 * pmax(xr - psi_i[as.integer(id)], 0)
  dr <- data.frame(y = mu + stats::rnorm(m * ni, sd = 0.4), x = xr, id = id)
  # the sharp construction rejects this model outright
  expect_error(
    statmod(y ~ jseg(x, psi ~ random(~ 1 | id)),
            distributions7::gaussian1_distrib(), dr),
    "development")
  fit <- statmod(y ~ jseg(x, psi ~ random(~ 1 | id),
                          smoothed = penalties7::smooth_probit()),
                 distributions7::gaussian1_distrib(), dr)
  psi_hat <- modelterms7::seg_psi(fitted_term(fit, "jseg("))
  per_id <- vapply(seq_len(m), function(i)
    mean(psi_hat[as.integer(id) == i, 1]), numeric(1))
  expect_gt(stats::cor(per_id, psi_i), 0.9)
  # the summary names the smoother and the width, and prints the apparent
  # scale beside the corrected one -- the probit's convolution identity
  txt <- paste(utils::capture.output(print(summary(fit))), collapse = "\n")
  expect_match(txt, "smoothed \\(probit")
  expect_match(txt, "corrected")
})
