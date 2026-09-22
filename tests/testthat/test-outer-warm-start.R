## The outer search warm-starts every trial from the INCUMBENT, the best
## usable point so far. Its accept/reject happens inside optimizers7 and is
## not visible to evaluate(), so writing the warm start at every usable point
## let a rejected trial move it; where the inner problem is multimodal the
## chain left the incumbent's basin and could not come back. Measured on the
## reference battery's gas-panel-alpha: the first evaluation read -1367.360
## at the starting hyperparameter, every later one -1370.656 from a rejected
## trial's filter state, all three routes backtracked to exhaustion and the
## fit reported the start at the WORSE of the two values.

test_that("the reported criterion is not worse than the best the search saw", {
  skip_on_cran()
  set.seed(3)
  m <- 12; ni <- 40
  pan <- data.frame(g = factor(rep(seq_len(m), each = ni)), x = runif(m * ni))
  form <- y ~ gas(p = 1, q = 1, by = g, alpha1 ~ 1 + random(~1 | g))
  G <- distributions7::gaussian1_distrib()
  set.seed(13)
  d <- rstatmod(form, G, pan)$data
  fit <- statmod(form, G, d, inner_optimizer = iwls(maxit = 10000))
  h <- fit@history$outer
  expect_gt(nrow(h), 1L)
  ## the invariant the defect broke: the point reported is the one the search
  ## stands on, so its criterion is the best the search evaluated
  expect_gte(fit@criterion, max(h$criterion) - 1e-6)
})
