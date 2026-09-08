# A covariance class inside a structural term ------------------------------
#
# A filter's own parameters are estimated beside the coefficients, so a
# labelled random effect developing one of them carries a prior like any
# other. What differs is only where its coordinates are: among the term's
# own parameters, which contribute no design column, rather than in the
# stacked coefficient vector.

# a panel whose loading varies by group, which is the model the case exists
# for: alpha_i lognormal around 0.25, a persistence of 0.6 and a random
# intercept per group
panel <- function(seed = 51L, m = 8L, ni = 40L) {
  set.seed(seed)
  g <- factor(rep(seq_len(m), each = ni))
  al <- 0.25 * exp(stats::rnorm(m, 0, 0.5))
  bi <- stats::rnorm(m, 0, 0.6)
  y <- unlist(lapply(seq_len(m), function(i) {
    f <- 0
    out <- numeric(ni)
    for (t in seq_len(ni)) {
      mu <- bi[i] + f
      out[t] <- stats::rnorm(1, mu, 1)
      f <- al[i] * (out[t] - mu) + 0.6 * f
    }
    out
  }))
  data.frame(g = g, y = y)
}

test_that("a class of one structural member is the unlabelled fit exactly", {
  # THE CONTROL THAT SAYS THE ADDRESSING IS RIGHT. A label collecting one
  # member describes the same prior over the same coordinates as no label at
  # all, so every number has to come back the same -- and it does so through
  # the whole of the new route: a class is assembled, its penalty built from
  # the class's own default, and its value read at interleaved positions in
  # the term's parameter vector.
  dd <- panel()
  one <- statmod(y ~ gas(p = 1, q = 1, by = g,
                         alpha1 ~ 1 + random(~ 1 | u | g)),
                 distributions7::gaussian1_distrib(), dd,
                 outer_criterion = reml())
  bare <- statmod(y ~ gas(p = 1, q = 1, by = g,
                          alpha1 ~ 1 + random(~ 1 | g)),
                  distributions7::gaussian1_distrib(), dd,
                  outer_criterion = reml())
  expect_equal(as.numeric(stats::logLik(one)),
               as.numeric(stats::logLik(bare)), tolerance = 1e-12)
  expect_equal(one@coefficients, bare@coefficients, tolerance = 1e-12)
  expect_equal(sum(one@edf$edf), sum(bare@edf$edf), tolerance = 1e-12)
  expect_equal(unname(one@structural[[1L]]$parameter),
               unname(bare@structural[[1L]]$parameter), tolerance = 1e-12)
  # the hyperparameter is the same number under a different key: the class's,
  # since the prior belongs to the class and not to the member
  expect_equal(hyper(one)$estimate, hyper(bare)$estimate, tolerance = 1e-12)
  expect_identical(hyper(one)$term, "u | g")
})

test_that("the penalty a one-member class builds is the member's own", {
  # the mechanism behind the control above, checked directly so that a
  # failure there says which half moved
  dd <- panel()
  bare <- modelterms7::term_build(
    modelterms7::gas(p = 1, q = 1, by = g, alpha1 ~ 1 + random(~ 1 | g)), dd)
  lab <- modelterms7::term_build(
    modelterms7::gas(p = 1, q = 1, by = g, alpha1 ~ 1 + random(~ 1 | u | g)),
    dd)
  # a LABELLED sub-term declares no penalty of its own: the class carries it
  expect_length(modelterms7::term_penalties(lab), 0L)
  own <- modelterms7::term_penalties(bare)[[1L]]$penalty
  cls <- statmod_classes(list(mu = list(T = lab)))[[1L]]$penalty
  b <- stats::rnorm(nlevels(dd$g))
  th <- list(sigma = 0.7)
  expect_identical(penalties7::penalty_value(own, b, th),
                   penalties7::penalty_value(cls, b, th))
  expect_identical(penalties7::penalty_gradient(own, b, th),
                   penalties7::penalty_gradient(cls, b, th))
})

test_that("two developed parameters of one filter share a covariance", {
  # the case the work exists for: the level and the loading of one filter
  # vary by group and the two deviations are correlated
  dd <- panel()
  fit <- statmod(y ~ gas(p = 1, q = 1, by = g,
                         alpha1 ~ 1 + random(~ 1 | u | g),
                         omega ~ 1 + random(~ 1 | u | g)),
                 distributions7::gaussian1_distrib(), dd,
                 outer_criterion = reml())
  des <- statmod_design(fit@spec)
  u <- Filter(function(z) !is.null(z$pieces),
              statmod_penalized(fit@spec, des))[[1L]]
  expect_true(u$structural)
  expect_identical(u$class$dim, 2L)
  expect_length(u$cols, 2L * nlevels(dd$g))
  # INTERLEAVED GROUP BY GROUP, which is the order the prior reads its
  # argument in: for each group, that group's coordinate from each member
  expect_identical(u$cols[1:2],
                   c(u$pieces[[1L]]$cols[[1L]], u$pieces[[2L]]$cols[[1L]]))
  # the prior is over the two coordinates of one group, repeated over groups
  expect_identical(as.integer(u$penalty@block), 2L)

  # and the report names each coordinate for the parameter it develops, not
  # for the equation alone: both are effects inside one term of `mu`
  cd <- class_coords(u)
  expect_identical(cd$column,
                   c("omega:(Intercept)", "alpha1:(Intercept)"))
  s <- summary(fit)
  expect_length(s@classes, 1L)
  # the class's own count is the term's, taken ONCE: both members are
  # developments of the same term, and a sum over the pieces read its row
  # twice
  bits <- s@classes[[1L]]$bits
  ed <- as.numeric(sub("^edf ", "", grep("^edf ", bits, value = TRUE)))
  expect_lte(ed, sum(fit@edf$edf))
})

test_that("two classes in one filter accumulate into one vector", {
  # Both are addressed in the SAME term's parameters, and the accumulator is
  # built once per term, so the second must add to it rather than reset it.
  dd <- panel(m = 8L, ni = 20L)
  spec <- statmod_spec(y ~ 0 + gas(p = 1, q = 1, by = g,
                                   omega ~ 1 + random(~ 1 | u | g),
                                   alpha1 ~ 1 + random(~ 1 | v | g)),
                       distributions7::gaussian1_distrib(), dd)
  des <- statmod_design(spec)
  us <- Filter(function(z) !is.null(z$pieces), statmod_penalized(spec, des))
  expect_length(us, 2L)
  expect_true(all(vapply(us, function(z) isTRUE(z$structural), TRUE)))
  # two classes, disjoint coordinates, one term
  expect_length(intersect(us[[1L]]$cols, us[[2L]]$cols), 0L)
  expect_identical(us[[1L]]$term, us[[2L]]$term)

  # displaced by hand, because the START puts one class's coordinates at
  # exactly zero -- where a centered gaussian prior's gradient IS zero -- and
  # a check read there would pass with the second class doing all the work
  sst <- statmod_structural_state(des)
  v <- sst$zeta[[1L]]
  v[us[[1L]]$cols] <- seq(-0.4, 0.4, length.out = length(us[[1L]]$cols))
  v[us[[2L]]$cols] <- seq(0.3, -0.3, length.out = length(us[[2L]]$cols))
  sst$zeta[[1L]] <- v
  sst$key <- NULL
  sst$value <- NULL
  hy <- statmod_hyper_start(spec, des)
  g <- statmod_structural_penalty(spec, des, hy, "gradient")[[1L]]
  h <- statmod_structural_penalty(spec, des, hy, "hessian")[[1L]]
  expect_true(all(g[us[[1L]]$cols] != 0))
  expect_true(all(g[us[[2L]]$cols] != 0))
  expect_identical(sort(unname(which(rowSums(abs(h)) > 0))),
                   sort(c(us[[1L]]$cols, us[[2L]]$cols)))
  # at sigma = 1 a centered gaussian's gradient is the coordinate itself
  expect_equal(unname(g[us[[1L]]$cols]), v[us[[1L]]$cols],
               ignore_attr = TRUE, tolerance = 1e-12)
})

test_that("a model with no label is untouched by the class route", {
  # the negative control: the unit a class produces is new, but the two sites
  # a class reaches -- the penalty's value and the structural gradient -- are
  # on the path of every fit
  set.seed(7)
  n <- 300L
  g <- factor(rep(seq_len(15L), each = 20L))
  x <- stats::rnorm(n)
  b <- stats::rnorm(15L, 0, 0.5)
  d2 <- data.frame(y = 1 + 0.8 * x + b[as.integer(g)] + stats::rnorm(n, 0, 0.7),
                   x = x, g = g)
  fit <- statmod(y ~ x + random(~ 1 | g),
                 distributions7::gaussian1_distrib(), d2,
                 outer_criterion = reml())
  us <- statmod_penalized(fit@spec, statmod_design(fit@spec))
  expect_length(Filter(function(z) !is.null(z$pieces), us), 0L)
  expect_false(any(vapply(us, function(z) isTRUE(z$structural), TRUE)))
  expect_true(fit@converged)
})


# A class split between the coefficients and a filter's own parameters ------
#
# Its two halves are positions in two different vectors, and the vector that
# holds both is the one the inner step, the marginal criterion and the
# variance already assemble: the stacked coefficients, then the filter's free
# parameters. What such a class adds to those three is the cross block.

mixed_panel <- function(seed = 77L, m = 8L, ni = 30L, rho = 0) {
  set.seed(seed)
  S <- matrix(c(0.36, rho * 0.3, rho * 0.3, 0.25), 2, 2)
  z <- matrix(stats::rnorm(2 * m), m, 2) %*% chol(S)
  bi <- z[, 1L]
  al <- 0.25 * exp(z[, 2L])
  g <- factor(rep(seq_len(m), each = ni))
  y <- unlist(lapply(seq_len(m), function(i) {
    f <- 0
    out <- numeric(ni)
    for (t in seq_len(ni)) {
      mu <- bi[i] + f
      out[t] <- stats::rnorm(1, mu, 1)
      f <- al[i] * (out[t] - mu) + 0.6 * f
    }
    out
  }))
  data.frame(g = g, y = y)
}

mixed_formula <- y ~ random(~ 1 | u | g) +
  gas(p = 1, q = 1, by = g, alpha1 ~ 1 + random(~ 1 | u | g))

test_that("a mixed class is addressed in the joint vector, interleaved", {
  dd <- mixed_panel()
  spec <- statmod_spec(mixed_formula, distributions7::gaussian1_distrib(), dd)
  des <- statmod_design(spec)
  u <- Filter(function(z) isTRUE(z$mixed), statmod_penalized(spec, des))[[1L]]
  nb <- sum(vapply(des, function(d) d$npar, integer(1)))
  sst <- statmod_structural_state(des)
  key <- u$class$sterm
  free <- setdiff(names(sst$zeta[[key]]), sst$held[[key]])

  # neither vector on its own, so every reader written for one of the two
  # skips it and only what carries `joint` sees it
  expect_false(isTRUE(u$structural))
  expect_null(u$index)
  expect_null(u$cols)
  expect_length(u$joint, 2L * nlevels(dd$g))
  # half among the coefficients, half among the filter's FREE parameters
  expect_true(all(u$joint[c(TRUE, FALSE)] <= nb))
  expect_true(all(u$joint[c(FALSE, TRUE)] > nb))
  expect_true(all(u$joint <= nb + length(free)))
  # INTERLEAVED group by group, which is the order the prior reads
  expect_identical(u$joint[1:2],
                   c(u$pieces[[1L]]$joint[[1L]], u$pieces[[2L]]$joint[[1L]]))
  # a member inside the filter carries no design column, so a reader marking
  # design columns from the pieces cannot mark the wrong ones
  st <- Filter(function(z) isTRUE(z$structural), u$pieces)[[1L]]
  expect_null(st$cols)
})

test_that("the joint penalty matches a prior written out by hand", {
  dd <- mixed_panel()
  spec <- statmod_spec(mixed_formula, distributions7::gaussian1_distrib(), dd)
  des <- statmod_design(spec)
  u <- Filter(function(z) isTRUE(z$mixed), statmod_penalized(spec, des))[[1L]]
  nb <- sum(vapply(des, function(d) d$npar, integer(1)))
  sst <- statmod_structural_state(des)
  key <- u$class$sterm
  free <- setdiff(names(sst$zeta[[key]]), sst$held[[key]])
  njoint <- nb + length(free)

  # a point of the joint vector, written into both halves
  set.seed(5)
  w <- stats::rnorm(njoint, 0, 0.4)
  cf <- list(mu = w[seq_len(des$mu$npar)],
             sigma = w[des$mu$npar + seq_len(des$sigma$npar)])
  z <- sst$zeta[[key]]
  z[free] <- w[nb + seq_along(free)]
  sst$zeta[[key]] <- z
  sst$key <- NULL
  sst$value <- NULL

  hy <- statmod_hyper_start(spec, des)
  th <- hy[[u$param]][[u$key]]
  Sig <- parameters7::param_value(parameters7::dr_prod(2L), as.numeric(th))
  B <- matrix(w[u$joint], ncol = 2L, byrow = TRUE)
  Si <- solve(Sig)

  # THE NORMALIZING CONSTANT IS KEPT, which is what makes a large prior scale
  # expensive and is the only term standing between the criterion and a
  # runaway
  hand_v <- sum(0.5 * rowSums((B %*% Si) * B)) +
    nrow(B) * 0.5 * log(det(2 * pi * Sig))
  expect_equal(joint_penalty_at(spec, des, cf, hy, "value"), hand_v,
               tolerance = 1e-10)

  hand_g <- numeric(njoint)
  hand_g[u$joint] <- as.numeric(t(B %*% Si))
  expect_equal(joint_penalty_at(spec, des, cf, hy, "gradient", njoint), hand_g,
               tolerance = 1e-10)

  hand_h <- matrix(0, njoint, njoint)
  for (i in seq_len(nrow(B))) {
    ii <- u$joint[(i - 1L) * 2L + 1:2]
    hand_h[ii, ii] <- Si
  }
  got_h <- joint_penalty_at(spec, des, cf, hy, "hessian", njoint)
  expect_equal(got_h, hand_h, tolerance = 1e-10)

  # THE CROSS BLOCK is what did not exist before, and it is not zero: it is
  # the only place the correlation between a coefficient and a filter own
  # parameter enters at all
  jb <- u$joint[c(TRUE, FALSE)]
  jz <- u$joint[c(FALSE, TRUE)]
  expect_gt(max(abs(got_h[jb, jz])), 0.01)
  expect_equal(got_h[jb, jz], hand_h[jb, jz], tolerance = 1e-10)
})

test_that("at zero correlation the mixed model IS the independent one", {
  # The strongest control available, and an identity rather than a tolerance:
  # a block-diagonal prior is exactly the two separate priors, so the two
  # models are the same model and every piece of the criterion has to agree.
  dd <- mixed_panel()
  s1 <- 0.55
  s2 <- 0.40
  mt <- reml()
  crit <- function(form, set) {
    sp <- statmod_spec(form, distributions7::gaussian1_distrib(), dd)
    hy <- set(statmod_hyper_start(sp, statmod_design(sp)), sp)
    a <- fit_at_hyper(form, distributions7::gaussian1_distrib(), dd, hy)
    statmod_marginal(a$spec, a$design, a$coefficients, hy, mt,
                     basis = integrated_basis(a$spec, a$design, mt@kind))
  }
  mix <- crit(mixed_formula, function(hy, sp) {
    u <- Filter(function(z) isTRUE(z$mixed),
                statmod_penalized(sp, statmod_design(sp)))[[1L]]
    hy[[u$param]][[u$key]] <- c(sigma_log_sd1 = log(s1),
                                sigma_log_sd2 = log(s2), sigma_z2.1 = 0)
    hy
  })
  ind <- crit(y ~ random(~ 1 | g) +
                gas(p = 1, q = 1, by = g, alpha1 ~ 1 + random(~ 1 | g)),
              function(hy, sp) {
                for (u in statmod_penalized(sp, statmod_design(sp))) {
                  hy[[u$param]][[u$key]][["sigma"]] <-
                    if (isTRUE(u$structural)) s2 else s1
                }
                hy
              })
  expect_false(is.null(mix))
  expect_false(is.null(ind))
  expect_equal(mix$value, ind$value, tolerance = 1e-9)
  expect_equal(mix$loglik, ind$loglik, tolerance = 1e-9)
  expect_equal(mix$penalty, ind$penalty, tolerance = 1e-9)
  expect_equal(mix$logdet, ind$logdet, tolerance = 1e-9)
  expect_identical(mix$q, ind$q)
})

test_that("a class coordinate that is held is refused", {
  # The prior dimension is fixed when the class is assembled. A coordinate
  # held afterwards is not in the joint vector at all, so the block would have
  # one fewer than the prior describes. Not reachable from the formula
  # language today -- what a linear intercept holds is the DEVELOPMENT own
  # intercept, never one of the collected effects -- so the guard is exercised
  # by holding one here.
  dd <- mixed_panel()
  spec <- statmod_spec(mixed_formula, distributions7::gaussian1_distrib(), dd)
  des <- statmod_design(spec)
  u <- Filter(function(z) isTRUE(z$mixed), statmod_penalized(spec, des))[[1L]]
  sst <- statmod_structural_state(des)
  key <- u$class$sterm
  hit <- names(sst$zeta[[key]])[
    Filter(function(z) isTRUE(z$structural), u$pieces)[[1L]]$zcols[[1L]]]
  sst$held[[key]] <- c(sst$held[[key]], hit)
  err <- tryCatch(statmod_penalized(spec, des), error = conditionMessage)
  expect_type(err, "character")
  expect_match(err, hit, fixed = TRUE)
  expect_match(err, "which is held", fixed = TRUE)
})

test_that("the exact outer gradient of a mixed class is answered at order 1", {
  # It used to be refused at both orders, because the member loop of
  # statmod_structural_grad() addresses a unit by `index` -- which a mixed
  # class does not have -- so the member passed with an EMPTY scatter and the
  # gradient came back exactly zero, which reads as stationarity. That
  # function assembles on the joint vector already, so what such a class
  # needed was its positions in it. Order 2 stays refused, and not for a
  # reason of its own: statmod_marginal_hess() is written over the stacked
  # coefficients and has no joint twin, so every structural model is refused
  # there.
  dd <- mixed_panel()
  spec <- statmod_spec(mixed_formula, distributions7::gaussian1_distrib(), dd)
  des <- statmod_design(spec)
  idx <- outer_hyper_index(spec, statmod_blocks(spec, des))
  mt <- reml()
  expect_identical(nrow(idx), 3L)
  expect_true(mixed_penalized(spec, des))
  expect_true(outer_gradient_ok(spec, des, idx, mt, 1L))
  expect_false(outer_gradient_ok(spec, des, idx, mt, 2L))
  # the class's structural half has to answer term_third() like any other, and
  # it is named on the class rather than on the unit
  u <- Filter(function(z) isTRUE(z$mixed), statmod_penalized(spec, des))[[1L]]
  expect_false(is.null(find_term(spec, u$class$sterm)))

  expect_identical(class(outer_default_optimizer(TRUE, FALSE, TRUE))[[1L]],
                   class(optimizers7::lbfgs())[[1L]])
  # a mixed class whose gradient is unavailable for some further reason still
  # keeps lbfgs, and a model with no mixed class keeps the simplex it had
  expect_identical(class(outer_default_optimizer(FALSE, FALSE, TRUE))[[1L]],
                   class(optimizers7::lbfgs())[[1L]])
  expect_identical(class(outer_default_optimizer(FALSE, FALSE, FALSE))[[1L]],
                   class(optimizers7::nelder_mead())[[1L]])
})


test_that("that gradient converges on a difference of the criterion", {
  skip_on_cran()
  # THE REFERENCE SHARES NO ARITHMETIC WITH IT: a central difference of the
  # criterion with the mode refitted from a fresh design at every point, the
  # design being rebuilt because a structural term's own parameters live in
  # its state and reusing one reads the criterion off a non-stationary point.
  # O(h^2) is what separates a correct gradient from one missing a term, which
  # is flat in h, and from a badly located mode, which grows as 1/h.
  dd <- mixed_panel()
  spec0 <- statmod_spec(mixed_formula, distributions7::gaussian1_distrib(), dd)
  des0 <- statmod_design(spec0)
  idx <- outer_hyper_index(spec0, statmod_blocks(spec0, des0))
  mt <- reml()
  hy0 <- statmod_hyper_start(spec0, des0)

  at <- function(eta) {
    hy <- eta_to_hyper(eta, idx, hy0)
    a <- fit_at_hyper(mixed_formula, distributions7::gaussian1_distrib(), dd,
                      hy)
    b <- integrated_basis(a$spec, a$design, mt@kind)
    list(spec = a$spec, design = a$design, coef = a$coefficients, hy = hy,
         basis = b)
  }
  crit <- function(eta) {
    r <- at(eta)
    statmod_marginal(r$spec, r$design, r$coef, r$hy, mt, basis = r$basis)$value
  }
  r0 <- at(rep(0, nrow(idx)))
  g <- statmod_marginal_grad(r0$spec, r0$design, r0$coef, r0$hy, mt, idx,
                             r0$basis)
  # it is not the zero vector it used to be
  expect_true(all(abs(g) > 1e-3))

  rel <- vapply(c(3e-3, 1e-3), function(h) {
    cd <- vapply(seq_len(nrow(idx)), function(k) {
      ep <- rep(0, nrow(idx)); em <- ep
      ep[k] <- h; em[k] <- -h
      (crit(ep) - crit(em)) / (2 * h)
    }, 0)
    max(abs(g - cd)) / max(abs(cd))
  }, 0)
  expect_lt(rel[[2L]], 1e-5)
  # and it FALLS as the step does, which a missing term would not
  expect_lt(rel[[2L]], rel[[1L]] / 3)
})
