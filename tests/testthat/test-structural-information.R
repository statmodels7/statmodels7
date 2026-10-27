# The expected information of a model with a structural term. A filter has
# the expected information conditional on the past,
#   sum_t J_t' I(eta_t) J_t,
# which the variance objects and the inner fit read where expected_info() is
# named; an outer criterion rejects it, and a term of the likelihood shape
# (regime(), a marginal break-point term) has none. The references are
# written in the test: J_t by numDeriv over the predictors, the per-observation
# scores by numDeriv over the per-observation log-likelihood.

sinf_gas_data <- function(n = 120, omega = 0.3, kappa1 = 0.4, b1 = 0.7,
                          seed = 31) {
  set.seed(seed)
  f <- numeric(n)
  y <- numeric(n)
  s <- 0
  f0 <- omega / (1 - b1)
  for (t in seq_len(n)) {
    f[t] <- omega + kappa1 * (if (t > 1) s else 0) +
      b1 * (if (t > 1) f[t - 1] else f0)
    y[t] <- f[t] + stats::rnorm(1, sd = 0.8)
    s <- y[t] - f[t]
  }
  data.frame(t = seq_len(n), y = y, x = stats::rnorm(n))
}

# the predictors, or the per-observation log-likelihood, as a function of
# the joint vector (coefficients, then the filter's free parameters)
sinf_at <- function(fit) {
  spec <- fit@spec
  design <- statmod_design(spec)
  obj <- statmod_objective(spec, fit@hyper, design, expected = FALSE)
  jp <- statmod_joint_pieces(spec, design, obj, fit@hyper)
  u0 <- c(unlist(fit@coefficients, use.names = FALSE), jp$zeta())
  ev_at <- function(u) {
    jp$setz(u[jp$ix])
    statmod_eta(spec, design, obj$split(u[seq_len(jp$nb)]))
  }
  list(spec = spec, design = design, jp = jp, u0 = u0, ev_at = ev_at,
       eta = function(u, p) ev_at(u)$eta[[p]],
       loglik = function(u) {
         ev <- ev_at(u)
         distributions7::distrib_pdf(spec@distrib, spec@response, ev$theta,
                                     log = TRUE)
       })
}

test_that("the conditional information is sum_t J_t' I(eta_t) J_t", {
  dd <- sinf_gas_data()
  fit <- statmod(y ~ x + gas(p = 1, q = 1, time = t),
                 distrib = distributions7::gaussian1_distrib(), data = dd)
  a <- sinf_at(fit)
  Jm <- numDeriv::jacobian(function(u) a$eta(u, "mu"), a$u0)
  Js <- numDeriv::jacobian(function(u) a$eta(u, "sigma"), a$u0)
  sig <- a$ev_at(a$u0)$theta$sigma
  sig <- rep_len(sig, nrow(Jm))
  # the gaussian on (identity mu, log sigma): I = diag(1/sigma^2, 2)
  ref <- crossprod(Jm / sig) + 2 * crossprod(Js)
  got <- statmod_full_information(a$spec, fit@coefficients, a$design,
                                  expected = TRUE)
  expect_equal(unname(got), unname(ref), tolerance = 1e-6)
  expect_gt(min(eigen(got, symmetric = TRUE, only.values = TRUE)$values), 0)
  # and it is not the observed information, whose second sum it drops
  obs <- statmod_full_information(a$spec, fit@coefficients, a$design)
  expect_gt(max(abs(obs - got)) / max(abs(obs)), 1e-3)
})

test_that("expect_opg() gives the outer product of the per-observation scores", {
  dd <- sinf_gas_data()
  fit <- statmod(y ~ x + gas(p = 1, q = 1, time = t),
                 distrib = distributions7::gaussian1_distrib(), data = dd)
  a <- sinf_at(fit)
  G <- numDeriv::jacobian(a$loglik, a$u0)
  a$ev_at(a$u0)
  got <- statmod_full_information(a$spec, fit@coefficients, a$design,
                                  expected = TRUE,
                                  emethod = distributions7::expect_opg())
  expect_equal(unname(got), unname(crossprod(G)), tolerance = 1e-6)
})

test_that("vcov() inverts the conditional information where expected_info() is named", {
  dd <- sinf_gas_data()
  fit <- statmod(y ~ x + gas(p = 1, q = 1, time = t),
                 distrib = distributions7::gaussian1_distrib(), data = dd)
  design <- statmod_design(fit@spec)
  K <- statmod_full_information(fit@spec, fit@coefficients, design,
                                expected = TRUE)
  ve <- vcov(fit, vcov_type = bayesian_vcov(
    information = distributions7::expected_info()), readable = FALSE)
  vo <- vcov(fit, readable = FALSE)
  # vcov() lists the joint vector by equation; the reference is in the
  # order of the information (coefficients, then the filter's parameters)
  rn <- c(rownames(coef_labels(fit@spec, design)),
          unname(structural_tail_names(fit@spec, design)))
  expect_equal(unname(ve[rn, rn]), unname(solve(K)), tolerance = 1e-8)
  expect_gt(max(abs(sqrt(diag(ve)) / sqrt(diag(vo)) - 1)), 1e-3)
  # the frequentist convention on the same matrix is the same inverse
  vf <- vcov(fit, vcov_type = frequentist_vcov(
    information = distributions7::expected_info()), readable = FALSE)
  expect_equal(vf, ve, tolerance = 1e-8)
})

test_that("summary() reads the structural table on the variance it names", {
  dd <- sinf_gas_data()
  fit <- statmod(y ~ x + gas(p = 1, q = 1, time = t),
                 distrib = distributions7::gaussian1_distrib(), data = dd)
  vt <- bayesian_vcov(information = distributions7::expected_info())
  tab <- statmod_structural_table(fit, vcov_type = vt)
  V <- vcov(fit, vcov_type = vt, readable = FALSE)
  design <- statmod_design(fit@spec)
  st <- statmod_structural_state(design)
  tm_name <- names(fit@structural)
  tm <- fit@spec@terms$mu[[tm_name]]
  nm <- modelterms7::term_params(tm)
  free <- setdiff(nm, st$held[[tm_name]])
  rn <- unname(structural_tail_names(fit@spec, design))
  rd <- modelterms7::term_readable(tm, st$zeta[[tm_name]])
  J <- rd$jacobian[, match(free, nm), drop = FALSE]
  ref <- sqrt(diag(J %*% V[rn, rn, drop = FALSE] %*% t(J)))
  got <- tab$se[match(rd$name, tab$name)]
  # a level that the intercept carries is held and reported without one
  ok <- is.finite(ref) & !tab$held[match(rd$name, tab$name)]
  expect_true(any(ok))
  expect_equal(got[ok], unname(ref[ok]), tolerance = 1e-6)
})

test_that("inner_information = expected_info() fits the joint system by scoring", {
  dd <- sinf_gas_data()
  f_auto <- statmod(y ~ x + gas(p = 1, q = 1, time = t),
                    distrib = distributions7::gaussian1_distrib(), data = dd)
  f_exp <- statmod(y ~ x + gas(p = 1, q = 1, time = t),
                   distrib = distributions7::gaussian1_distrib(), data = dd,
                   inner_information = distributions7::expected_info())
  expect_true(f_exp@converged)
  expect_equal(as.numeric(logLik(f_exp)), as.numeric(logLik(f_auto)),
               tolerance = 1e-8)
  # auto_info() settles the gaussian on the expected information for iwls(),
  # and leaves the joint system on the observed one
  expect_true(S7::S7_inherits(joint_inner_info(f_auto@methods$smooth),
                              distributions7::ObservedInfo))
  expect_true(S7::S7_inherits(joint_inner_info(f_exp@methods$smooth),
                              distributions7::ExpectedInfo))
  it <- function(f) {
    b <- f@history$blocks
    sum(b$iterations[b$block == "joint"])
  }
  expect_false(identical(it(f_exp), it(f_auto)))
  # an optimizer settled from auto_info() records it, a named one does not
  m_auto <- inner_bind(optimizers7::newton(), distributions7::gaussian1_distrib(),
                       auto_info())
  m_exp <- inner_bind(optimizers7::newton(), distributions7::gaussian1_distrib(),
                      distributions7::expected_info())
  expect_true(S7::S7_inherits(joint_inner_info(m_auto),
                              distributions7::ObservedInfo))
  expect_true(S7::S7_inherits(joint_inner_info(m_exp),
                              distributions7::ExpectedInfo))
})

test_that("an outer criterion beside a filter rejects the expected information", {
  dd <- sinf_gas_data(n = 100)
  dd$g <- factor(rep(1:5, length.out = nrow(dd)))
  expect_error(
    statmod(y ~ x + random(~ 1 | g) + gas(p = 1, q = 1, time = t),
            distrib = distributions7::gaussian1_distrib(), data = dd,
            outer_criterion = reml(information = distributions7::expected_info())),
    "not available for\n  an outer criterion")
  # auto_info() is settled as the observed information beside a structural
  # term, even for a family whose rule is the expected information
  fit <- statmod(y ~ gas(p = 1, q = 1, time = t),
                 distrib = distributions7::student_t1_distrib(), data = dd)
  design <- statmod_design(fit@spec)
  m <- outer_resolve(reml(information = auto_info()),
                     distributions7::student_t1_distrib(), design)
  expect_true(S7::S7_inherits(m@hessian, distributions7::ObservedInfo))
  m0 <- outer_resolve(reml(information = auto_info()),
                      distributions7::student_t1_distrib())
  expect_true(S7::S7_inherits(m0@hessian, distributions7::ExpectedInfo))
})

test_that("a term of the likelihood shape has no expected information", {
  set.seed(5)
  st <- rep(c(1, 2), each = 30)[sample.int(60)]
  dd <- data.frame(y = stats::rnorm(60, mean = c(0, 3)[st]))
  fit <- statmod(y ~ regime(k = 2), distrib = distributions7::gaussian1_distrib(),
                 data = dd)
  expect_error(vcov(fit, vcov_type = bayesian_vcov(
    information = distributions7::expected_info())),
    "mixture over latent components")
  expect_error(statmod(y ~ regime(k = 2),
                       distrib = distributions7::gaussian1_distrib(), data = dd,
                       inner_information = distributions7::expected_info()),
               "'inner_information' is distributions7::expected_info\\(\\)")
  # the default resolves to the observed information there
  expect_identical(vcov(fit), vcov(fit, vcov_type = bayesian_vcov(
    information = distributions7::observed_info())))
  design <- statmod_design(fit@spec)
  expect_true(S7::S7_inherits(
    vcov_information(bayesian_vcov(), fit@spec@distrib, design),
    distributions7::ObservedInfo))
})
