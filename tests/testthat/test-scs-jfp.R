# A Jeffreys block beside scs(): the outer gradient and Hessian read the
# third derivative of the scs() block through the movement of the Jeffreys
# penalty, and agree with central differences of the criterion and of the
# gradient with the mode refitted.

test_that("the outer gradient and Hessian are exact beside scs()", {
  skip_if_not_installed("numDeriv")
  set.seed(6)
  n <- 200
  d <- data.frame(z = stats::rbinom(n, 1, 0.3), x = stats::runif(n))
  d$y <- stats::rpois(n, exp(0.3 + 1.2 * d$x^2 + 0.5 * d$z))
  P <- distributions7::poisson_distrib()
  fm <- y ~ modelterms7::scs(x, modelterms7::shape_increasing(),
                             basis7::bspline_smooth(k = 5)) +
    modelterms7::jfp(~ z)
  method <- ml(information = observed_info())
  fit0 <- statmod(fm, distrib = P, data = d, outer_criterion = NULL)
  spec <- fit0@spec
  design <- statmod_design(spec)
  idx <- outer_hyper_index(spec, statmod_blocks(spec, design))
  basis <- integrated_basis(spec, design, method@kind)
  at <- function(eta) {
    hy <- eta_to_hyper(eta, idx, fit0@hyper)
    obj <- statmod_objective(spec, hy, design, FALSE)
    par <- obj$stack(fit_at_hyper_spec(spec, design, hy, P))
    for (k in 1:8) {
      gg <- obj$gr(par)
      if (max(abs(gg)) < 1e-10) break
      J <- numDeriv::jacobian(obj$gr, par)
      par <- par - as.numeric(solve((J + t(J)) / 2, gg))
    }
    list(hy = hy, cf = obj$split(par))
  }
  fn <- function(eta) {
    a <- at(eta)
    statmod_marginal(spec, design, a$cf, a$hy, method, basis = basis)$value
  }
  gr <- function(eta) {
    a <- at(eta)
    statmod_marginal_grad(spec, design, a$cf, a$hy, method, idx, basis)
  }
  eta <- hyper_to_eta(fit0@hyper, idx) + 0.3
  a <- at(eta)
  g <- gr(eta)
  gn <- (fn(eta + 1e-4) - fn(eta - 1e-4)) / 2e-4
  expect_lt(abs(g - gn) / abs(gn), 1e-7)
  H <- statmod_marginal_hess(spec, design, a$cf, a$hy, method, idx, basis)
  Hn <- (gr(eta + 1e-4) - gr(eta - 1e-4)) / 2e-4
  expect_lt(abs(H - Hn) / abs(Hn), 1e-7)
})
