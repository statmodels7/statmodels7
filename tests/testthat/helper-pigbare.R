# A discrete family that APPROXIMATES its expected information: the mass, the
# score, the Hessian and the draws of pig1_distrib(), and no expected method of
# its own, so distributions7's base class answers with the outer product of the
# scores. Every shipped family computes its expected information exactly since
# distributions7 0.65.0, and the refusal of a criterion on that approximation
# is tested on this one.
pig_bare_distrib <- local({
  distrib_pdf <- distributions7::distrib_pdf
  distrib_gradient <- distributions7::distrib_gradient
  distrib_hessian <- distributions7::distrib_hessian
  distrib_rng <- distributions7::distrib_rng
  PigBare <- S7::new_class("PigBareSm", parent = distributions7::discrete_distrib,
                           package = NULL)
  S7::method(distrib_pdf, PigBare) <- function(distrib, y, theta, log = FALSE, ...) {
    distributions7::distrib_pdf(distributions7::pig1_distrib(), y, theta, log = log)
  }
  S7::method(distrib_gradient, PigBare) <- function(distrib, y, theta,
                                                    scale = c("parameter", "link"), ...) {
    distributions7::distrib_gradient(distributions7::pig1_distrib(), y, theta)
  }
  S7::method(distrib_hessian, PigBare) <- function(distrib, y, theta,
                                                   scale = c("parameter", "link"), ...) {
    distributions7::distrib_hessian(distributions7::pig1_distrib(), y, theta)
  }
  S7::method(distrib_rng, PigBare) <- function(distrib, n, theta, ...) {
    distributions7::distrib_rng(distributions7::pig1_distrib(), n, theta)
  }
  p <- distributions7::pig1_distrib()
  function() {
    PigBare(distrib_name = "pig bare", dimension = "univariate",
            bounds = c(0, Inf), params = c("mu", "sigma"), n_params = 2,
            params_bounds = p@params_bounds, link_params = p@link_params,
            params_smooth = c(mu = TRUE, sigma = TRUE),
            params_interpretation = p@params_interpretation)
  }
})
