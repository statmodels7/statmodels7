# Fitting the coefficients at GIVEN hyperparameters.
#
# statmod() has no argument for that any more: which hyperparameters are held
# is said by the terms, and a caller writes it where the penalty is named. The
# checks below need the other thing -- the same model refitted at a sequence
# of hyperparameters, the formula being a parameter of the harness rather than
# something written out -- so they call what statmod() itself calls once the
# criteria are out of the way. Rebuilding the specification is part of that:
# starting from a previous fit's coefficients leaves a refreshable term's
# state where that fit left it, and the mode reached then is not the one the
# criterion is read at.

#
# `polish = TRUE` finishes the fit with Newton steps on the exact penalized
# Hessian until the score is at rounding. The inner fit stops on its
# objective-stall guard with the score still well above rounding (1.6e-4
# summed on a negative binomial with a modelled dispersion), and that
# residual is a smooth function of the hyperparameter, so a difference taken
# across refits turns it into a bias that does not shrink with the step: the
# exact outer Hessian read 4.1e-06 out at every step from 1e-2 to 1e-3, and
# 4.1e-07, 3.7e-08, 4.1e-09 with the mode polished. Where the objective's
# Hessian is not that of the penalized likelihood -- a block that moves with
# its coefficients, a structural term, a kinked penalty -- it is refused
# rather than done approximately.
fit_at_hyper <- function(formula, distrib, data, hy, inner = iwls(),
                         polish = FALSE) {
  spec <- statmod_spec(formula, distrib, data, NULL, NULL)
  design <- statmod_design(spec)
  blocks <- statmod_blocks(spec, design)
  cfg <- inner_settings(inner, distrib)
  obj <- statmod_objective(spec, hy, design, cfg$expected, cfg$approx)
  beta <- statmod_start(spec, design, obj, NULL)
  r <- statmod_alternate(spec, design, blocks, hy, inner, beta, cfg$expected,
                         cfg$approx, cfg$maxit, cfg$tol, verbosity(0))
  par <- r$par
  if (polish) {
    if (length(attr(design, "structural")) || length(attr(design, "refresh")) ||
        length(blocks$sparse)) {
      stop("fit_at_hyper(polish = TRUE) covers a fixed design with smooth ",
           "penalties only", call. = FALSE)
    }
    obs <- statmod_objective(spec, hy, design, FALSE, cfg$approx)
    for (k in seq_len(8L)) {
      g <- obs$gr(par)
      if (max(abs(g)) < 1e-11) break
      par <- par - as.numeric(solve(as.matrix(obs$he(par)), g))
    }
    if (max(abs(obs$gr(par))) > 1e-8) {
      stop("fit_at_hyper(polish = TRUE) did not reach the mode", call. = FALSE)
    }
  }
  list(coefficients = r$obj$split(par), spec = spec, design = design)
}
