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

fit_at_hyper <- function(formula, distrib, data, hy, inner = iwls()) {
  spec <- statmod_spec(formula, distrib, data, NULL, NULL)
  design <- statmod_design(spec)
  blocks <- statmod_blocks(spec, design)
  cfg <- inner_settings(inner)
  obj <- statmod_objective(spec, hy, design, cfg$expected, cfg$approx)
  beta <- statmod_start(spec, design, obj, NULL)
  r <- statmod_alternate(spec, design, blocks, hy, inner, beta, cfg$expected,
                         cfg$approx, cfg$maxit, cfg$tol, verbosity(0))
  list(coefficients = r$obj$split(r$par), spec = spec, design = design)
}
