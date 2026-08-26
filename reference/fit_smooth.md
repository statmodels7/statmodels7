# Fit the Smooth Block

Runs `inner_optimizer` on the jointly fitted coefficients, holding every
kinked block at its current values. This is the smooth half of one pass
of
[`statmod_alternate()`](https://statmodels7.github.io/statmodels7/reference/statmod_alternate.md).

## Usage

``` r
fit_smooth(obj, beta, idx, spec, design, hyper, method, vb)
```

## Arguments

- obj:

  The objective, as
  [`statmod_objective()`](https://statmodels7.github.io/statmodels7/reference/statmod_objective.md)
  returns it.

- beta:

  The current stacked coefficients, all of them.

- idx:

  The smooth block's positions within `beta`, an integer vector.

- spec:

  The specification.

- design:

  The design.

- hyper:

  The hyperparameters, held fixed.

- method:

  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
  or an optimizers7 optimizer.

- vb:

  The resolved verbosity.

## Value

A list of five:

- `par`:

  the full stacked coefficient vector, with the block's positions
  replaced.

- `value`:

  the objective there.

- `converged`:

  a single logical, the method's own verdict.

- `iterations`:

  how many the method took.

- `history`:

  a data frame, one row per iteration.

## Details

The objective handed to the optimizer is the full one restricted to
`idx`: the coefficients outside the block enter it as constants, so the
value it reports is the model's objective, never a block's own.

With
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
the step is solved through `fit_smooth()`'s own scoring loop, which
reads the exact information. With an optimizers7 optimizer the
objective's `fn`, `gr` and `he` are all supplied, so `newton()` gets the
exact second derivative instead of differencing the gradient: measured,
that is 2.5x on one smooth and 5.5x on three smooths with a random
effect, at identical answers.

## See also

[`statmod_alternate()`](https://statmodels7.github.io/statmodels7/reference/statmod_alternate.md),
the caller,
[`sparse_fit()`](https://statmodels7.github.io/statmodels7/reference/sparse_fit.md)
for the kinked half of the same pass,
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
for the default method.
