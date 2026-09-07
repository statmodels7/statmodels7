# Which Coefficients a Restricted Fit Can Hold

The stacked positions of the coefficients the three restricted
statistics are defined for: every coordinate less the ones under a
kinked penalty, less the aliased ones, and none at all where the model
carries a structural term.

## Usage

``` r
testable_coords(fit, spec = fit@spec, design = statmod_design(spec))
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- spec, design:

  The specification and its design, where the caller has them already.

## Value

An integer vector of stacked positions, possibly empty, named by the
labels
[`coef_labels()`](https://statmodels7.github.io/statmodels7/reference/coef_labels.md)
gives them.

## Details

It is the enumeration
[`statmod_restrict()`](https://statmodels7.github.io/statmodels7/reference/statmod_restrict.md)
refuses one coordinate at a time, asked once for the whole vector
instead. A caller reporting a table needs it that way round, since each
restricted statistic is a refit and an interval by inversion is a few
dozen: a row that cannot be tested has to be recognized before it is
paid for rather than after.

A coordinate under a penalty that is twice differentiable IS included –
a smooth's linear column and its rotated coordinates, a ridge, a random
effect – and the restricted fit holds it exactly as it holds an
unpenalized one,
[`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md)
enforcing the hold in the solve. Only a kinked penalty is excluded, for
the reasons
[`statmod_restrict()`](https://statmodels7.github.io/statmodels7/reference/statmod_restrict.md)
gives.

## See also

[`statmod_restrict()`](https://statmodels7.github.io/statmodels7/reference/statmod_restrict.md),
which refuses the same coordinates one at a time, and
[`kinked_coords()`](https://statmodels7.github.io/statmodels7/reference/kinked_coords.md).
