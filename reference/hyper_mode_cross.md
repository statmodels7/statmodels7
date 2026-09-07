# The Mixed Derivative of the Penalty in the Coefficients and the Hyperparameters

\\\partial^2\rho / \partial\beta \partial\theta\\, written into the
stacked coefficient vector with one column per estimated hyperparameter.

## Usage

``` r
hyper_mode_cross(spec, design, coef, hyper, idx, n)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  The coefficients.

- hyper:

  The hyperparameters.

- idx:

  The hyperparameter index, from
  [`outer_hyper_index()`](https://statmodels7.github.io/statmodels7/reference/outer_hyper_index.md).

- n:

  How many stacked coefficients the design carries.

## Value

A list with `cross`, an `n` by `nrow(idx)` matrix, and `skipped`, how
many penalties contributed nothing to it.

## Details

This is the one ingredient of the penalized mode's movement that nothing
else computes, and both consumers of that movement read it here rather
than assembling it each:
[`statmod_edf_correction()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf_correction.md),
which contracts it against the information to price what estimating a
hyperparameter cost, and
[`hyper_correction()`](https://statmodels7.github.io/statmodels7/reference/hyper_correction.md),
which keeps the matrix and adds it to a variance.

A shared hyperparameter is ONE column standing for several penalties, so
each member writes into the group's column and they accumulate. Where
nothing is shared each member is its own row and the lookup is what was
here before the groups existed.

A penalty over a STRUCTURAL term's own parameters covers positions among
those parameters rather than columns of a design, so there is nothing to
write and it is skipped. It is counted rather than passed over in
silence: a correction assembled without it is incomplete, and a caller
reporting to a reader has to be able to say so.

## See also

[`statmod_edf_correction()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf_correction.md),
[`hyper_correction()`](https://statmodels7.github.io/statmodels7/reference/hyper_correction.md),
[`penalties7::penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.html)
