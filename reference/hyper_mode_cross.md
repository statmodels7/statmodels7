# The Mixed Derivative of the Penalty in the Coefficients and the Hyperparameters

\\\partial^2\rho / \partial\beta \partial\theta\\, written into the
stacked coefficient vector with one column per estimated hyperparameter.

## Usage

``` r
hyper_mode_cross(spec, design, coef, hyper, idx, n, joint = FALSE)
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

  How many rows the matrix carries: the stacked coefficients, and with
  `joint` a structural term's free parameters after them.

- joint:

  Whether those rows include that tail. `FALSE`, the default, is the
  coefficient-only matrix, which is what
  [`statmod_edf_correction()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf_correction.md)
  contracts.

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
those parameters rather than columns of a design, so where the caller's
matrix spans the coefficients alone there is nowhere to write it and it
is skipped. It is counted rather than passed over in silence: a
correction assembled without it is incomplete, and a caller reporting to
a reader has to be able to say so.

With `joint` the matrix spans the vector the mode really moves in for
such a model, the coefficients followed by the term's own free
parameters, and those penalties have rows after all. What they lacked
was an address and not a derivative:
[`unit_joint_positions()`](https://statmodels7.github.io/statmodels7/reference/unit_joint_positions.md)
says where each unit's coordinates live in that vector and
[`unit_joint_beta()`](https://statmodels7.github.io/statmodels7/reference/unit_joint_beta.md)
reads their values, so nothing here is derived.
[`hyper_correction()`](https://statmodels7.github.io/statmodels7/reference/hyper_correction.md)
asks for it because
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
inverts the joint penalized information; the coefficient-space consumer
does not, and its answer is unchanged.

## See also

[`statmod_edf_correction()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf_correction.md),
[`hyper_correction()`](https://statmodels7.github.io/statmodels7/reference/hyper_correction.md),
[`penalties7::penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.html)
