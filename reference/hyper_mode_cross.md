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

  Whether those rows include that tail. `FALSE`, the default, gives the
  matrix over the coefficients alone.

## Value

A list with `cross`, an `n` by `nrow(idx)` matrix, and `skipped`, how
many penalties contributed nothing to it.

## Details

This function computes the cross derivative of the penalty that the
movement of the penalized mode requires. Three functions read it here
instead of assembling it separately:
[`statmod_edf_correction()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf_correction.md),
which contracts the movement with the information to obtain the
additional effective degrees of freedom,
[`hyper_correction()`](https://statmodels7.github.io/statmodels7/reference/hyper_correction.md),
which adds the resulting matrix to a variance, and
[`marginal_blocks()`](https://statmodels7.github.io/statmodels7/reference/marginal_blocks.md),
which builds the movement of the integrated coefficients.

A shared hyperparameter is one column that stands for several penalties,
so each member writes into the column of its group and the contributions
accumulate. Where nothing is shared, each member has its own column.

A penalty over the own parameters of a structural term covers positions
among those parameters, which are not columns of a design. Where the
matrix spans the coefficients alone, there is no position to write to
and the penalty is skipped. The skipped penalties are counted in
`skipped`, because a correction assembled without them is incomplete.

With `joint` the matrix spans the vector in which the mode of such a
model moves, the coefficients followed by the free parameters of the
term, and the penalties over those parameters have rows.
[`unit_joint_positions()`](https://statmodels7.github.io/statmodels7/reference/unit_joint_positions.md)
returns the positions of the coordinates of each unit in that vector and
[`unit_joint_beta()`](https://statmodels7.github.io/statmodels7/reference/unit_joint_beta.md)
returns their values.
[`hyper_correction()`](https://statmodels7.github.io/statmodels7/reference/hyper_correction.md)
and
[`statmod_edf_correction()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf_correction.md)
pass `joint = TRUE` for a model that carries a filter, because
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
and
[`statmod_edf()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf.md)
use the joint penalized information.
[`marginal_blocks()`](https://statmodels7.github.io/statmodels7/reference/marginal_blocks.md)
uses the matrix over the coefficients alone.

## See also

[`statmod_edf_correction()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf_correction.md),
[`hyper_correction()`](https://statmodels7.github.io/statmodels7/reference/hyper_correction.md),
[`penalties7::penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.html)
