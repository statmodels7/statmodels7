# The Constraints That the Kinks Place on the Stacked Coefficients

Reads, through
[`modelterms7::term_kink_constraints()`](https://statmodels7.github.io/modelterms7/reference/term_kink_constraints.html),
the linear constraints that hold the objective on its kinks at the
coefficients given, for every term that recomputes its own block, and
places their rows in the columns of the stacked vector.

## Usage

``` r
kink_constraints(spec, design, coef, split, total)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  The coefficients, a named list with one numeric vector per parameter.

- split:

  The stacked positions of each parameter's coefficients, a named list
  as the objective's `split` returns it for
  [`seq_along()`](https://rdrr.io/r/base/seq.html) of the stacked
  vector.

- total:

  The length of the stacked vector.

## Value

A numeric matrix with `total` columns and one row per constraint, with
no rows where no term has a kink.

## Details

A break-point that sits on an observation is a kink of the objective.
Where its position is a single coefficient the row is that coefficient's
unit vector; where the position is developed over a design the row is
the development's row of that observation, and the coefficients can
still move along the hyperplane on which the break-point stays on the
observation.
[`iwls_fit()`](https://statmodels7.github.io/statmodels7/reference/iwls_fit.md)
solves on that subspace where a line search rejects every step length,
and
[`mode_solver()`](https://statmodels7.github.io/statmodels7/reference/mode_solver.md)
moves the mode in it.

## See also

[`iwls_fit()`](https://statmodels7.github.io/statmodels7/reference/iwls_fit.md),
[`mode_solver()`](https://statmodels7.github.io/statmodels7/reference/mode_solver.md),
[`refresh_units()`](https://statmodels7.github.io/statmodels7/reference/refresh_units.md).
