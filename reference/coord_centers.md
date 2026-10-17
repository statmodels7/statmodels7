# Whether a Penalized Block Is Solved with Its Equation's Intercept Profiled

Returns `TRUE` where the equation carries an unpenalized intercept that
is free to move, in which case
[`coord_fit()`](https://statmodels7.github.io/statmodels7/reference/coord_fit.md)
centers the block's columns with the working weights and solves it with
the intercept profiled out.

## Usage

``` r
coord_centers(spec, design, p, obj, beta)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- p:

  The equation, a parameter name.

- obj:

  The objective, as
  [`statmod_objective()`](https://statmodels7.github.io/statmodels7/reference/statmod_objective.md)
  returns it.

- beta:

  The stacked coefficients.

## Value

`TRUE` or `FALSE`.

## Details

The block is fitted while the other columns of its equation are held,
and the intercept is one of them. Where the columns are not centered,
each change of a coefficient moves the mean of the fit, which only the
intercept can take back, and the intercept is updated in another block.
The alternation between the two then converges at a rate set by how
close each column is to the constant. Centering the columns removes this
dependence, and the alternation then converges to the same point that a
fit with centered predictors reaches.

Solving with the intercept profiled out is the Frisch-Waugh-Lovell
reading of an unpenalized constant: the coefficients of the block are
those of the centered columns against the centered working response, and
they do not depend on the value that the intercept holds.
[`coord_fit()`](https://statmodels7.github.io/statmodels7/reference/coord_fit.md)
then sets the intercept to the value that goes with them, the weighted
mean of the working response net of the block, so its step is a joint
step in the block and the intercept. A held intercept cannot take that
value, so an intercept named in `held_coef` leaves the block uncentered,
and so does an equation with no intercept at all.

## See also

[`coord_fit()`](https://statmodels7.github.io/statmodels7/reference/coord_fit.md),
[`coord_colsq()`](https://statmodels7.github.io/statmodels7/reference/coord_colsq.md),
[`parametric_intercept()`](https://statmodels7.github.io/statmodels7/reference/parametric_intercept.md)
