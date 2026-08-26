# The Rows of the Joint Predictor Derivative

\\V\_{a,t} = \partial\eta\_{a,t}/\partial u\\ over the coefficients of
every equation followed by a structural term's own parameters: the
equation's design placed in its own columns, and for the equation
carrying a filter the forward Jacobian of the recursion, which is dense.

## Usage

``` r
joint_design_rows(spec, design, coef)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  The coefficients.

## Value

A list with `V` (one matrix per distribution parameter), `ap` (which
equation carries the filter), `keep` (the columns that survive a held
level), `ev`, `f` and the sizes.

## Details

It exists so that the contraction
[`u_vector()`](https://statmodels7.github.io/statmodels7/reference/u_vector.md)
performs is written once. The formula there does not change where a
filter is present; only its operand does, \\X\\ becoming \\\[X_p \mid
D\]\\. Everything downstream, the third derivative against the diagonal
of \\M\\ and the movement of the mode, then reads the joint vector with
no special case.

The rows are returned at the full width, a held level included, and the
caller drops it exactly as
[`statmod_full_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_full_information.md)
does.

## See also

[`u_vector()`](https://statmodels7.github.io/statmodels7/reference/u_vector.md),
[`statmod_full_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_full_information.md)
