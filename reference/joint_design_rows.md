# The Rows of the Joint Predictor Derivative

Computes \\V\_{a,t} = \partial\eta\_{a,t}/\partial u\\ over the
coefficients of every equation followed by a structural term's own
parameters. The function returns each equation's static design placed in
its own columns, with zeros in the columns of the structural term's
parameters. The callers replace the matrix of the equation that carries
the filter by the forward Jacobian of the recursion, which is dense.

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

`NULL` where the model has no filter. Otherwise it is a list with `V`
(one matrix per distribution parameter), `ap` (which equation carries
the filter), `nb`, `np` and `mf` (the numbers of coefficients, of the
term's parameters and the joint width), `keep` (the columns that survive
a held level), `zn` (the term's parameter names), `ev`, `f`, `params`
and `n`.

## Details

It exists so that the contraction that
[`u_vector()`](https://statmodels7.github.io/statmodels7/reference/u_vector.md)
carries out is written once. The formula there does not change where a
filter is present; only its operand does: \\X\\ becomes \\\[X_p \mid
D\]\\. Everything downstream, the third derivative against the diagonal
of \\M\\ and the movement of the mode, then reads the joint vector with
no special case.

The rows are returned at the full width, including a held level, and the
caller drops that level exactly as
[`statmod_full_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_full_information.md)
does.

## See also

[`u_vector()`](https://statmodels7.github.io/statmodels7/reference/u_vector.md),
[`statmod_full_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_full_information.md)
