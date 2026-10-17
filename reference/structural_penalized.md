# Does a Structural Term Carry a Penalty of Its Own?

Reports whether any penalty of the model covers the own parameters of a
structural term, as opposed to a block of design columns, as with the
deviations of a score-driven filter over a panel. A covariance class
that spans both the coefficients and those parameters counts as well.

## Usage

``` r
structural_penalized(spec, design)

mixed_penalized(spec, design)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design, whose structural state holds the terms.

## Value

A single logical. It is `TRUE` where a penalty covers a structural
term's own parameters, or a covariance class spans both the coefficients
and those parameters, and `FALSE` for a model with no structural term
and for one whose structural term has no such penalty.
`mixed_penalized()` returns `TRUE` for the second case only.

## Details

The answer determines the matrix over which the Laplace determinant is
taken. `TRUE` sends the criterion to
[`statmod_marginal_full()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_full.md),
which spans the coefficients and the term's free parameters; `FALSE`
leaves it over the coefficients only.

## See also

[`statmod_marginal_full()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_full.md)
for the matrix that it selects,
[`statmod_structural_penalty()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_penalty.md)
for the penalties that it seeks.
