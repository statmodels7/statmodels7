# Does a Structural Term Carry a Penalty of Its Own?

Reports whether any penalty of the model covers the own parameters of a
structural term, as opposed to a block of design columns. The deviations
of a score-driven filter over a panel are the case.

## Usage

``` r
structural_penalized(spec, design)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design, whose structural state holds the terms.

## Value

A single logical. `FALSE` for a model with no structural term, and for
one whose structural term declares no penalty.

## Details

The answer decides which matrix the Laplace determinant is taken over.
`TRUE` sends the criterion to
[`statmod_marginal_full()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_full.md),
which spans the coefficients **and** the term's free parameters; `FALSE`
leaves it over the coefficients alone.

## See also

[`statmod_marginal_full()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_full.md)
for the matrix it selects,
[`statmod_structural_penalty()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_penalty.md)
for the penalties it looks for.
