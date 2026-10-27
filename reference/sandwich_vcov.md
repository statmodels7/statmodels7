# The Sandwich Variance of the Coefficients

Computes \\c\\V_b B V_b\\ for the coordinates in `cols`, with \\V_b\\
the inverse penalized information that
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
has built, \\B\\ the sum over observations (or over clusters) of the
outer product of the scores, and \\c\\ the correction of the type.

## Usage

``` r
sandwich_vcov(fit, design, Vb, cols, sw)
```

## Arguments

- fit:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- design:

  The design of the fit.

- Vb:

  The inverse penalized information over the coordinates `cols`.

- cols:

  The positions, among the coefficients stacked by equation, of the
  coordinates that `Vb` covers.

- sw:

  The
  [SandwichMethod](https://statmodels7.github.io/statmodels7/reference/SandwichMethod.md)
  object.

## Value

A symmetric matrix of the dimension of `Vb`.

## See also

[`statmod_score_obs()`](https://statmodels7.github.io/statmodels7/reference/statmod_score_obs.md)
for the scores,
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md),
which calls it.
