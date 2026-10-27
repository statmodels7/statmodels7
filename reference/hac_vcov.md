# The HAC Variance of the Coefficients

Computes \\V_b B V_b\\ for the coordinates in `cols`, with \\V_b\\ the
inverse penalized information that
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
has built and \\B\\ the HAC meat of the scores.

## Usage

``` r
hac_vcov(fit, design, Vb, cols, sw)
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

  The positions, among the unknowns of the joint information, of the
  coordinates that `Vb` covers.

- sw:

  The
  [HacSandwich](https://statmodels7.github.io/statmodels7/reference/HacSandwich.md)
  object.

## Value

A symmetric matrix of the dimension of `Vb`.

## See also

[`hac_sandwich()`](https://statmodels7.github.io/statmodels7/reference/hac_sandwich.md),
[`joint_scores()`](https://statmodels7.github.io/statmodels7/reference/joint_scores.md)
