# What Each Distribution Parameter Reached

A single string that gives, for each parameter of the distribution, the
range of its fitted values.
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
adds it as a note where the fit did not converge.

## Usage

``` r
fitted_ranges(x)
```

## Arguments

- x:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

## Value

A single string, empty when the parameters cannot be read. It is a note
of
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
and not a line of
[`print.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/print.StatmodFit.md),
because it qualifies the fit instead of describing it.

## Details

The note is written for a lasso at a fixed hyperparameter, with a free
scale and a design that the model can interpolate. Fitting the
coefficients shrinks the residuals, which shrinks the scale and raises
the working weights, so that the penalty weighs relatively less and more
coefficients enter. A scale near zero in the note points to this case.

The function reports where the parameters ended up and does not name a
cause, because naming one would require a threshold for what counts as
running away.

## See also

[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
