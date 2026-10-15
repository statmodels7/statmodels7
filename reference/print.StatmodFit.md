# Print a Fitted Model

Prints the call, the distribution and the number of rows of a fitted
model. For each distribution parameter it prints the formula and one
line per term, giving the number of coefficients (or of structural
parameters) and, where it differs from that number, the effective
degrees of freedom of the term. Where an outer criterion ran, it prints
the value of the criterion with the number of hyperparameter
evaluations. It then prints the elapsed time and, where the fit made
several passes over blocks, the number of passes. It prints no
log-likelihood and no convergence verdict.

[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
is the long form, with coefficients, standard errors and the
hyperparameters.

## Usage

``` r
# S3 method for class 'StatmodFit'
print(x, ...)
```

## Arguments

- x:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- ...:

  Unused.

## Value

`x`, invisibly, as a print method should.

## See also

[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
for the full report,
[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md)
for the reading of the point that the fit reached.
