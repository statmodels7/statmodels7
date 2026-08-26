# Print a Fitted Model

Prints a compact report of a fitted model: the call, the distribution,
one line per equation naming its terms and the effective degrees of
freedom each spends, the conditional log-likelihood, the elapsed time
and whether every loop stopped on its own rule.

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
for the verdict behind the convergence line.
