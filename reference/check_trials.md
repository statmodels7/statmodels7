# Reject a Number of Trials That Cannot Follow These Rows

Signals an error where the family carries one number of trials per
observation and the rows to which it is being applied are a different
number.

## Usage

``` r
check_trials(distrib, n)
```

## Arguments

- distrib:

  The family.

- n:

  The number of rows to which it is about to be applied.

## Value

`distrib`, unchanged.

## Details

A `size` given to `binomial_distrib()` is a vector fixed when the fit
was written, so at other rows it is recycled and the density returns a
value without any signal: a size of length 200 against 20 observations
returns 200 log-densities, so the family and not the data determines the
number of terms.

The function signals an error because it cannot return a correct number:
the fit does not record to which rows the trials belonged. The caller
can write the response as `cbind(successes, failures)` instead, whose
row sums are recomputed wherever the expression is evaluated, and the
error message gives that remedy.

## See also

[`split_binomial_matrix()`](https://statmodels7.github.io/statmodels7/reference/split_binomial_matrix.md),
[`statmod_respec()`](https://statmodels7.github.io/statmodels7/reference/statmod_respec.md)
