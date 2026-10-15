# Reject a Number of Trials That Cannot Follow These Rows

Signals an error where the family carries one number of trials per
observation and the rows it is being applied to are a different number.

## Usage

``` r
check_trials(distrib, n)
```

## Arguments

- distrib:

  The family.

- n:

  The number of rows it is about to be applied to.

## Value

`distrib`, unchanged.

## Details

A `size` given to `binomial_distrib()` is a vector fixed when the fit
was written, so at other rows it is recycled and the density returns a
value without any signal: a size of length 200 against 20 observations
returns 200 log-densities, the number of terms being decided by the
family and not by the data.

The function signals an error because no correct number can be returned:
the fit does not record which rows the trials belonged to. The caller
can write the response as `cbind(successes, failures)` instead, whose
row sums are recomputed wherever the expression is evaluated, and the
error message gives that remedy.

## See also

[`split_binomial_matrix()`](https://statmodels7.github.io/statmodels7/reference/split_binomial_matrix.md),
[`statmod_respec()`](https://statmodels7.github.io/statmodels7/reference/statmod_respec.md)
