# Refuse a Number of Trials That Cannot Follow These Rows

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

A `size` handed to `binomial_distrib()` is a vector fixed when the fit
was written, so at other rows it is neither right nor obviously wrong:
it is recycled, and the answer comes back with no complaint. Measured
directly on the density, a size of length 200 against 20 observations
returns **200** log-densities summing to -1125.63 and signals nothing –
the response is recycled against the size, so the number of terms is
decided by the family rather than by the data.

The alternative to erroring is not a better number, it is a guess:
nothing on the fit says which rows those trials belonged to. What the
caller can do instead is write the response as
`cbind(successes, failures)`, whose row sums are recomputed wherever the
expression is evaluated, and the message says so.

## See also

[`split_binomial_matrix()`](https://statmodels7.github.io/statmodels7/reference/split_binomial_matrix.md),
[`statmod_respec()`](https://statmodels7.github.io/statmodels7/reference/statmod_respec.md)
