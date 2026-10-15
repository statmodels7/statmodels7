# Validate Prior Weights

Returns a vector of prior weights of the right length, defaulting to one
per observation.

## Usage

``` r
check_weights(weights, n)
```

## Arguments

- weights:

  The weights, or `NULL`.

- n:

  The number of observations.

## Value

A numeric vector of length `n`.

## Details

They are not normalized. Making them sum to one would turn \\\sum_i w_i
\ell_i\\ into a mean: every standard error would be multiplied by
\\\sqrt{n}\\ and the information criteria would stop being comparable
with those of an unweighted fit. The weights act as frequency weights.

The weights must be numeric, of length `n`, non-negative, finite and
complete, and not all zero; anything else is an error.
