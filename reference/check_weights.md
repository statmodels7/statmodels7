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
\ell_i\\ into a mean, which is the averaging trap the objective's own
scale already has to avoid: every standard error would shrink by
\\\sqrt{n}\\ and the information criteria would stop being comparable.
