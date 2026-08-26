# Assign the Observations to Folds

A fold number per observation, either the one the method carries or a
fresh permutation.

## Usage

``` r
cv_folds(n, k, folds)
```

## Arguments

- n:

  The number of observations.

- k:

  How many folds.

- folds:

  The method's own assignment, or `numeric(0)`.

## Value

An integer vector of length `n`.

## Details

The permutation is drawn from the caller's stream and put back, so a fit
is not silently reproducible only when the caller happens to have
seeded. Passing `folds` explicitly is how two criteria are made
comparable on the same partition.
