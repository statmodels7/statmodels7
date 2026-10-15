# Assign the Observations to Folds

A fold number per observation, either the assignment the method carries
or a permutation drawn with a fixed seed.

## Usage

``` r
cv_folds(n, k, folds)
```

## Arguments

- n:

  The number of observations.

- k:

  How many folds, at most `n`.

- folds:

  The method's own assignment, or `numeric(0)`.

## Value

An integer vector of length `n`.

## Details

The permutation is drawn after `set.seed(20260811L)`, and the caller's
random-number state is restored afterwards, so the partition is the same
whatever the caller's seed. Passing `folds` explicitly is how two
criteria are made comparable on a partition of the caller's own. The
values of `folds` are relabeled to the integers `1, 2, ...` in the order
of their factor levels, and a length other than `n` signals an error.
