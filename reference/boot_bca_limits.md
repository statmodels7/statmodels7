# The BCa Limits of One Quantity

Computes the two limits of the BCa interval of one coordinate from its
replicates and from the refits without one unit.

## Usage

``` r
boot_bca_limits(q, est, jk, alpha)
```

## Arguments

- q:

  The replicates of the quantity.

- est:

  Its estimate.

- jk:

  Its values in the refits without one unit.

- alpha:

  One minus the confidence level.

## Value

A numeric vector of length two, with the attribute `extreme`.

## Details

The bias correction is the normal quantile of the proportion of
replicates below the estimate, with the replicates equal to it counted
by halves. The acceleration is the skewness form of the jackknife values
(Efron 1987). A jackknife value that is not finite is left out. The
adjusted levels are undefined where \\1 - a (z_0 + z) \le 0\\ or where
the bias correction is infinite, and the limits are then `NA`. The
attribute `extreme` of the result is `TRUE` when an adjusted level
selects the first or the last replicate.

## References

Efron, B. (1987). Better bootstrap confidence intervals. *Journal of the
American Statistical Association*, 82(397), 171–185.
