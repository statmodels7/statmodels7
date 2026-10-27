# Reject a Resampling Scheme That the BCa Interval Cannot Use

Signals an error for a scheme whose entry has no influence of the
observations to accelerate by.

## Usage

``` r
boot_bca_check(object, resample)
```

## Arguments

- object:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- resample:

  A
  [resample_method](https://statmodels7.github.io/statmodels7/reference/resample_method.md)
  object, or `NULL`, which passes.

## Value

`NULL`, invisibly.

## Details

The acceleration is computed from the refits of the fit without one
observation or one cluster. The parametric scheme draws new responses
from the model, so it does not resample the data. The cluster scheme
beside a
[`modelterms7::random()`](https://statmodels7.github.io/modelterms7/reference/random.html)
term relabels the groups drawn more than once, and deleting a cluster
there has no counterpart in the stored coordinates.
