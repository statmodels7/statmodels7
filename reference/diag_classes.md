# Classes of a Covariate with Equal Counts

Divides the observations by a covariate into `n_inter` classes of about
equal count (the quantile intervals of a numeric variable), or into the
levels of a factor.

## Usage

``` r
diag_classes(v, n_inter)
```

## Arguments

- v:

  The covariate.

- n_inter:

  The number of classes of a numeric variable.

## Value

A factor with one element per observation.
