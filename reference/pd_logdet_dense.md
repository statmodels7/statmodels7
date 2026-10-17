# The Dense Route of pd_logdet

Applies the three routes described at
[`pd_logdet()`](https://statmodels7.github.io/statmodels7/reference/pd_logdet.md)
to a dense matrix.

## Usage

``` r
pd_logdet_dense(M, scale = NULL)
```

## Arguments

- M:

  A dense symmetric matrix.

- scale:

  A reference magnitude.

## Value

A list with `logdet` and `ok`; `factor`, the Cholesky factor, when the
Cholesky route accepted the matrix; and `min_ev` and `max_ev` when the
eigendecomposition rejected it.

## Details

This function is split out so that
[`pd_factor()`](https://statmodels7.github.io/statmodels7/reference/pd_factor.md)
can reach it as the fallback of the sparse route without restating the
verdict: it is the one place that determines whether a matrix is
positive definite, and it holds the one set of thresholds.

## See also

[`pd_logdet()`](https://statmodels7.github.io/statmodels7/reference/pd_logdet.md)
