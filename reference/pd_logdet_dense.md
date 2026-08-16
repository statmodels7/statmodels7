# The Dense Route of pd_logdet

The three routes described at
[`pd_logdet`](https://statmodels7.github.io/statmodels7/reference/pd_logdet.md),
on a dense matrix.

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

A list with `logdet`, `ok` and, on a refusal reached through the
eigendecomposition, `min_ev` and `max_ev`.

## Details

Split out so that
[`pd_factor`](https://statmodels7.github.io/statmodels7/reference/pd_factor.md)
can reach it as the fallback of the sparse route without restating the
verdict: there is one place that decides whether a matrix is positive
definite, and one set of thresholds.

## See also

[`pd_logdet`](https://statmodels7.github.io/statmodels7/reference/pd_logdet.md)
