# The Diagonal of the Information a Step Uses

Returns the diagonal of the information
[`iwls_score()`](https://statmodels7.github.io/statmodels7/reference/iwls_score.md)
normalizes by, one entry per stacked coefficient.

## Usage

``` r
iwls_info_diag(pieces)
```

## Arguments

- pieces:

  The pieces, as
  [`iwls_pieces()`](https://statmodels7.github.io/statmodels7/reference/iwls_pieces.md)
  builds them, carrying `R` or `A`.

## Value

An unnamed numeric vector, one entry per stacked coefficient.
[`iwls_pieces()`](https://statmodels7.github.io/statmodels7/reference/iwls_pieces.md)
always builds one of the two matrices, so there is no third branch.

## Details

On the augmented route it is the column sums of the squares of the
square-root design, whose cross-product is the unpenalized information,
so the penalty does not enter.

On the assembled route it is the diagonal of the matrix itself, which
folds the penalty in. That is a different normalizer, and it is accepted
because the assembled route is already the fallback: it is reached only
where no square root exists, and the reading it feeds arbitrates a
verdict rather than driving the loop.

## See also

[`iwls_score()`](https://statmodels7.github.io/statmodels7/reference/iwls_score.md),
the only caller.
