# The Curvature's Own Scale at One Point

The largest diagonal entry of the penalized information, read off the
pieces without assembling it.

## Usage

``` r
iwls_scale(pieces)
```

## Arguments

- pieces:

  What
  [`iwls_pieces`](https://statmodels7.github.io/statmodels7/reference/iwls_pieces.md)
  built.

## Value

A positive number, or 1 where the pieces say nothing.

## Details

It is what the Levenberg damping is measured against, so that
[`iwls_escalate`](https://statmodels7.github.io/statmodels7/reference/iwls_escalate.md)
carries no constant with units. On the augmented route the diagonal of
\\K = R'R + C'C\\ is the column sums of the squares, which keeps a
sparse design sparse; on the assembled route it is the diagonal itself.
