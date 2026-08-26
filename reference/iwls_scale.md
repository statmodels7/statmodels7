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
  [`iwls_pieces()`](https://statmodels7.github.io/statmodels7/reference/iwls_pieces.md)
  built, carrying `R` and `C`, or `A`.

## Value

A single positive number: the largest diagonal entry. `1` where the
pieces carry nothing usable, which makes the damping schedule absolute
there instead of relative.

## Details

The Levenberg damping is measured against this, so
[`iwls_escalate()`](https://statmodels7.github.io/statmodels7/reference/iwls_escalate.md)
carries no constant with units of its own.

On the augmented route the diagonal of \\K = R'R + C'C\\ is the column
sums of the squares of \\R\\ and \\C\\, which reads a sparse design
without densifying it. On the assembled route it is the matrix's own
diagonal.

## See also

[`iwls_escalate()`](https://statmodels7.github.io/statmodels7/reference/iwls_escalate.md),
which divides by this,
[`iwls_pieces()`](https://statmodels7.github.io/statmodels7/reference/iwls_pieces.md)
for the input.
