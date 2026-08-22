# Raise the Levenberg Damping

The next \\\lambda\\ to try after a step that failed or moved nothing.

## Usage

``` r
iwls_escalate(damp, pieces)
```

## Arguments

- damp:

  The current damping.

- pieces:

  What
  [`iwls_pieces`](https://statmodels7.github.io/statmodels7/reference/iwls_pieces.md)
  built, for the scale.

## Value

The next damping.

## Details

The first value is \\10^{-8}\\ of the curvature's own scale, which is
negligible against a well-curved coordinate and already large against
one whose information has fallen by eight orders; each further try
multiplies by a hundred, so eight tries span sixteen orders and reach a
damping that dominates the largest diagonal. Measured on the case this
exists for, the coordinate needed a \\\lambda\\ of about 100 against a
largest diagonal of 2328, which the sixth escalation passes.

## See also

[`iwls_scale`](https://statmodels7.github.io/statmodels7/reference/iwls_scale.md),
[`iwls_solve`](https://statmodels7.github.io/statmodels7/reference/iwls_solve.md)
