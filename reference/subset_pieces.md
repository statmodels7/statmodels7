# The Scoring Pieces of a Subset of the Coefficients

The penalized information over the coordinates `idx`, in the assembled
form
[`iwls_fit()`](https://statmodels7.github.io/statmodels7/reference/iwls_fit.md)
reads when a block is solved with the others held.

## Usage

``` r
subset_pieces(spec, design, coef, hyper, method, idx)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  A named list of coefficient vectors.

- hyper:

  The hyperparameters.

- method:

  An
  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
  method, its `hessian` read.

- idx:

  Integer positions of the subset in the stacked coefficients.

## Value

A list with `R = NULL`, `C = NULL` and `A`, the penalized information
over `idx`.

## Details

A subset has no square-root route of its own: that route needs the whole
design, and the solve on a subset reads the assembled matrix anyway. So
the information is formed over the subset's own columns through
[`statmod_information_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_information_at.md)'s
`index`, and the penalty's Hessian is added on the same rows and
columns. Forming the whole \\X'WX\\ to keep its submatrix was measured
at 54 per cent of a lasso path at \\n = 5000\\ and 200 columns, the
subset being the two intercepts.

Where the whole route went through the square roots it read \\R'R +
C'C\\, which equals \\X'WX + S\\ up to rounding wherever the factors
exist and was replaced by exactly this sum where they do not.

## See also

[`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md),
[`iwls_pieces()`](https://statmodels7.github.io/statmodels7/reference/iwls_pieces.md)
