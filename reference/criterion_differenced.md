# The Outer Criterion Differenced by Refitting

The gradient and the curvature of the outer criterion in its free
coordinates, read by central differences of the criterion itself: the
model is refitted at each displaced point, warm-started at the fitted
coefficients, and the criterion is read at the refitted mode.

## Usage

``` r
criterion_differenced(fit, spec, design, blocks, method, idx, cf, hy, h = 0.01)
```

## Arguments

- fit:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- spec, design, blocks:

  The specification, its design and its blocks.

- method:

  The outer criterion.

- idx:

  The hyperparameter index,
  [`outer_hyper_index()`](https://statmodels7.github.io/statmodels7/reference/outer_hyper_index.md).

- cf, hy:

  The fitted coefficients and hyperparameters.

- h:

  The step on the free scale.

## Value

A list with the gradient `g` and the curvature `A` (minus the Hessian),
or `NULL` where a refit did not reach a mode.

## Details

With one coordinate the readings are \\g = (c\_+ - c\_-)/(2h)\\ and \\A
= -(c\_+ - 2c_0 + c\_-)/h^2\\; with two the four axial points give the
diagonal and the four diagonal points the cross term, \\(c\_{++} -
c\_{+-} - c\_{-+} + c\_{--})/(4h^2)\\. The centre \\c_0\\ is refitted
too, so every reading comes from the same inner iteration. The step is
\\h = 10^{-2}\\ on the free scale, where a hyperparameter's logarithm
moves by one per cent; the truncation error is of order \\h^2\\ and the
effect of an inner tolerance \\\epsilon\\ on the curvature of order
\\\epsilon/h^2\\.
