# The Criterion a Fit Maximized, Rebuilt From the Fit

The specification and the design on which a fit's marginal criterion was
maximized: the coefficients
[`marginal_coords()`](https://statmodels7.github.io/statmodels7/reference/marginal_coords.md)
names are held at their fitted values, and the coordinates left out of
the determinant are recorded on the design.

## Usage

``` r
fit_marginal_context(spec, design, coef, method, pinned = NULL)
```

## Arguments

- spec:

  The fit's specification.

- design:

  Its design.

- coef:

  The fitted coefficients, a named list.

- method:

  The
  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md)
  the fit ran, or `NULL`.

- pinned:

  The coordinates left out of the determinant, as the fit records them
  in `methods$pinned`.

## Value

A list with `spec`, `design` and `gam`, the last as
[`marginal_coords()`](https://statmodels7.github.io/statmodels7/reference/marginal_coords.md)
returns it.

## Details

A consumer that reads the criterion's derivatives at a fitted object has
to read them on this criterion. On the fit's own specification, which
holds nothing, the same functions return the derivatives of a criterion
that integrates the estimated coefficients, which is another function.
