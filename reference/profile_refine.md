# A Profile Over Two Labels on a Finer Grid

Interpolates the criterion of a profile over two labels, given on a
complete rectangular grid, to a grid of 60 by 60 values by cubic splines
along each axis in turn, so that its contours are smooth curves and not
the polygons of the coarse grid. A grid that is not complete, or that
has fewer than four values of a label, is returned as it is.

## Usage

``` r
profile_refine(d)
```

## Arguments

- d:

  A data frame with the columns `a`, `b` and `profile`.

## Value

A data frame with the same columns.
