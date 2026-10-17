# The Degree of Each Free Value of a Chart in the Scale of Its Matrix

Computes \\k_j \in \\-1, 0, 1\\\\ for each free value: 1 when
multiplying the matrix by \\\lambda^2\\ multiplies the free value by
\\\lambda\\, -1 when it divides it by \\\lambda\\, and 0 otherwise (a
log-scale coordinate, which shifts, or a correlation, which does not
move).

## Usage

``` r
chart_homogeneity(chart)
```

## Arguments

- chart:

  A parameters7 matrix chart.

## Value

An integer vector of length `chart@n_free`.
