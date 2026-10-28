# Fitted Positions of the Break-Points Next to One

Returns the fitted positions of the break-points of the same term that
lie immediately below and above the break-point whose position is the
coefficient `j`, which bound the search for an interval by inversion.

## Usage

``` r
neighbor_positions(fit, spec, design, param, j)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- spec, design:

  The specification and its design.

- param:

  The distribution parameter of the coefficient.

- j:

  The position of the coefficient in the parameter's equation.

## Value

A numeric vector of two, the nearest position below and the nearest
above, `NA` where there is none.

## Details

Only a position of a
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html)
term with two or more break-points and positions that carry no
development has neighbors. A neighbor is a position strictly below or
above the estimate, read at the confinement of the term.

## See also

[`statmod_invert()`](https://statmodels7.github.io/statmodels7/reference/statmod_invert.md),
which stops its search there, and
[`break_edges()`](https://statmodels7.github.io/statmodels7/reference/break_edges.md).
