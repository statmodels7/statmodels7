# The Range of a Break-Point Position

Returns the smallest and the largest value that the position of a
break-point can be held at before the term stops depending on it: the
range of the variable the break-point is on.

## Usage

``` r
break_edges(spec, design, param, j)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  Its design.

- param:

  The distribution parameter of the coefficient.

- j:

  The position of the coefficient in the parameter's equation.

## Value

A numeric vector of two, `-Inf` and `Inf` where there is no range.

## Details

A position at or beyond the largest observation makes every hinge \\(x -
\psi)\_+\\ and every step \\1(x \> \psi)\\ of the term zero, and a
position at or below the smallest makes the hinge a line that the term's
linear part absorbs. Either way the restricted refit is the same for
every such value, so that an interval by inversion can stop at the end
of the range. A coefficient that is not the position of a break-point of
a
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html)
term, or a term whose positions are developed over covariates, has no
range.

## See also

[`statmod_invert()`](https://statmodels7.github.io/statmodels7/reference/statmod_invert.md),
which stops its search there.
