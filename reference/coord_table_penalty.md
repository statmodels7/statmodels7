# The Penalty from Which a Coordinate Descent Builds Its Table

Restricts a kinked penalty to the coordinates that the strong rule kept,
and writes into a scaled SCAD or MCP the curvature of the current step.

## Usage

``` r
coord_table_penalty(pen, keep, v, prev = NULL)
```

## Arguments

- pen:

  The block's penalty.

- keep:

  The kept coordinates, one-based.

- v:

  The column curvatures, one per coordinate of the block, where at least
  the kept ones are filled in.

- prev:

  `NULL`, or the curvature that the previous table used, on the scale of
  \\u\\, one per coordinate of the block; the new one is its geometric
  mean with the current step's.

## Value

A penalty over `length(keep)` coordinates, or `pen` itself where the
penalty can stay as it is.

## Details

The compiled descent reads row \\a\\ of the table for coordinate
\\k_a\\, the \\a\\-th kept one. A table built from the whole penalty
pairs row \\a\\ with coordinate \\a\\ instead, so the two agree once the
penalty is restricted: every per-coordinate property is subset,
including the map's diagonal and the curvature. The restriction is done
only where a coordinate was screened out or a curvature is written, and
a penalty under a map that is not diagonal is returned as it stands,
having no table.

A SCAD or MCP carrying a curvature is scaled self-consistently: at the
point at which the descent settles, its curvature is the one with which
the step is taken, \\c_j = v_j\\, the (centered) weighted sum of squares
of the column at the current working weights, divided by \\d_j^2\\ under
a diagonal map, because the penalty is written on \\u = D\beta\\. On the
way there the curvature is damped, \\c_j \leftarrow
\sqrt{c_j^{\mathrm{prev}} v_j}\\. Taken undamped, the coefficient and
the curvature of a coordinate between the two knees of the penalty can
alternate between two values without settling, whereas the damped update
converges. The product of the step and the curvature is
\\\sqrt{c^{\mathrm{prev}}\_j / v_j}\\, which is one where the weights
have not moved, so the step condition \\t c_j \< a - 1\\ (SCAD) or \\t
c_j \< \gamma\\ (MCP) holds unless the weights moved by a factor of
\\(a-1)^2\\.
[`statmod_curv()`](https://statmodels7.github.io/statmodels7/reference/statmod_curv.md)
writes the undamped quantity into the specification at every pass of the
alternation, so the objective and the degrees of freedom read it too.

## See also

[`coord_fit()`](https://statmodels7.github.io/statmodels7/reference/coord_fit.md),
[`statmod_curv()`](https://statmodels7.github.io/statmodels7/reference/statmod_curv.md)
