# The Penalty a Coordinate Descent Builds Its Table From

Restricts a kinked penalty to the coordinates the strong rule kept, and
writes into a scaled SCAD or MCP the curvature of the current step.

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

  The column curvatures, one per coordinate of the block, with at least
  the kept ones filled in.

- prev:

  `NULL`, or the curvature the previous table used, on the scale of
  \\u\\, one per coordinate of the block; the new one is its geometric
  mean with the current step's.

## Value

A penalty over `length(keep)` coordinates, or `pen` itself where nothing
needs to change.

## Details

The compiled descent reads row \\a\\ of the table for coordinate
\\k_a\\, the \\a\\-th kept one. A table built from the whole penalty
pairs row \\a\\ with coordinate \\a\\ instead, so the restriction is
what makes the two agree: every per-coordinate property is subset, the
map's diagonal and the curvature among them. It is done only where a
coordinate was screened out or a curvature is written, and a penalty
under a map that is not diagonal is returned as it stands, having no
table.

A SCAD or MCP carrying a curvature is scaled SELF-CONSISTENTLY: at the
point the descent settles at, its curvature is the one the step is taken
with, \\c_j = v_j\\, the (centered) weighted sum of squares of the
column at the current working weights, divided by \\d_j^2\\ under a
diagonal map, the penalty being written on \\u = D\beta\\. On the way
there it is DAMPED, \\c_j \leftarrow \sqrt{c_j^{\mathrm{prev}} v_j}\\:
taken undamped, a coefficient between the two knees of a logistic SCAD
alternated between 1.009 and 1.496 with its curvature between 21.8 and
30.7 and never settled, 3 fits of 16 not converging on a strong-effect
probe; damped, all 16 converge with the KKT conditions at the fit's
curvature met to \\5 \times 10^{-8}\\. The scaled step is
\\\sqrt{c^{\mathrm{prev}}\_j / v_j}\\, one where the weights have not
moved, so the step condition \\t c_j \< a - 1\\ (SCAD) or \\t c_j \<
\gamma\\ (MCP) holds unless they moved by a factor of \\(a-1)^2\\.
[`statmod_curv()`](https://statmodels7.github.io/statmodels7/reference/statmod_curv.md)
writes the undamped quantity into the specification at every pass of the
alternation, so the objective and the degrees of freedom read it too.

## See also

[`coord_fit()`](https://statmodels7.github.io/statmodels7/reference/coord_fit.md),
[`statmod_curv()`](https://statmodels7.github.io/statmodels7/reference/statmod_curv.md)
