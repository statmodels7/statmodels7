# Name the Aliased Coefficients of a Fit

Converts the coordinates that the pivot of the scoring step dropped into
coefficient labels.

## Usage

``` r
aliased_labels(spec, design, idx)
```

## Arguments

- spec:

  A `statmod_spec`.

- design:

  The design, as
  [`statmod_design()`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
  builds it.

- idx:

  Integer coordinates of the coefficient vector, as the alternation
  reports them. `NULL` or empty gives `character(0)`.

## Value

A character vector of coefficient labels, possibly empty.

## Details

A design of less than full rank has no estimate for the columns that
repeat information already carried: only combinations are estimable, and
the pivot determines which column is left out, as in
[`lm()`](https://rdrr.io/r/stats/lm.html) and
[`glm()`](https://rdrr.io/r/stats/glm.html). Those columns come back
from the solve as coordinates of the coefficient vector; this function
converts them into the labels that
[`coef_labels()`](https://statmodels7.github.io/statmodels7/reference/coef_labels.md)
supplies, so that no later step has to rebuild a design to identify a
coefficient.

A coordinate that a penalty covers is never among them, because the
pivot runs on the augmented system (design and penalty factor together),
in which a column that the design alone does not identify is identified.
Under `ridge(~ 0 + x + v)` with two identical columns, the effect is
split evenly between the columns and the variance matrix is complete.

## See also

[`coef_labels()`](https://statmodels7.github.io/statmodels7/reference/coef_labels.md),
which supplies the names, and
[`vcov()`](https://rdrr.io/r/stats/vcov.html), which reports these
coordinates as missing and inverts the rest of the matrix.
