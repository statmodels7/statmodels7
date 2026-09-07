# Name the Aliased Coefficients of a Fit

Turns the coordinates the scoring step's pivot dropped into the labels a
reader sees.

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
which column is left out is settled by the pivot, exactly as
[`lm()`](https://rdrr.io/r/stats/lm.html) and
[`glm()`](https://rdrr.io/r/stats/glm.html) settle it. Those columns
come back from the solve as coordinates of the coefficient vector; here
they become the labels
[`coef_labels()`](https://statmodels7.github.io/statmodels7/reference/coef_labels.md)
gives, so that nothing downstream has to rebuild a design to say which
coefficient is which.

A coordinate a penalty covers is never among them, and that falls out of
where the test is made rather than being written as a rule: the pivot
runs on the AUGMENTED system, design and penalty factor together, so a
column the design alone does not identify is identified there and is not
dropped. Measured on two identical columns, `ridge(~ 0 + x + v)` fits
and splits the effect evenly between them, at 0.421897 each, with a full
variance matrix.

## See also

[`coef_labels()`](https://statmodels7.github.io/statmodels7/reference/coef_labels.md),
which supplies the names, and
[`vcov()`](https://rdrr.io/r/stats/vcov.html), which holds these
coordinates rather than refusing the whole matrix.
