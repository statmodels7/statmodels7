# Where a Specification's Held Coefficients Sit in the Stacked Vector

Translates `spec@held_coef`, which names coefficients, into positions in
the stacked coefficient vector the objective works on, with the values
beside them.

## Usage

``` r
held_positions(spec, design, obj, beta)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md),
  read for `held_coef`.

- design:

  Its design, read for each equation's coefficient names.

- obj:

  The objective, read for its split of the stacked vector.

- beta:

  The stacked coefficients, read for their length alone.

## Value

A list of two numeric vectors of equal length, `where` (integer
positions) and `value`. Both are empty where nothing is held.

## Details

A hold is written by name because a name is what a caller can say and a
position is not: the design's column order is the formula's business,
and a term added to an equation moves every coordinate after it. The
translation is the objective's own split of the stacked vector, so the
answer is in the numbering every consumer of that vector already uses.

A name that matches no coefficient of its equation is an error rather
than a hold nothing enforces, which is what a silent `NA` would have
been.

## See also

[`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md),
which enforces the hold.
