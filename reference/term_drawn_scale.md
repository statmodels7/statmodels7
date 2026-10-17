# The Width at Which a Term Drew a Penalty's Coordinates

Returns one number where every coordinate that a penalty covers was
written by the term that owns it, and `NULL` otherwise.

## Usage

``` r
term_drawn_scale(tg, ok, owned)
```

## Arguments

- tg:

  The unit's targets, as
  [`unit_draw_targets()`](https://statmodels7.github.io/statmodels7/reference/unit_draw_targets.md)
  returns them.

- ok:

  Which of them the prior's draw reached.

- owned:

  The `owned` element of
  [`rstatmod_term_draw()`](https://statmodels7.github.io/statmodels7/reference/rstatmod_term_draw.md).

## Value

A single number, or `NULL`.

## Details

A penalty whose coordinates were overwritten by a term no longer
describes them, so the truth reported for it is the width at which the
term drew. The test is that every covered coordinate was owned: a
penalty spanning owned and unowned ones together is described by neither
width, and reporting either would name a truth that part of the
coefficients do not have.

## See also

[`rstatmod_term_draw()`](https://statmodels7.github.io/statmodels7/reference/rstatmod_term_draw.md)
