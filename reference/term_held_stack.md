# The Stacked Positions the Terms Hold

The positions, in the stacked coefficient vector, of the coefficients
the model's terms hold
([`modelterms7::term_held()`](https://statmodels7.github.io/modelterms7/reference/term_held.html)):
the slot of a held break-point. The terms are read from the design's
state where the fit carries one, so a term held during the fit is seen
at once, and from the specification otherwise.

## Usage

``` r
term_held_stack(spec, design)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  Its design.

## Value

An integer vector of positions, possibly empty.
