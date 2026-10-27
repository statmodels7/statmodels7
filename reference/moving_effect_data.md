# The Data of the Effect of a Term with a Moving Block

The counterpart of
[`effect_data()`](https://statmodels7.github.io/statmodels7/reference/effect_data.md)
for
[modelterms7::SegTerm](https://statmodels7.github.io/modelterms7/reference/SegTerm.html)
and
[modelterms7::NlTerm](https://statmodels7.github.io/modelterms7/reference/NlTerm.html).
The contribution on the grid is
[`modelterms7::term_value()`](https://statmodels7.github.io/modelterms7/reference/term_value.html),
not the block times the coefficients, and the band is the delta method
with the block of
[`modelterms7::term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.html),
which is the derivative of the contribution in the coefficients.

## Usage

``` r
moving_effect_data(type, tm, fit, ctx)
```

## Arguments

- type:

  A `PlotEffect` object.

- tm:

  The term.

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- ctx:

  The context described at
  [`plot_data()`](https://statmodels7.github.io/statmodels7/reference/plot_data.md).

## Value

As
[`effect_data()`](https://statmodels7.github.io/statmodels7/reference/effect_data.md).
