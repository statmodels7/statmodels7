# The Parameters a Structural Term Starts From

What the term itself declares through
[`modelterms7::term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.html):
as near the model without the term as its charts allow, which only the
term can determine. A score loading on the log chart has no coordinate
for zero, so it starts at a weak response instead.

## Usage

``` r
structural_zeta_start(term, target = NULL)
```

## Arguments

- term:

  A built structural term.

- target:

  The response on the scale of the term's equation, as
  [`predictor_target()`](https://statmodels7.github.io/statmodels7/reference/predictor_target.md)
  returns it, or `NULL`. A term whose start is read off the data uses
  it: the marginal break-point term (its exact profile) and the regime
  term (quantiles of the target) do so, and the other methods ignore it
  through the dots.

## Value

A named numeric vector.
