# The Parameters a Structural Term Starts From

What the term itself declares through
[`term_start`](https://statmodels7.github.io/modelterms7/reference/term_start.html):
as near the model without the term as its charts allow, which only the
term can say – a score loading on the log chart has no coordinate for
zero, and starts at a weak response instead.

## Usage

``` r
structural_zeta_start(term, target = NULL)
```

## Arguments

- term:

  A built structural term.

- target:

  The response on the scale of the term's equation, from
  [`predictor_target`](https://statmodels7.github.io/statmodels7/reference/predictor_target.md),
  or `NULL`. A term whose start is read off the data – the marginal
  break-point term's exact profile – consumes it; every other method
  ignores it through the dots.

## Value

A named numeric vector.
