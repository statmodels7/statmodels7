# Where an Intercept on a Non-Identity Chart Stops Being an Estimate

The absolute value, on the link scale, past which
[`statmod_intercepts()`](https://statmodels7.github.io/statmodels7/reference/statmod_intercepts.md)
reads an intercept-only estimate as a limit the fit ran to: 16, a
parameter beyond \\e^{16} \approx 8.9 \times 10^6\\ or below its
reciprocal on a log chart.

## Usage

``` r
intercept_chart_limit()
```

## Value

A single positive number.
