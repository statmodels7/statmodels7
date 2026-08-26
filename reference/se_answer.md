# The Shape a Prediction With Its Uncertainty Comes Back In

The rows
[`predict_se()`](https://statmodels7.github.io/statmodels7/reference/predict_se.md)
computed, reduced to what was asked for.

## Usage

``` r
se_answer(su, what, params, spec)
```

## Arguments

- su:

  The per-parameter tables.

- what:

  What was asked for.

- params:

  The distribution's parameters.

- spec:

  The specification, for the message.

## Value

A data frame, or a named list of them.

## Details

The vocabulary is the one a prediction without uncertainty answers,
minus the moments: a moment's delta method needs its derivative in every
parameter, which distributions7 does not offer as a generic, and an
interval assembled from what is at hand instead would be a number nobody
could check.
