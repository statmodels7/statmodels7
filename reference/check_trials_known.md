# Refuse a Quantity That Needs Trials the Rows Do Not Carry

Signals an error where the family's numbers of trials are missing, which
is what
[`statmod_respec()`](https://statmodels7.github.io/statmodels7/reference/statmod_respec.md)
leaves when a prediction is asked at rows with no response.

## Usage

``` r
check_trials_known(distrib, what)
```

## Arguments

- distrib:

  The family, as
  [`statmod_respec()`](https://statmodels7.github.io/statmodels7/reference/statmod_respec.md)
  returns it.

- what:

  The quantity asked for, in words, for the message.

## Value

`distrib`, unchanged.

## Details

A parameter of a binomial family, its probability, is a function of the
linear predictor alone and is predicted at any rows. A moment of the
response or a prediction interval for a new observation depends on the
number of trials, and at rows that carry no `cbind(successes, failures)`
there is none to read.

## See also

[`check_trials()`](https://statmodels7.github.io/statmodels7/reference/check_trials.md),
[`statmod_respec()`](https://statmodels7.github.io/statmodels7/reference/statmod_respec.md),
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
