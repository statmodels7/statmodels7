# S7 Class for a Starting-Value Strategy

Where a fit begins, as an object rather than as a vector of numbers.

## Usage

``` r
start_strategy(label = character(0))
```

## Arguments

- label:

  A short name, used when printing.

## Value

An S7 object of class `start_strategy`.

## Details

The reason it is an object is that the interesting answers are not
values but PROCEDURES: draw around the intercept-only fit, or search the
likelihood for a basin. A strategy is asked once, before the alternation
between the coefficients and the hyperparameters begins, which is what
separates it from an optimizer: an optimizer runs at every step of the
fit, a strategy runs once at its start.

## See also

[`start_intercepts`](https://statmodels7.github.io/statmodels7/reference/start_intercepts.md),
[`start_origin`](https://statmodels7.github.io/statmodels7/reference/start_origin.md),
[`start_random`](https://statmodels7.github.io/statmodels7/reference/start_random.md)

## Examples

``` r
start_origin()
#> <start> zero
```
