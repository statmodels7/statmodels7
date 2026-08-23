# The Column a Simulated Response Is Written To

The name on the formula's left-hand side, which must be a symbol.

## Usage

``` r
rstatmod_response_name(response)
```

## Arguments

- response:

  The formula's left-hand side.

## Value

A single string.

## Details

A transformed response is rejected rather than answered. Under
`log(y) ~ x` the model generates values of `log(y)`, so a column called
`y` would hold the wrong quantity and one called `"log(y)"` would not be
the name the formula reads back; the earlier version wrote the first of
those silently. A censored response is rejected for the reason
[`statmod`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
rejects one – there is no censored likelihood to fit it back with.

## See also

[`rstatmod`](https://statmodels7.github.io/statmodels7/reference/rstatmod.md)
