# A Structural Term's Parameters for a Simulation

The term's starting values with whatever the caller named written over
them.

## Usage

``` r
rstatmod_psi(tm, psi, given)
```

## Arguments

- tm:

  The built structural term.

- psi:

  Its parameters at the specification's own state.

- given:

  A named list, or `NULL`.

## Value

A named list.

## Details

The values are on the scale
[`modelterms7::term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.html)
names, which is the one a reader knows: a loading is the loading rather
than its logarithm, and a persistence is the partial autocorrelation its
chart carries rather than the autoregressive coefficient that chart
produces. A name the term does not have is reported with the ones it
does, since a misspelled parameter would otherwise leave the term at a
starting value and simulate a model with almost no dynamics.

## See also

[`rstatmod()`](https://statmodels7.github.io/statmodels7/reference/rstatmod.md)
