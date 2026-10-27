# The Status of a Replicate of a Special Model

Returns 0 for a replicate that is kept and 2 for a replicate that did
not reach a mode.

## Usage

``` r
boot_special_status(mode, spec_b, r)
```

## Arguments

- mode:

  `"structural"` or `"refresh"`.

- spec_b:

  The specification of the replicate.

- r:

  The refit.

## Value

A single integer.

## Details

A model with a score-driven term is checked on the Newton decrement over
the coefficients and the term's parameters,
[`filter_decrement()`](https://statmodels7.github.io/statmodels7/reference/filter_decrement.md),
against
[`mode_error_limit()`](https://statmodels7.github.io/statmodels7/reference/mode_error_limit.md).
A model with a regime term, or with a term whose block moves with its
coefficients, has no information that reads the mode (the block is a
working linearization), and its replicate is kept when every coefficient
and the objective are finite.
