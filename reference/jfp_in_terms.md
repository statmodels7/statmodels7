# Whether a Model Carries a Jeffreys Penalty

`jfp_in_terms()` returns whether a term of an equation, or a sub-term of
one of its subformulas, carries a
[`modelterms7::JeffreysPenalty()`](https://statmodels7.github.io/modelterms7/reference/JeffreysPenalty.html),
read from
[`modelterms7::term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.html).
`spec_jeffreys()` returns the property `jeffreys` of the specification,
which
[`statmod_spec()`](https://statmodels7.github.io/statmodels7/reference/statmod_spec.md)
sets with `jfp_in_terms()`, and reads the terms where the property is
`NA`.

## Usage

``` r
jfp_in_terms(terms)

spec_jeffreys(spec)
```

## Arguments

- terms:

  A named list with one entry per parameter, each a named list of built
  terms.

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

## Value

A single logical.

## Details

The answer depends on which penalties the terms declare and not on their
coefficients, so a refresh, which replaces a term by the same term at
other coefficients, leaves it valid.
[`jfp_blocks()`](https://statmodels7.github.io/statmodels7/reference/jfp_blocks.md)
reads it before enumerating the penalized units, which the objective
would otherwise do at every evaluation on a design that keeps no memo of
them.
