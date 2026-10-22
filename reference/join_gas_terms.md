# Join the Score-Driven Terms of Several Equations

Replaces the
[`modelterms7::gas()`](https://statmodels7.github.io/modelterms7/reference/gas.html)
terms of a model, when there is more than one, with a single
[`modelterms7::GasJoint()`](https://statmodels7.github.io/modelterms7/reference/GasJoint.html)
placed in the equation of the first of them under the key `"gas_joint"`.

## Usage

``` r
join_gas_terms(terms)
```

## Arguments

- terms:

  The built terms, a named list of named lists, one per distribution
  parameter.

## Value

`terms`, with the score-driven terms joined where there are two or more.

## Details

[`reject_incompatible()`](https://statmodels7.github.io/statmodels7/reference/reject_incompatible.md)
has already checked that each equation carries at most one such term and
that no other structural term sits beside them.
[`modelterms7::gas_joint()`](https://statmodels7.github.io/modelterms7/reference/gas_joint.html)
checks that the terms walk the same series and share one `scaling`, and
fills the cross loadings of the terms that share a `cross` label.
