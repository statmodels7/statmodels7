# The Design of the Terms Holding Nested Effects, at Other Coefficients

Recomputes the columns of every term holding a nested random effect at
the coefficients given, which is the Jacobian there for a term whose
block is one.

## Usage

``` r
nested_design_at(spec, design, coef, terms)
```

## Arguments

- spec, design:

  The specification and its design.

- coef:

  The coefficients, a named list.

- terms:

  A data frame with columns `param` and `term`.

## Value

The design with those columns replaced.
