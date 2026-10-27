# The Quantities of a Structural Term at Each Replicate

Evaluates, for every bootstrap replicate, the quantities that a
structural term reports
([`modelterms7::term_readable()`](https://statmodels7.github.io/modelterms7/reference/term_readable.html))
from the replicate's free parameters of the term.

## Usage

``` r
boot_term_values(object, B, spec, design)
```

## Arguments

- object:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- B:

  The matrix of replicates of
  [`boot_replicate_matrix()`](https://statmodels7.github.io/statmodels7/reference/boot_replicate_matrix.md).

- spec, design:

  The specification and design of the fit.

## Value

A named list, one entry for each structural term, holding a matrix with
one row for each replicate and one column for each quantity of the term,
in the order of
[`modelterms7::term_readable()`](https://statmodels7.github.io/modelterms7/reference/term_readable.html).
A replicate with a missing parameter of the term gives a row of missing
values.
