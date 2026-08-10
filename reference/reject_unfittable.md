# Reject a Term the Fitting Scheme Does Not Cover

Signals an error naming any term whose block is not a fixed design,
which is what the alternation of
[`statmod`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
assembles.

## Usage

``` r
reject_unfittable(terms)
```

## Arguments

- terms:

  The built terms, a named list of named lists, one per distribution
  parameter.

## Value

`NULL`, invisibly; called for the error.

## Details

Two shapes are outside that assembly, and both are read off the term
rather than from a list of class names, so a term written later is
covered without an edit here.

A **structural** term rewrites the likelihood instead of contributing a
predictor, so it has no design block at all and answers neither
`term_matrix()` nor `term_npar()`. Reaching it through the design
produced an error naming one of those generics, which says nothing about
the cause.

A term whose block **depends on its own coefficients** registers a
[`term_refresh`](https://statmodels7.github.io/modelterms7/reference/term_refresh.html)
method of its own, the base method on `model_term` being the identity;
the class a method was registered on is `attr(m, "signature")[[1]]`. For
those the block is a Jacobian and the working solution is an increment,
so assembling it once and solving for the coefficients estimates
something else while reporting convergence: measured on `seg()`, the
break-point stays at its starting value and the fitted mean of a
continuous construction carries a step. Rejecting them is what keeps
that out of a returned object until the alternation refreshes a block
between inner fits.

Every equation is examined before the error is raised, so a model
carrying one such term in the mean and another in the scale reports both
rather than the first.

## See also

[`statmod_terms`](https://statmodels7.github.io/statmodels7/reference/statmod_terms.md),
[`statmod`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
