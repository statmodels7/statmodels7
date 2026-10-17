# Why a Term Is Outside the Fitting Scheme

Returns the reason why a term cannot be assembled as a fixed design
block, or the empty string when it can.

## Usage

``` r
unfittable_reason(term)
```

## Arguments

- term:

  One built term.

## Value

A single string.

## Details

There are two reasons, and both are read from the term rather than from
its class.

A covariance label
([`modelterms7::term_tag()`](https://statmodels7.github.io/modelterms7/reference/term_tag.html))
under a structural term of the likelihood shape would make a latent that
the likelihood integrates out share a prior with coefficients that are
estimated. The two are integrated by different routes and there is no
one prior to share, so the term is rejected. Fitting it as though the
label were absent would report a different model under the name of the
one requested. A label under a filter is fitted: its parameters are
numbers estimated beside the coefficients, and the block is read among
them.

Whether a class is admissible depends on all its members, and
[`class_space()`](https://statmodels7.github.io/statmodels7/reference/class_space.md)
determines it where every member is visible.

A structural term that implements neither `term_filter()` nor
`term_loglik()` is rejected for the reason that its message gives.

## See also

[`reject_unfittable()`](https://statmodels7.github.io/statmodels7/reference/reject_unfittable.md)
