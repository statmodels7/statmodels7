# Why a Term Is Outside the Fitting Scheme

Returns the reason a term cannot be assembled as a fixed design block,
or the empty string when it can.

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

Two reasons, and both are read from the term rather than from its class.

A covariance label
([`modelterms7::term_tag()`](https://statmodels7.github.io/modelterms7/reference/term_tag.html))
under a structural term of the **likelihood** shape says that a latent
the likelihood integrates out shares a prior with coefficients that are
estimated. The two are integrated by different routes and there is no
one prior to share, so the term is rejected rather than fitted as though
the label were absent, which would be a different model reported under
the name of the one that was asked for. A label under a **filter** is
fitted: its parameters are numbers estimated beside the coefficients,
and the block is read among them.

Whether a class is admissible is not a question about one term, and is
asked where every member is visible, by
[`class_space()`](https://statmodels7.github.io/statmodels7/reference/class_space.md).

A structural term implementing neither shape of the contract is rejected
for the reason its message gives.

## See also

[`reject_unfittable()`](https://statmodels7.github.io/statmodels7/reference/reject_unfittable.md)
