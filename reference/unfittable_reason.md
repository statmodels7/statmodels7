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

A term carrying a covariance label
([`modelterms7::term_tag()`](https://statmodels7.github.io/modelterms7/reference/term_tag.html))
says that its coefficients share a block with those of other terms,
which this layer cannot yet build: the penalty would have to read
columns from more than one equation, and every enumeration here
addresses a penalty by the pair of its parameter and its key. The term
is rejected rather than fitted as though the label were absent, which
would be a different model reported under the name of the one that was
asked for.

A structural term implementing neither shape of the contract is rejected
for the reason its message gives.

## See also

[`reject_unfittable()`](https://statmodels7.github.io/statmodels7/reference/reject_unfittable.md)
