# The Positions of One Break-Point Term Among the Held Coefficients

Finds, for each continuous
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html)
term without subformulas that has two or more of its positions among the
labelled coefficients, the indices into `which` of those positions in
the order of the term's break-points.

## Usage

``` r
profile_psi_groups(spec, design, param, coefname)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  Its design.

- param, coefname:

  The equation and the name of each labelled coefficient.

## Value

A list of integer vectors, possibly empty.
