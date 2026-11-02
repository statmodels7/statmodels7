# The Fourth-Order Blocks of One Filter

Returns the callback, with its data as the attribute `"data"`, through
which
[`modelterms7::term_fourth()`](https://statmodels7.github.io/modelterms7/reference/term_fourth.html)
reads the family's derivatives to the fifth order in
[`structural_chain_extra2()`](https://statmodels7.github.io/statmodels7/reference/structural_chain_extra2.md).

## Usage

``` r
structural_blocks4(spec, jd, st)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- jd:

  The joint rows, from
  [`joint_design_rows()`](https://statmodels7.github.io/statmodels7/reference/joint_design_rows.md).

- st:

  The shared quantities, from
  [`structural_grad_parts()`](https://statmodels7.github.io/statmodels7/reference/structural_grad_parts.md).

## Value

The callback, or `NULL` where the family has no fifth derivative.
