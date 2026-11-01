# The Weighted Cross Product of the Joint Rows

Assembles \\\sum\_{c,d}V_c^\top\mathrm{diag}(w\\W\_{cd})V_d\\ over the
joint vector of a model with a filter, where \\V_c\\ are the rows of
every equation in the joint vector
([`jfp_joint_rows()`](https://statmodels7.github.io/statmodels7/reference/jfp_joint_rows.md))
and the weights are keyed as
[`jfp_terms3()`](https://statmodels7.github.io/statmodels7/reference/jfp_terms3.md)
returns them.

## Usage

``` r
jfp_joint_sparse(st, sparse)
```

## Arguments

- st:

  A state, from
  [`jfp_state()`](https://statmodels7.github.io/statmodels7/reference/jfp_state.md),
  of a model with a filter.

- sparse:

  The list of per-observation weights.

## Value

A base matrix over the joint vector.
