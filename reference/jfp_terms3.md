# The Terms of the Third and Fourth Derivatives

Returns the terms of \\D_3(v)\\ (`jfp_terms3()`) and \\D_4(v, w)\\
(`jfp_terms4()`), summed over the blocks with their powers, in the form
that
[`jfp_terms_matrix()`](https://statmodels7.github.io/statmodels7/reference/jfp_terms_matrix.md)
assembles and
[`jfp_terms_trace()`](https://statmodels7.github.io/statmodels7/reference/jfp_terms_matrix.md)
traces.

## Usage

``` r
jfp_terms3(st, v)

jfp_terms4(st, v, w)
```

## Arguments

- st:

  A state, from
  [`jfp_state()`](https://statmodels7.github.io/statmodels7/reference/jfp_state.md).

- v, w:

  Directions over the joint vector.

## Value

A list with `sparse`, a list of per-observation weights keyed as
[`statmod_information_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_information_at.md)
reads its `H` (the matrix is
\\\sum\_{c,d}X_c^\top\mathrm{diag}(w\\W\_{cd})X_d\\), and `lr`, a list
of triples `X`, `K`, `Y` (the matrix is \\XKY^\top\\).
