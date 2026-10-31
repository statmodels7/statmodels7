# Assemble or Trace a List of Derivative Terms

Assemble or Trace a List of Derivative Terms

## Usage

``` r
jfp_terms_matrix(st, terms)

jfp_terms_trace(st, terms, M)
```

## Arguments

- st:

  A state, from
  [`jfp_state()`](https://statmodels7.github.io/statmodels7/reference/jfp_state.md).

- terms:

  What
  [`jfp_terms3()`](https://statmodels7.github.io/statmodels7/reference/jfp_terms3.md)
  or
  [`jfp_terms4()`](https://statmodels7.github.io/statmodels7/reference/jfp_terms3.md)
  returns.

- M:

  A symmetric matrix over the joint vector.

## Value

`jfp_terms_matrix()` a symmetric base matrix over the joint vector;
`jfp_terms_trace()` the trace of that matrix against `M`, a single
number.
