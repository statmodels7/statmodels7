# The Leibniz Rule Over Three Factors

Returns \\\sum\_{R = S_1\sqcup S_2\sqcup S_3}\partial\_{S_1}f\\ P(S_2,
S_3)\\ for a set \\R\\ of at most two directions, written out case by
case: one term for no direction, three for one, nine for two.

## Usage

``` r
jfp_leibniz(jq, R, P)
```

## Arguments

- jq:

  A function of a set of directions returning the derivative of the
  scalar factor along them, a numeric vector over the observations.

- R:

  A character vector of at most two direction names.

- P:

  A function of two sets of directions returning the product of the
  other two factors differentiated along them, or `NULL` where it is
  zero.

## Value

The sum, or `NULL` where every term is zero.
