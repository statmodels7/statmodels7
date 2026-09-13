# The Trace Against a Contraction, Without Forming It

\\\mathrm{tr}(M\\T)\\ where \\T\\ is the third derivative of the
penalized objective contracted once, or the fourth contracted twice.

## Usage

``` r
trace_design_form(spec, G, deriv, params, npar, tv, tu = NULL, key = NULL)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- G:

  The per-observation diagonals, from
  [`block_leverage()`](https://statmodels7.github.io/statmodels7/reference/block_leverage.md).

- deriv:

  The third or fourth derivatives on the link scale.

- params, npar:

  The parameter names and block sizes.

- tv:

  The direction the derivative is contracted in.

- tu:

  A second direction, for a fourth derivative; `NULL` for a third.

- key:

  A function of the three or four parameter positions returning the
  component of `deriv` to read, or `NULL` for the observed route's keys,
  which are symmetric in every position. On the expected route `deriv`
  is the derivative of the expected information, symmetric in its first
  pair and in its second separately, and carries its own builder.

## Value

A single number.

## Details

[`contract3()`](https://statmodels7.github.io/statmodels7/reference/contract3.md)
and
[`contract4()`](https://statmodels7.github.io/statmodels7/reference/contract4.md)
assemble \\-X_a'\\\mathrm{diag}(\omega_i w_i)\\X_b\\ block by block, and
where the result is only ever traced against \\M\\ the assembly is
waste: the trace of that block against \\M\_{ab}\\ is \\-\sum_i \omega_i
w_i G\_{ab,i}\\. The block and its transpose both contribute, so an
off-diagonal pair counts twice.

The weights \\w_i\\ are built exactly as those two functions build them,
and that keeps the two routes the same arithmetic in the same order
where it matters.

## See also

[`contract3()`](https://statmodels7.github.io/statmodels7/reference/contract3.md),
[`contract4()`](https://statmodels7.github.io/statmodels7/reference/contract4.md)
