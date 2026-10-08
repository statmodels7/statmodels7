# Where the Break-Point Term of a Fit Is

Finds the break-point term a test reads: its equation, the prefix of its
coefficient names, the built term, and whether its step is sharp. The
term is looked for among the terms of every equation and among the
sub-terms of the parameters of a
[`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html)
term.

## Usage

``` r
bp_test_locate(spec, param = NULL, term = NULL)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- param, term:

  As
  [`statmod_breakpoint_test()`](https://statmodels7.github.io/statmodels7/reference/statmod_breakpoint_test.md).

## Value

A list with `param`, `name` (the term's name in its equation), `prefix`,
`term` and `sharp`.
