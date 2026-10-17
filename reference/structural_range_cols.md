# Which of a Structural Term's Free Parameters a Penalty Covers

Returns the positions, among a structural term's free parameters, that
some penalty of the model shrinks. Those are the directions over which
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
integrates; the rest are profiled.

## Usage

``` r
structural_range_cols(spec, design, key, free)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- key:

  The term's key in the specification, its call as written.

- free:

  The term's free parameter names, in the order in which the term holds
  them.

## Value

An integer vector of positions into `free`, sorted and without repeats.
The vector is empty if the term is not covered by any penalty.

## Details

A term declares its penalties through
[`modelterms7::term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.html),
and each entry names a subset of the term's own parameters. This
collects those subsets and maps them onto positions in `free`.

## See also

[`statmod_marginal_full()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_full.md),
the caller.
