# Hold What a Simulation's Caller Named

Writes the entries of `par` over the drawn truth, each addressed by a
name the result reports it under.

## Usage

``` r
rstatmod_named(coef, psi, design, params, par)
```

## Arguments

- coef:

  The drawn coefficients.

- psi:

  The drawn structural parameters, or `NULL`.

- design:

  The design.

- params:

  The distribution parameter names.

- par:

  A named list, or `NULL`.

## Value

A list with `coef` and `psi`.

## Details

A single namespace covers the whole model, so the caller names the
quantities to fix and everything else stays drawn. A key may be a
distribution parameter, which is that whole equation; one of its
coefficients or a group of them, written `parameter.coefficient`; a
structural term's own parameter; or a group of those. A group is a name
that its members extend at a dot, so `omega.random` is every deviation
in a development of `omega` and `mu.s(x)` is every coordinate of that
smooth.

A value may be a vector of the group's own length, a single number used
for all of them, or a function of the count. The function is how a
structured truth is written: `function(k) rnorm(k, 0, 0.4)` is a random
effect at a scale of its own and `function(k) c(2, -1.5, rep(0, k - 2))`
is a sparse truth for a lasso to find.

A structural parameter is named on the scale of
[`modelterms7::term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.html):
a loading is the loading and not its logarithm.

## See also

[`rstatmod()`](https://statmodels7.github.io/statmodels7/reference/rstatmod.md),
[`rstatmod_truth()`](https://statmodels7.github.io/statmodels7/reference/rstatmod_truth.md)
