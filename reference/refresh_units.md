# The Terms Whose Block Moves With Their Coefficients

Returns one entry per refreshable term, carrying everything that its
derivative needs: the distribution parameter of the equation in which
the term sits, its columns within that equation and in the stacked
coefficient vector, its own coefficients and the term rebuilt at them.

## Usage

``` r
refresh_units(spec, design, coef, params, npar, offs)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design, already refreshed at `coef`.

- coef:

  The coefficients.

- params, npar, offs:

  The block bookkeeping.

## Value

A list, which is empty if the design has no block that moves. Each
element is a list with `a`, the index of the term's equation; `param`,
the distribution parameter of that equation; `term`, the name of the
term within the equation; `cols`, its columns within the equation; `ra`,
the same columns in the stacked coefficient vector; `bt`, the term's
coefficients; and `built`, the term refreshed at `bt`.

## Details

The entries are resolved once and passed to every consumer, so the
gradient's contraction and the Hessian's three corrections cannot
disagree about which terms move or where their columns are.

## See also

[`u_refresh()`](https://statmodels7.github.io/statmodels7/reference/u_refresh.md),
[`contract3_refresh()`](https://statmodels7.github.io/statmodels7/reference/contract3_refresh.md)
