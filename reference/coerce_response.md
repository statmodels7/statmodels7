# Carry a Categorical Response onto the Numeric Scale

Returns the response as a number, converting a two-level `factor`,
`character` or `logical` to 0 and 1 and leaving anything else alone.

## Usage

``` r
coerce_response(response, levels = NULL)
```

## Arguments

- response:

  The evaluated left-hand side.

- levels:

  The mapping to reuse, from a previous call's attribute, or `NULL` to
  read it off this response.

## Value

`response` unchanged, or a numeric vector of 0 and 1 carrying the
attribute `"response_levels"`.

## Details

A binary outcome is ordinarily recorded as a factor, and
[`glm()`](https://rdrr.io/r/stats/glm.html) accepts one: its binomial
family reads the first level as the failure. The conversion makes the
response a number before it reaches the density, which would otherwise
receive a factor.

The first level is the failure, which is
[`glm()`](https://rdrr.io/r/stats/glm.html)'s rule, so a model moved
from [`glm()`](https://rdrr.io/r/stats/glm.html) keeps its signs. The
conversion applies to every family and not only to the binary ones,
which differs from base R, where `glm(f ~ x, gaussian)` on a factor `f`
signals *"NA/NaN/Inf in 'y'"* after three warnings from `Ops.factor`. A
two-level factor has one numeric reading and this is it; on a gaussian
family the result is the linear probability model for the response
`as.numeric(f == "b")`.

The levels used are returned as the attribute `"response_levels"`, which
lets
[`statmod_respec()`](https://statmodels7.github.io/statmodels7/reference/statmod_respec.md)
reproduce the same mapping on other rows. A factor carries its levels
through a subset, but a character vector does not, so a cross-validation
fold holding only one of the two values would otherwise be coded the
other way round, and the deviance would be wrong with no signal.

## See also

[`statmod_spec()`](https://statmodels7.github.io/statmodels7/reference/statmod_spec.md),
[`statmod_respec()`](https://statmodels7.github.io/statmodels7/reference/statmod_respec.md)
