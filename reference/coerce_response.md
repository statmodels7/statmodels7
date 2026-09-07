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

A binary outcome recorded as a factor is how these data are ordinarily
written, and [`glm()`](https://rdrr.io/r/stats/glm.html) has always
accepted one: its binomial family reads the first level as the failure.
Without this the response reached
[`stats::dbinom()`](https://rdrr.io/r/stats/Binomial.html) untouched and
the run died on *"Non-numeric argument to mathematical function"*, which
names neither the variable nor the cause.

The first level is the failure, which is
[`glm()`](https://rdrr.io/r/stats/glm.html)'s rule, so a model moved
from [`glm()`](https://rdrr.io/r/stats/glm.html) keeps its signs. The
conversion is made for EVERY family and not only for the binary ones,
which is where this departs from base R – and it departs in the
direction of working: `glm(f ~ x, gaussian)` on a factor `f` does not
refuse cleanly, it signals *"NA/NaN/Inf in 'y'"* after three warnings
from `Ops.factor`. A two-level factor has one numeric reading and this
is it; on a gaussian family the result is the linear probability model,
which is what `as.numeric(f == "b")` would have given.

The levels used are returned as the attribute `"response_levels"`. That
is what lets
[`statmod_respec()`](https://statmodels7.github.io/statmodels7/reference/statmod_respec.md)
reproduce the SAME mapping on other rows: a factor carries its levels
through a subset, but a character vector does not, so a cross-validation
fold holding only one of the two values would otherwise be coded the
other way round and the deviance would come back wrong without a word.

## See also

[`statmod_spec()`](https://statmodels7.github.io/statmodels7/reference/statmod_spec.md),
[`statmod_respec()`](https://statmodels7.github.io/statmodels7/reference/statmod_respec.md)
