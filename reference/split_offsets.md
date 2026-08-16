# Take the Offsets Out of an Equation

Splits a one-sided formula into the
[`offset()`](https://rdrr.io/r/stats/offset.html) terms it names and
whatever is left, which is what the interpreter is given.

## Usage

``` r
split_offsets(eq)
```

## Arguments

- eq:

  A one-sided formula.

## Value

A list with `formula`, the equation with the offsets removed, and
`offsets`, a list of the expressions taken out.

## Details

An offset is a column of the linear predictor whose coefficient is known
to be one, and `y ~ x + offset(log_n)` is how R has always written it.
Without this the term reached `model.matrix` through
[`linpar`](https://statmodels7.github.io/modelterms7/reference/linpar.html),
where [`terms()`](https://rdrr.io/r/stats/terms.html) marks it in the
`"offset"` attribute and the design EXCLUDES it: the term contributed no
column, no offset and no message, and the model fitted was the one
without it. On a count model over person-years that moved the intercept
from -7.5 to -0.6, which is the difference between a log rate and a log
count.

Only a top-level additive term is taken, which is where R recognizes one
too; `stats::offset(x)` is recognized beside `offset(x)`. Several are
summed, as [`glm`](https://rdrr.io/r/stats/glm.html) sums them.

## See also

[`statmod_spec`](https://statmodels7.github.io/statmodels7/reference/statmod_spec.md),
[`eval_offsets`](https://statmodels7.github.io/statmodels7/reference/eval_offsets.md)
