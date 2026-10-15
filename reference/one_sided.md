# Build a One-Sided Formula From an Expression

Wraps a term expression as `~ expr` and attaches `env` to it, so that
the symbols in a term resolve where the caller wrote them. Used by
[`statmod_equations()`](https://statmodels7.github.io/statmodels7/reference/statmod_equations.md)
on each piece it collects, and on `quote(1)` for a parameter with no
equation.

## Usage

``` r
one_sided(expr, env)
```

## Arguments

- expr:

  A language object, the right-hand side to wrap: a symbol, a call, or
  the literal `1`. Substituted unevaluated, so a term call is not run
  here.

- env:

  The environment to attach. It is not checked, and `NULL` gives a
  formula with a `NULL` environment.

## Value

A one-sided formula of length 2, whose `[[2]]` is `expr` and whose
[`environment()`](https://rdrr.io/r/base/environment.html) is `env`.

## Details

The formula is built with `eval(bquote(~ .(expr)), envir = env)`, and
the environment is then assigned again.
[`eval()`](https://rdrr.io/r/base/eval.html) gives the formula the
environment in which it was evaluated, and the assignment states the
intended environment explicitly, so the result carries `env` whatever
the construction does.

## See also

[`statmod_equations()`](https://statmodels7.github.io/statmodels7/reference/statmod_equations.md),
its only caller.
