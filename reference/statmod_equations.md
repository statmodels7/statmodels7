# Split a Multi-Parameter Formula Into One Equation Per Parameter

Takes the single formula
[`statmod`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
accepts, in which the equations of the distribution's parameters are
separated by `|`, and returns the response together with one one-sided
formula per parameter.

## Usage

``` r
statmod_equations(formula, params)
```

## Arguments

- formula:

  The model formula.

- params:

  The distribution's parameter names, in the family's order.

## Value

A list with `response` (the unevaluated left-hand side), `equations` (a
named list of one-sided formulas, one per element of `params`, in that
order) and `given` (the names the formula actually supplied).

## Details

The syntax is

        y ~ <terms>  |  p2 ~ <terms>  |  p3 ~ <terms>

with the first equation carrying the response and modelling the family's
first parameter, and each `|` introducing the next.

The recovery is not the obvious one, and the reason is R's own
precedence. `~` binds looser than `|` and associates to the left, so
`y ~ a | p2 ~ b | p3 ~ c` is the tree `((y ~ (a | p2)) ~ (b | p3)) ~ c`:
the right-hand side of the whole formula is `c` alone, and splitting it
on `|` returns one piece and silently drops two equations. What the
associativity calls for instead is a walk down the left spine of the
nested `~` calls, collecting each level's right-hand side; the innermost
left operand is the response, and each collected piece is either
`<terms> | <name of the next parameter>` or, for the last, `<terms>`
alone.

A `|` inside a call is untouched, the walk descending only through `~`
and the top-level `|`, so `random(1 | id)` and
`gas(by = ~ random(1 | id))` survive intact.

## See also

[`statmod`](https://statmodels7.github.io/statmodels7/reference/statmod.md)

## Examples

``` r
statmod_equations(y ~ x1 + x2 | sigma ~ z, c("mu", "sigma"))
#> $response
#> y
#> 
#> $equations
#> $equations$mu
#> ~x1 + x2
#> <environment: 0x55c80ab23ec8>
#> 
#> $equations$sigma
#> ~z
#> <environment: 0x55c80ab23ec8>
#> 
#> 
#> $given
#> [1] "mu"    "sigma"
#> 
```
