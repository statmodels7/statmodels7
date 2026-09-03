# Split a Multi-Parameter Formula Into One Equation Per Parameter

Takes the single formula
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
accepts, in which the equations of the distribution's parameters are
separated by `|`, and returns the response together with one one-sided
formula per parameter of the family. Every parameter gets an equation:
one a caller did not write is filled in as `~ 1`, an intercept.

## Usage

``` r
statmod_equations(formula, params)
```

## Arguments

- formula:

  A two-sided formula, with the parameters' equations separated by `|`
  as above. Its environment is carried onto every equation returned, so
  a term's symbols resolve where the caller wrote them; a formula with a
  `NULL` environment gets
  [`baseenv()`](https://rdrr.io/r/base/environment.html).

- params:

  The distribution's parameter names in the family's own order, as
  `distrib@params` gives them. A character vector. The first element is
  the parameter the response's own equation belongs to.

## Value

A list of three:

- `response`:

  the unevaluated left-hand side, a language object, usually a symbol
  but any expression the caller wrote.

- `equations`:

  a named list of one-sided formulas, one per element of `params` and in
  that order, whatever order the caller wrote them in. A parameter with
  no equation of its own holds `~ 1`.

- `given`:

  the parameter names the formula supplied, in the order it supplied
  them. A subset of `params`, possibly of length one.

## The syntax

        y ~ <terms>  |  p2 ~ <terms>  |  p3 ~ <terms>

The first equation carries the response and models the family's first
parameter, and each `|` introduces the next. The parameters may be named
in any order and any of them may be omitted.

## Why the recovery is a walk

R's own precedence decides this. `~` binds looser than `|` and
associates to the left, so `y ~ a | p2 ~ b | p3 ~ c` parses as

        ((y ~ (a | p2)) ~ (b | p3)) ~ c

The right-hand side of the whole formula is `c` alone. Splitting it on
`|` returns one piece and drops two equations without an error.

What the associativity calls for is a walk down the left spine of the
nested `~` calls, collecting each level's right-hand side and then
reversing. The innermost left operand is the response. Each collected
piece is either `<terms> | <name of the next parameter>`, or, for the
last, `<terms>` alone.

## Bars inside a call survive

The walk descends only through `~` and through a `|` at the top of a
collected piece. A `|` anywhere inside a call is left as it stands, so
`random(1 | id)` and `gas(by = ~ random(1 | id))` reach the term
constructor whole.

## What is refused

Five conditions, each with an error naming the offending name:

- `formula` is not a formula.

- `formula` has no left-hand side, so there is no response.

- A `|` is not followed by an equation, as in `y ~ x | sigma`.

- A name is not a parameter of this distribution. The message lists the
  ones that are.

- A parameter is given two equations.

## See also

[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md),
which calls this on the formula it is given,
[`statmod_spec()`](https://statmodels7.github.io/statmodels7/reference/statmod_spec.md)
for what is built from the result.

## Examples

``` r
e <- statmod_equations(y ~ x1 + x2 | sigma ~ z, c("mu", "sigma"))
e$response
#> y
e$equations
#> $mu
#> ~x1 + x2
#> <environment: 0x55872a693628>
#> 
#> $sigma
#> ~z
#> <environment: 0x55872a693628>
#> 

# Every parameter gets an equation; the ones not written get an intercept.
f <- statmod_equations(y ~ x, c("mu", "sigma", "nu"))
vapply(f$equations, function(q) deparse(q[[2]]), character(1))
#>    mu sigma    nu 
#>   "x"   "1"   "1" 
f$given
#> [1] "mu"

# Three equations survive, where splitting the right-hand side on "|"
# would keep only the last.
g <- y ~ a | sigma ~ b | nu ~ cc
deparse(g[[3]])                       # the whole right-hand side: just "cc"
#> [1] "cc"
statmod_equations(g, c("mu", "sigma", "nu"))$given
#> [1] "mu"    "sigma" "nu"   

# A bar inside a call is not a separator.
h <- statmod_equations(y ~ random(1 | id) | sigma ~ z, c("mu", "sigma"))
deparse(h$equations$mu[[2]])
#> [1] "random(1 | id)"

# A name that is not a parameter is refused, and the message says which
# names are.
try(statmod_equations(y ~ x | zz ~ z, c("mu", "sigma")))
#> Error : 'zz' is not a parameter of this distribution.
#>   Its parameters are: mu, sigma.
```
