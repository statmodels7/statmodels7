# Split a Multi-Parameter Formula Into One Equation Per Parameter

Takes the single formula that
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
accepts, in which the equations of the distribution's parameters are
separated by `|`, and returns the response together with one one-sided
formula per parameter of the family. Every parameter receives an
equation, and a parameter for which the caller did not write an equation
receives `~ 1`, an intercept.

## Usage

``` r
statmod_equations(formula, params)
```

## Arguments

- formula:

  A two-sided formula in which the equations of the parameters are
  separated by `|` as above. Its environment is carried onto every
  equation returned, so a term's symbols resolve where the caller wrote
  them; if the formula has a `NULL` environment,
  [`baseenv()`](https://rdrr.io/r/base/environment.html) is used.

- params:

  The distribution's parameter names in the family's own order, as
  `distrib@params` gives them. It is a character vector. The first
  element is the parameter to which the response's own equation belongs.

## Value

A list of three:

- `response`:

  the unevaluated left-hand side, a language object, usually a symbol
  but any expression that the caller wrote.

- `equations`:

  a named list of one-sided formulas, one per element of `params` and in
  that order, regardless of the order in which the caller wrote them. A
  parameter with no equation of its own holds `~ 1`.

- `given`:

  the parameter names that the formula supplied, in the order in which
  it supplied them. It is a subset of `params`, possibly of length one.

## The syntax

        y ~ <terms>  |  p2 ~ <terms>  |  p3 ~ <terms>

The first equation carries the response and models the family's first
parameter, and each `|` introduces the next. The parameters may be named
in any order, and any of them may be omitted.

## Splitting the formula

The split follows R's own precedence. `~` binds looser than `|` and
associates to the left, so `y ~ a | p2 ~ b | p3 ~ c` parses as

        ((y ~ (a | p2)) ~ (b | p3)) ~ c

The right-hand side of the whole formula is `c` alone. Splitting it on
`|` returns one piece and drops two equations without an error.

The function therefore walks down the left spine of the nested `~`
calls, collects the right-hand side of each level and reverses the
collection. The innermost left operand is the response. Each collected
piece is either `<terms> | <name of the next parameter>` or, for the
last piece, `<terms>` alone.

## Bars inside a call

The walk descends only through `~` and through a `|` at the top of a
collected piece. A `|` anywhere inside a call is left as it stands, so
`random(~ 1 | id)` and `gas(omega ~ random(~ 1 | id), by = id)` reach
the term constructor whole.

## Invalid formulas

Five conditions signal an error:

- `formula` is not a formula.

- `formula` has no left-hand side, so there is no response.

- A `|` is not followed by an equation, as in `y ~ x | sigma`. The
  message names the parameter whose equation is missing.

- A name is not a parameter of this distribution. The message names it
  and lists the parameters of the distribution.

- A parameter is given two equations. The message names the parameter.

## See also

[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md),
which calls this on the formula that it is given,
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
#> <environment: 0x55fddf2317a0>
#> 
#> $sigma
#> ~z
#> <environment: 0x55fddf2317a0>
#> 

# Every parameter receives an equation; a parameter without one receives an
# intercept.
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
h <- statmod_equations(y ~ random(~ 1 | id) | sigma ~ z, c("mu", "sigma"))
deparse(h$equations$mu[[2]])
#> [1] "random(~1 | id)"

# A name that is not a parameter is an error, and the message lists the
# parameters of the distribution.
try(statmod_equations(y ~ x | zz ~ z, c("mu", "sigma")))
#> Error : 'zz' is not a parameter of this distribution.
#>   Its parameters are: mu, sigma.
```
