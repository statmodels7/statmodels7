# Evaluate a Formula's Terms With modelterms7 in Front

Returns a fresh environment holding every function that modelterms7
exports, whose parent is `env`. A term call evaluated there reaches
modelterms7's function whatever the caller has attached, and every other
name that the caller wrote resolves as usual, one step further out.

## Usage

``` r
terms_first(env)
```

## Arguments

- env:

  The environment that the formula carried, which becomes the parent of
  the result. Everything visible from `env` stays visible.

## Value

A new environment holding one binding per exported function of
modelterms7, with `env` as its parent.

## Details

mgcv exports `s()` and `te()`, and segmented exports `seg()`. If either
package is attached, a
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
formula that uses one of these names would call that package's function,
whose value is not a model term. The failure then occurs inside
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) as
`invalid type (list) for variable 's(x)'`, a message that does not
mention the masking package.

Putting modelterms7 in front removes the ambiguity at the point where
the formula is read. A caller who wants the other package's function
writes `mgcv::s(x)`, which reaches past the shim; the formula
interpreter of modelterms7 then rejects the value, naming the class
returned and the call as written.

The function binds each export of modelterms7 for which
[`is.function()`](https://rdrr.io/r/base/is.function.html) is `TRUE`.
The S7 class generators of that package are functions and are copied
with the others.

## See also

[`statmod_equations()`](https://statmodels7.github.io/statmodels7/reference/statmod_equations.md)
for the formula that is interpreted with this function,
[`modelterms7::interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.html)
for the reading that it feeds.
