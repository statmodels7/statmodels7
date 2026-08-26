# Evaluate a Formula's Terms With modelterms7 in Front

Returns a fresh environment holding every function modelterms7 exports,
whose parent is `env`. A term call evaluated there reaches modelterms7's
function whatever the caller has attached, and every other name the
caller wrote resolves as usual, one step further out.

## Usage

``` r
terms_first(env)
```

## Arguments

- env:

  The environment the formula carried, which becomes the parent of the
  result. Everything visible from `env` stays visible.

## Value

A new environment holding one binding per exported function of
modelterms7, currently 82 of them, with `env` as its parent.

## Details

The masking this removes is real and its symptom names nothing. mgcv
exports `s()` and `te()`, and segmented exports `seg()`. With either
package attached, a caller writing a
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
formula gets that package's function, whose value is not a model term;
the failure then surfaces inside
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) as
`invalid type (list) for variable 's(x)'`, which names neither the call
nor the mask.

Putting modelterms7 in front removes the ambiguity at the point where
the formula is read. A caller who wants the other package's function
writes `mgcv::s(x)`, which reaches past the shim; modelterms7's formula
interpreter then rejects it, naming the class returned and the package
that supplied it.

Only functions are copied. modelterms7 exports objects that are not
functions, S7 classes among them, and shadowing those would change what
a name means without any of the benefit.

## See also

[`statmod_equations()`](https://statmodels7.github.io/statmodels7/reference/statmod_equations.md)
for the formula this is used to interpret,
[`modelterms7::interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.html)
for the reading it feeds.
