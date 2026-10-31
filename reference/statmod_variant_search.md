# A Fit Searched over the Versions of a Term

Fits the model once with each term that carries a choice replaced by its
pilot
([`modelterms7::term_pilot()`](https://statmodels7.github.io/modelterms7/reference/term_pilot.html)),
then, one such term at a time, once for each of its versions
([`modelterms7::term_variants()`](https://statmodels7.github.io/modelterms7/reference/term_variants.html)),
and returns the fit whose score
([`variant_score()`](https://statmodels7.github.io/statmodels7/reference/variant_score.md))
is the smallest.

## Usage

``` r
statmod_variants

statmod_variant_search(cl, caller, auto)
```

## Format

`statmod_variants` is an environment holding, during a search, whether
the search is running (`active`) and which terms replace which in the
formula (`override`, a list by parameter of terms keyed by the call as
written, which
[`statmod_terms()`](https://statmodels7.github.io/statmodels7/reference/statmod_terms.md)
reads); it is empty outside a search.

## Arguments

- cl:

  The call of
  [`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md).

- caller:

  The environment in which it was made.

- auto:

  The terms that carry a choice, as
  [`variant_terms()`](https://statmodels7.github.io/statmodels7/reference/variant_terms.md)
  returns them.

## Value

A
[StatmodFit](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

## Details

Every fit is the call of
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
as written, evaluated where it was made, with the replacement recorded
in statmod_variants; the search is not entered again from inside. Where
several terms carry a choice, each is settled in turn with the earlier
ones at their chosen version and the later ones at their pilot. The
versions tried, with their scores, are in `fit@history$variants`, and
the fit's `elapsed` covers the whole search.
