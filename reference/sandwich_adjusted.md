# Rescale the Scores of a Fit by the Leverage

Computes the scores \\s_i^\*\\ of
[`hc2_sandwich()`](https://statmodels7.github.io/statmodels7/reference/hc2_sandwich.md)
(or of the CR2, CR3 and HC4 variants) and, for CR2, the closure that
evaluates the Satterthwaite degrees of freedom.

## Usage

``` r
sandwich_adjusted(
  fit,
  design,
  Vb,
  cols,
  g,
  id,
  sw,
  expected,
  emethod,
  lev_info = NULL
)
```

## Arguments

- fit:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- design:

  The design of the fit.

- Vb:

  The inverse penalized information over the coordinates `cols`.

- cols:

  The positions, among the coefficients stacked by equation, of the
  coordinates that `Vb` covers.

- g:

  The list of weighted scores of
  [`statmod_score_obs()`](https://statmodels7.github.io/statmodels7/reference/statmod_score_obs.md).

- id:

  `NULL` or the cluster numbers.

- sw:

  The
  [SandwichMethod](https://statmodels7.github.io/statmodels7/reference/SandwichMethod.md)
  object.

- expected, emethod:

  The information that the bread was built from.

- lev_info:

  The information of the leverage: the roots \\W_i^{1/2}\\ and the
  inverse \\K^{-1}\\ of the hat matrix are read on it. `NULL` takes the
  information of the bread.

## Value

A list with `g`, the list `g` with the rescaled scores in the place of
the originals, and `df_fn`, `NULL` unless the type is CR2 with a `t`
reference, otherwise a function of the positions (among `cols`) that
returns the degrees of freedom there.
