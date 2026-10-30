# Run the Gradient Check of optimizers7 Once per Site and Fit

Evaluates `expr`, a call of
[`optimizers7::minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.html),
with the check of the supplied gradient switched off (the option
`optimizers7.check_gradient`) when the same `site` has already been
checked within the current call of
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md).

## Usage

``` r
with_gradient_check(site, expr)
```

## Arguments

- site:

  A single string naming the routine and its block.

- expr:

  The call to evaluate.

## Value

The value of `expr`.

## Details

`minimize()` compares the gradient with a central difference of the
objective at the start of every call. Within one fit the objective and
the gradient at a site are the same code at every pass and every start,
so the first comparison is the only informative one. A site names the
block as well as the routine, since two blocks are two objectives.
Outside a call of
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
every check runs.
