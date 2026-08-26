# How Many Restarts the Terms Ask For

Returns the largest `n_boot` any break-point term of the specification
declares, and zero when the model carries no such term. This is the
budget
[`statmod_boot_restart()`](https://statmodels7.github.io/statmodels7/reference/statmod_boot_restart.md)
spends: how many proposals it may try before giving up on improving the
fit.

## Usage

``` r
seg_boot_total(spec)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md),
  whose terms are walked in every equation.

## Value

A single non-negative integer. Zero when no term declares `n_boot`, and
a budget of zero turns the restart loop off.

## Details

The number is declared on the term.
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html),
[`modelterms7::jump()`](https://statmodels7.github.io/modelterms7/reference/jump.html)
and
[`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html)
each take `n_boot`, with a default of 10, because those are the terms
whose objective has the spurious local optima the device exists for.
Running the restarts belongs here, in the layer that can refit the
model, and the split follows the one a penalty already uses: the term
says what it needs, this package does it.

A model with two break-point terms asking for 10 and 25 gets 25. The
restart loop works on all of them together, since a proposal moves every
break-point in the model at once, so the budget is one number and the
largest request is the one honored.

## See also

[`statmod_boot_restart()`](https://statmodels7.github.io/statmodels7/reference/statmod_boot_restart.md),
which spends this budget,
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html)
for where the number is set.
