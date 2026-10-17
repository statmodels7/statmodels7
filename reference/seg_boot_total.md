# Restart Budget of the Break-Point Terms

Returns the largest `n_boot` set in a break-point term that is written
directly in an equation of the specification, and zero when there is no
such term. A break-point term nested inside a parameter of
[`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html)
is not counted.
[`statmod_boot_restart()`](https://statmodels7.github.io/statmodels7/reference/statmod_boot_restart.md)
takes this number as the maximum number of proposals that it tries.

## Usage

``` r
seg_boot_total(spec)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).
  The terms of every equation are examined.

## Value

A single non-negative integer. It is zero when no term sets `n_boot`,
and
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
does not run the restart loop for a budget of zero.

## Details

The number is set in the term.
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html),
[`modelterms7::jump()`](https://statmodels7.github.io/modelterms7/reference/jump.html)
and
[`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html)
each take `n_boot`, with a default of 10, because the objective of these
terms has spurious local optima, and the restarts search for a better
one. The restarts are run by this package, which can refit the model.

A model with two break-point terms that set `n_boot` to 10 and to 25 has
a budget of 25. The restart loop works on all of the terms together,
since a proposal is a configuration of positions for every break-point
term in the model, so the budget is a single number and the largest of
the set values is used.

## See also

[`statmod_boot_restart()`](https://statmodels7.github.io/statmodels7/reference/statmod_boot_restart.md),
which spends this budget,
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html)
for where the number is set.
