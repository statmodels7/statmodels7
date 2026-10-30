# The Positions of the Break-Points of a Model

Lists every position of a break-point of a
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html),
[`modelterms7::jump()`](https://statmodels7.github.io/modelterms7/reference/jump.html)
or
[`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html)
term that is a single coefficient, that is a position without a
subformula.

## Usage

``` r
profile_breakpoints(spec, design)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  Its design.

## Value

A data frame with one row for each position and the columns `label` (the
label in `which`), `slot_label` (the label of the coefficient in the
slot), `param`, `term`, `k` (the break-point), `at` (the stacked
position of the slot), `sharp` (a sharp jump or jseg), `sharp_seg` (a
sharp continuous seg), `estimate`, `lower` and `upper` (the interval in
which the term confines the position), and the list column `covariate`,
the values of the covariate of the term.

## Details

The position of a continuous or smoothed term is a coefficient of the
design, and its label is that of the coefficient. The position of a
sharp
[`modelterms7::jump()`](https://statmodels7.github.io/modelterms7/reference/jump.html)
or
[`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html)
is not: the coefficient in its slot is the auxiliary \\g_k =
-\delta_k\psi_k\\ of the working linearization, labelled `g1`, `g2`, and
so on. Its label in `which` is then the one that
[`coef()`](https://rdrr.io/r/stats/coef.html) prints, with `psi1`,
`psi2` in place of `g1`, `g2`, and
[`statmod_profile()`](https://statmodels7.github.io/statmodels7/reference/statmod_profile.md)
holds it by moving the term to the value and holding it there
([`modelterms7::seg_relocate()`](https://statmodels7.github.io/modelterms7/reference/seg_relocate.html),
[`modelterms7::seg_hold()`](https://statmodels7.github.io/modelterms7/reference/seg_hold.html)).
