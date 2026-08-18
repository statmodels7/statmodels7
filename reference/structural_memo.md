# Reuse a Structural Quantity Computed at the Same Point

A depth-one exact memo on the design's structural state: where the last
call's key is [`identical()`](https://rdrr.io/r/base/identical.html) to
this one's, the stored value is returned; otherwise `compute()` runs and
replaces it.

## Usage

``` r
structural_memo(design, slot, key, compute)
```

## Arguments

- design:

  The design whose structural state holds the slots.

- slot:

  A short name, one cache per slot.

- key:

  The exact inputs the value depends on.

- compute:

  A function of no arguments.

## Value

`compute()`'s value, possibly from the cache.

## Details

Measured on the gas panel at 60 groups, 62 of the 154 curvature
recursions of one fit recompute a point already visited – the same
coefficients and the same term parameters, up to five times each,
because the criterion, its gradient and the joint step's curvature all
read the same mode – and the recursion is 35 per cent of the fit. The
cache returns the previously computed object itself, so a hit is
bit-identical to recomputing by construction, and the key is compared
with [`identical()`](https://rdrr.io/r/base/identical.html) on the full
numeric inputs, so a collision cannot happen.

It stands aside where the design carries REFRESHABLE terms: a
break-point block advances its rescaling schedule as the alternation
commits, so the same coefficients do not imply the same design there,
and a key that cannot see the schedule must not answer.
