# Store the Position of a Held Break-Point in Its Slot

A sharp
[`modelterms7::jump()`](https://statmodels7.github.io/modelterms7/reference/jump.html)
or
[`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html)
is settled on the profile and held: its slot in the coefficient vector
is a working column whose value is the same in every replicate, and the
position of the break-point is kept in the term. The function writes the
position of each held break-point into its slot, so that the replicates
carry it.

## Usage

``` r
boot_record_positions(spec_b, design, cf)
```

## Arguments

- spec_b:

  The specification of the replicate.

- design:

  The design that the refit ran on, whose state holds the terms as the
  refit left them.

- cf:

  The coefficients of the refit.

## Value

The coefficients, with the held slots replaced.
