# Bring an Entry to a Number of Replicates

Looks the entry for `resample` and `refit` up in the collection of the
fit, computes the replicates that are missing for `n` and stores the
result.

## Usage

``` r
boot_fill(fit, n, resample, refit, say = FALSE, cores = 1L)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- n:

  The total number of replicates wanted.

- resample, refit:

  The settings.

- say:

  Whether to report the replicates computed with a message.

## Value

The entry.

## Details

A new entry starts the stream at the current state of `.Random.seed`. An
existing entry continues from the state stored after its last replicate.
The write to the environment of the fit is made here, in the calling
process, after the engine has collected the refits.
