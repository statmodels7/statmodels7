# Check the Signatures of a Family's Own Methods

Signals an error naming the method and the signature it needs where a
family registers a method of its own, for the density or for a
derivative, whose formals have no `...`.

## Usage

``` r
assert_method_signatures(distrib)
```

## Arguments

- distrib:

  A univariate family.

## Value

`distrib`, invisibly.

## Details

[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
passes further arguments to these methods, the number of threads among
them, and a method without `...` stopped with "unused argument (threads
= 1)" from inside the fit. The derivative methods also take `scale`,
which the generic resolves and passes on. The methods of the families
shipped with distributions7 all carry `...`; a family written by a user
is where this is met, and where the message is read.
