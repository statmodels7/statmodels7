# How Many Restarts the Terms Ask For

The largest `n_boot` any break-point term of the specification declares,
or zero. The value is declared on
[`seg`](https://statmodels7.github.io/modelterms7/reference/seg.html),
[`jump`](https://statmodels7.github.io/modelterms7/reference/jump.html)
and
[`jseg`](https://statmodels7.github.io/modelterms7/reference/jseg.html),
the terms whose objective has the spurious local optima the device
exists for; running the restarts is this layer's, the way a penalty's
hyperparameters are declared by a term and estimated here.

## Usage

``` r
seg_boot_total(spec)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

## Value

A single non-negative integer.

## See also

[`statmod_boot_restart`](https://statmodels7.github.io/statmodels7/reference/statmod_boot_restart.md)
