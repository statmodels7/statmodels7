# A Term by Its Key, Whichever Equation It Sits In

The term a key names, searched over every distribution parameter.

## Usage

``` r
find_term(spec, key)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- key:

  A term's key, its call as written.

## Value

The term, or `NULL` where no equation carries that key.

## Details

A unit records the parameter under which its penalty is filed. For a
covariance class that is the parameter by which the class is keyed,
which is not necessarily the parameter carrying the structural term that
the class reaches into. A caller that has a term's key and needs the
term should use this function instead of assuming that the two agree.
