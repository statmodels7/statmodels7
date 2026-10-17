# A Term by Its Key in Any Equation

Returns the term that a key names, searching over every distribution
parameter.

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

The term, or `NULL` if no equation carries that key.

## Details

A unit records the parameter under which its penalty is filed. For a
covariance class that is the parameter by which the class is keyed,
which is not necessarily the parameter carrying the structural term into
which the class reaches. A caller that has a term's key and needs the
term should use this function instead of assuming that the two agree.
