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

A unit records the parameter its penalty is filed under, which for a
covariance class is the parameter the class is keyed by and not
necessarily the one carrying the structural term the class reaches into.
A caller that has a term's key and wants the term asks here rather than
assuming the two agree.
