# A Sentence for Each Direction That separation_diagnosis() Found

Describes each direction of recession and each direction of descent in
the result of
[`separation_diagnosis()`](https://statmodels7.github.io/statmodels7/reference/separation_diagnosis.md)
in one sentence, naming the coefficients that move and the number of
observations whose linear predictor moves with them.

## Usage

``` r
separation_message(sep)
```

## Arguments

- sep:

  The list that
  [`separation_diagnosis()`](https://statmodels7.github.io/statmodels7/reference/separation_diagnosis.md)
  returns.

## Value

A character vector with one sentence per direction, empty where there is
none.
