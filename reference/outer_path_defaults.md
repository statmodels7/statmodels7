# The Properties Every Criterion Carries

The ones
[`OuterMethod`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md)
needs whether or not a given criterion uses them, so that one class
serves every criterion.

## Usage

``` r
outer_path_defaults()
```

## Value

A named list.

## Details

What a PATH does is not among them. How many values it visits, how far
down it reaches and whether a term's own hyperparameters are combined or
swept one at a time all belong to the term, since the same criterion is
put to the smooth hyperparameters of the model as well and those are
read at the mode rather than swept.

## See also

[`path_fallbacks`](https://statmodels7.github.io/statmodels7/reference/path_fallbacks.md)
