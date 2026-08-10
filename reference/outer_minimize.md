# Is a Criterion Minimized?

`TRUE` for a prediction-error criterion, `FALSE` for a marginal
likelihood, which is maximized.

## Usage

``` r
outer_minimize(method)
```

## Arguments

- method:

  An
  [`OuterMethod`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md).

## Value

A single logical.

## Details

The search always minimizes; this is what says whether the criterion
goes in as it is or with its sign turned. It is asked of the method
rather than written into the search, so a criterion added later declares
its own orientation.
