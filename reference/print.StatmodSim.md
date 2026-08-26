# Printing a Simulation

The model that was written, the size of what came out, and the truth
behind it.

## Usage

``` r
# S3 method for class 'StatmodSim'
print(x, ...)
```

## Arguments

- x:

  A
  [`rstatmod()`](https://statmodels7.github.io/statmodels7/reference/rstatmod.md)
  result.

- ...:

  Unused.

## Value

`x`, invisibly.

## See also

[`rstatmod()`](https://statmodels7.github.io/statmodels7/reference/rstatmod.md)

## Examples

``` r
set.seed(1)
rstatmod(y ~ x, distributions7::gaussian1_distrib(),
         data.frame(x = runif(20)))
#> A simulation from a written model
#> 
#>   20 observations, 2 columns
#>   mu          1.3587  -0.1028
#>   sigma      0.3877
#> 
#>   the data is in $data, the truth in $par and $theta
```
