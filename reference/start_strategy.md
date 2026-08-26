# S7 Class for a Starting-Value Strategy

The abstract parent of the starting-value strategies. A strategy says
where a fit begins as a **procedure**, not as a vector of numbers: fit
the intercept-only model, draw around it, or search the likelihood for a
basin. `statmod(start =)` accepts one, or a plain named list of
coefficients.

## Usage

``` r
start_strategy(label = character(0))
```

## Arguments

- label:

  A short name, printed by
  [`print.start_strategy()`](https://statmodels7.github.io/statmodels7/reference/print.start_strategy.md).
  A single string. Set by each subclass's constructor, since the class
  itself cannot be built.

## Value

Nothing: the class is abstract and the constructor signals an error.
Used as a class object, for
[`S7::S7_inherits()`](https://rconsortium.github.io/S7/reference/S7_inherits.html)
and for registering methods, it has one property, `label`, which every
subclass inherits.

## Details

A strategy is asked exactly once, before the alternation between the
coefficients and the hyperparameters begins. An optimizer runs at every
step of that alternation, which is why a global search belongs here:
given to `inner_optimizer` it would rerun at every hyperparameter the
outer search tried and return the same answer each time.

**The class is abstract and cannot be instantiated.** Calling
`start_strategy()` signals an error; what the class is for is to be
subclassed and to make a membership test possible. To write a strategy
of your own, subclass it and register a
[`start_at()`](https://statmodels7.github.io/statmodels7/reference/start_at.md)
method on the subclass. That method receives the specification, the
design and the objective, and returns one vector per distribution
parameter.

## See also

[`start_intercepts()`](https://statmodels7.github.io/statmodels7/reference/start_intercepts.md)
(the default),
[`start_origin()`](https://statmodels7.github.io/statmodels7/reference/start_origin.md),
[`start_random()`](https://statmodels7.github.io/statmodels7/reference/start_random.md)
and
[`start_search()`](https://statmodels7.github.io/statmodels7/reference/start_search.md)
for the four shipped strategies,
[`start_at()`](https://statmodels7.github.io/statmodels7/reference/start_at.md)
for the generic they implement.

## Examples

``` r
# The four shipped strategies all inherit from this.
S7::S7_inherits(start_origin(), start_strategy)
#> [1] TRUE
S7::S7_inherits(start_search(), start_strategy)
#> [1] TRUE

# Every one of them carries the label the class defines.
vapply(list(start_intercepts(), start_origin(), start_random(),
            start_search()),
       function(s) s@label, character(1))
#> [1] "intercept-only fit"                       
#> [2] "zero"                                     
#> [3] "random around the intercept-only fit"     
#> [4] "search with simulated annealing (uniform)"

# The class itself is abstract.
try(start_strategy("my own"))
#> Error in S7::new_object(S7::S7_object(), label = label) : 
#>   Can't construct an object from abstract class <start_strategy>
```
