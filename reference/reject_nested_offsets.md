# Reject an Offset Buried Inside a Term

Stops where [`offset()`](https://rdrr.io/r/stats/offset.html) appears
anywhere other than as a term of the equation itself.

## Usage

``` r
reject_nested_offsets(eq, param)
```

## Arguments

- eq:

  A one-sided formula, with the equation's own offsets already taken out
  by
  [`split_offsets()`](https://statmodels7.github.io/statmodels7/reference/split_offsets.md).

- param:

  The parameter the equation belongs to, named in the message.

## Value

`NULL`, invisibly, when nothing is found. Signals an error naming
`param` and the offending call otherwise.

## Details

[`split_offsets()`](https://statmodels7.github.io/statmodels7/reference/split_offsets.md)
takes the top-level additive terms, which is where R recognizes an
offset. One written inside another term's formula, as in
`ridge(~ z + offset(o))`, `random(~ 1 + offset(o) | g)` or
`nl(a ~ 0 + ridge(~ g + offset(o)))`, reaches that term's own
`model.matrix` and is dropped there exactly as it used to be dropped at
the equation level. Measured: the fit ran, the block carried the columns
of the model without it, and the intercept came back at 1.33 against the
-5.01 the offset gives, a factor of 566.

Such a term is **refused**, never routed up to the equation, because the
meaning differs by where it sits. In a penalized term's formula an
offset would be a contribution to the equation's predictor, and writing
it at the equation level already says that. In a subformula it would be
a contribution to that parameter's own chart, a different quantity.
Picking one reading would fit the wrong model in half the cases; the
refusal names the place it belongs and costs a caller one edit.

The whole expression tree is walked, so a term written later is covered.

## See also

[`split_offsets()`](https://statmodels7.github.io/statmodels7/reference/split_offsets.md)
