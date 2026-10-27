# Checks of the Arguments of a Diagnostic Graph

Signal an error unless the argument has the stated form.

## Usage

``` r
diag_check_seed(seed)

diag_check_count(x, name, min = 1, null = FALSE)

diag_check_name(x, name, null = TRUE)

diag_check_level(level)
```

## Arguments

- seed:

  A single number or `NULL`.

- x:

  The value to check.

- name:

  The name of the argument in the message.

- min:

  The smallest admissible whole number.

- null:

  Whether `NULL` is admissible.

## Value

`NULL`, invisibly.
