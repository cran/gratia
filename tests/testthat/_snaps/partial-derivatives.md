# partial derivatives works with factor by

    Code
      pd
    Output
      # A tibble: 200 x 11
         .smooth    .focal .by   .fs   .partial_deriv    .se .crit .lower_ci .upper_ci
         <chr>      <chr>  <chr> <chr>          <dbl>  <dbl> <dbl>     <dbl>     <dbl>
       1 te(x1,x2)~ x1     z1    <NA>         -0.0553 0.0589  1.96    -0.171    0.0601
       2 te(x1,x2)~ x1     z1    <NA>         -0.0553 0.0589  1.96    -0.171    0.0601
       3 te(x1,x2)~ x1     z1    <NA>         -0.0553 0.0589  1.96    -0.171    0.0601
       4 te(x1,x2)~ x1     z1    <NA>         -0.0553 0.0589  1.96    -0.171    0.0601
       5 te(x1,x2)~ x1     z1    <NA>         -0.0553 0.0588  1.96    -0.171    0.0600
       6 te(x1,x2)~ x1     z1    <NA>         -0.0552 0.0587  1.96    -0.170    0.0599
       7 te(x1,x2)~ x1     z1    <NA>         -0.0552 0.0586  1.96    -0.170    0.0597
       8 te(x1,x2)~ x1     z1    <NA>         -0.0551 0.0585  1.96    -0.170    0.0595
       9 te(x1,x2)~ x1     z1    <NA>         -0.0550 0.0583  1.96    -0.169    0.0591
      10 te(x1,x2)~ x1     z1    <NA>         -0.0549 0.0580  1.96    -0.169    0.0587
      # i 190 more rows
      # i 2 more variables: x1 <dbl>, z1 <fct>

# partial derivatives works with selected factor by

    Code
      pd
    Output
      # A tibble: 100 x 11
         .smooth    .focal .by   .fs   .partial_deriv    .se .crit .lower_ci .upper_ci
         <chr>      <chr>  <chr> <chr>          <dbl>  <dbl> <dbl>     <dbl>     <dbl>
       1 te(x1,x2)~ x1     z1    <NA>         -0.0553 0.0589  1.96    -0.171    0.0601
       2 te(x1,x2)~ x1     z1    <NA>         -0.0553 0.0589  1.96    -0.171    0.0601
       3 te(x1,x2)~ x1     z1    <NA>         -0.0553 0.0589  1.96    -0.171    0.0601
       4 te(x1,x2)~ x1     z1    <NA>         -0.0553 0.0589  1.96    -0.171    0.0601
       5 te(x1,x2)~ x1     z1    <NA>         -0.0553 0.0588  1.96    -0.171    0.0600
       6 te(x1,x2)~ x1     z1    <NA>         -0.0552 0.0587  1.96    -0.170    0.0599
       7 te(x1,x2)~ x1     z1    <NA>         -0.0552 0.0586  1.96    -0.170    0.0597
       8 te(x1,x2)~ x1     z1    <NA>         -0.0551 0.0585  1.96    -0.170    0.0595
       9 te(x1,x2)~ x1     z1    <NA>         -0.0550 0.0583  1.96    -0.169    0.0591
      10 te(x1,x2)~ x1     z1    <NA>         -0.0549 0.0580  1.96    -0.169    0.0587
      # i 90 more rows
      # i 2 more variables: x1 <dbl>, z1 <fct>

# partial derivatives throws error with incorrect focal

    Code
      partial_derivatives(m_partial_deriv, focal = rep("x1", 1L))
    Message
      ! Ignoring univariate smooths & those involving random effects.
    Condition
      Error in `partial_derivatives()`:
      ! `focal` must have same length as number of chosen smooths
      i The relevant smooths are: te(x1,x2):z1A and te(x1,x2):z1B
      i The supplied `focal` should be length: 2
      x Your supplied `focal` was length: 1

