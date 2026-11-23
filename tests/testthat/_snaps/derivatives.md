# derivatives works for sz smooths

    Code
      print(d)
    Output
      # A tibble: 500 x 11
         .smooth .by   .fs   .derivative   .se .crit .lower_ci .upper_ci      x2 fac  
         <chr>   <chr> <chr>       <dbl> <dbl> <dbl>     <dbl>     <dbl>   <dbl> <fct>
       1 s(x2)   <NA>  <NA>         21.5  5.56  1.96      10.6      32.4 0.00131 <NA> 
       2 s(x2)   <NA>  <NA>         21.5  5.55  1.96      10.7      32.4 0.0114  <NA> 
       3 s(x2)   <NA>  <NA>         21.5  5.49  1.96      10.8      32.3 0.0215  <NA> 
       4 s(x2)   <NA>  <NA>         21.5  5.35  1.96      11.0      31.9 0.0316  <NA> 
       5 s(x2)   <NA>  <NA>         21.4  5.11  1.96      11.3      31.4 0.0417  <NA> 
       6 s(x2)   <NA>  <NA>         21.2  4.80  1.96      11.8      30.6 0.0517  <NA> 
       7 s(x2)   <NA>  <NA>         20.9  4.42  1.96      12.3      29.6 0.0618  <NA> 
       8 s(x2)   <NA>  <NA>         20.6  4.00  1.96      12.7      28.4 0.0719  <NA> 
       9 s(x2)   <NA>  <NA>         20.1  3.59  1.96      13.1      27.2 0.0820  <NA> 
      10 s(x2)   <NA>  <NA>         19.5  3.25  1.96      13.2      25.9 0.0921  <NA> 
      # i 490 more rows
      # i 1 more variable: x0 <dbl>

# derivatives gives nice alerts for multivariate smooths

    Code
      derivatives(gam(y ~ s(x0) + s(x1, x2), data = d_dat))
    Message
      i Can't compute derivatives of multivariate smooths
      i Ignoring: "s(x1,x2)"
      i See `partial_derivatives()` for a solution.
    Output
      # A tibble: 100 x 9
         .smooth .by   .fs   .derivative   .se .crit .lower_ci .upper_ci      x0
         <chr>   <chr> <chr>       <dbl> <dbl> <dbl>     <dbl>     <dbl>   <dbl>
       1 s(x0)   <NA>  <NA>         6.72  2.71  1.96      1.41      12.0 0.00206
       2 s(x0)   <NA>  <NA>         6.72  2.71  1.96      1.42      12.0 0.0121 
       3 s(x0)   <NA>  <NA>         6.71  2.68  1.96      1.45      12.0 0.0222 
       4 s(x0)   <NA>  <NA>         6.69  2.64  1.96      1.52      11.9 0.0322 
       5 s(x0)   <NA>  <NA>         6.66  2.57  1.96      1.62      11.7 0.0423 
       6 s(x0)   <NA>  <NA>         6.61  2.48  1.96      1.75      11.5 0.0524 
       7 s(x0)   <NA>  <NA>         6.55  2.38  1.96      1.89      11.2 0.0624 
       8 s(x0)   <NA>  <NA>         6.47  2.26  1.96      2.03      10.9 0.0725 
       9 s(x0)   <NA>  <NA>         6.37  2.15  1.96      2.16      10.6 0.0826 
      10 s(x0)   <NA>  <NA>         6.26  2.04  1.96      2.25      10.3 0.0926 
      # i 90 more rows

---

    Code
      derivatives(gam(y ~ s(x1, x2), data = d_dat))
    Condition
      Error in `derivatives()`:
      ! Can't compute derivatives for any smooths in `gam(y ~ s(x1, x2), data = d_dat)`.
      x All smooths are either random effects or multivariate
      ! See `partial_derivatives()` for one solution for multivariate smooths.

