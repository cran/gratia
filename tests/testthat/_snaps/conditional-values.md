# conditional_values works with vector condition

    Code
      print(cv)
    Output
      # A tibble: 500 x 9
          .row       x2       x1    x0    x3 .fitted   .se .lower_ci .upper_ci
         <int>    <dbl>    <dbl> <dbl> <dbl>   <dbl> <dbl>     <dbl>     <dbl>
       1     1 0.000571 0.000605 0.482 0.492    2.51 0.468      1.60      3.43
       2     2 0.000571 0.228    0.482 0.492    3.27 0.424      2.44      4.10
       3     3 0.000571 0.481    0.482 0.492    4.22 0.420      3.40      5.05
       4     4 0.000571 0.757    0.482 0.492    6.15 0.418      5.33      6.97
       5     5 0.000571 0.999    0.482 0.492    8.67 0.465      7.76      9.59
       6     6 0.0107   0.000605 0.482 0.492    2.88 0.424      2.05      3.71
       7     7 0.0107   0.228    0.482 0.492    3.64 0.375      2.90      4.37
       8     8 0.0107   0.481    0.482 0.492    4.59 0.371      3.86      5.32
       9     9 0.0107   0.757    0.482 0.492    6.52 0.369      5.80      7.24
      10    10 0.0107   0.999    0.482 0.492    9.04 0.422      8.22      9.87
      # i 490 more rows

# conditional_values works with complex list condition

    Code
      print(cv)
    Output
      # A tibble: 4,500 x 9
          .row       x2       x1    x0    x3 .fitted   .se .lower_ci .upper_ci
         <int>    <dbl>    <dbl> <dbl> <dbl>   <dbl> <dbl>     <dbl>     <dbl>
       1     1 0.000571 0.000605 0.258 0.197    1.98 0.479     1.04       2.92
       2     2 0.000571 0.000605 0.258 0.499    2.03 0.474     1.10       2.96
       3     3 0.000571 0.000605 0.258 0.800    2.08 0.477     1.15       3.02
       4     4 0.000571 0.000605 0.483 0.197    2.46 0.475     1.53       3.39
       5     5 0.000571 0.000605 0.483 0.499    2.51 0.468     1.60       3.43
       6     6 0.000571 0.000605 0.483 0.800    2.57 0.471     1.64       3.49
       7     7 0.000571 0.000605 0.747 0.197    1.91 0.485     0.963      2.87
       8     8 0.000571 0.000605 0.747 0.499    1.97 0.478     1.03       2.90
       9     9 0.000571 0.000605 0.747 0.800    2.02 0.479     1.08       2.96
      10    10 0.000571 0.228    0.258 0.197    2.73 0.433     1.88       3.58
      # i 4,490 more rows

# conditional_values works with factor by model

    Code
      print(cv)
    Output
      # A tibble: 15 x 8
          .row fac        x2    x0 .fitted   .se .lower_ci .upper_ci
         <int> <fct>   <dbl> <dbl>   <dbl> <dbl>     <dbl>     <dbl>
       1     1 1     0.00131 0.472   0.262 0.545    -0.806     1.33 
       2     2 1     0.227   0.472   1.43  0.295     0.855     2.01 
       3     3 1     0.446   0.472   1.99  0.292     1.42      2.57 
       4     4 1     0.744   0.472   1.40  0.299     0.816     1.99 
       5     5 1     1.00    0.472  -0.346 0.546    -1.42      0.724
       6     6 2     0.00131 0.472  -2.99  0.570    -4.11     -1.88 
       7     7 2     0.227   0.472  -2.40  0.309    -3.01     -1.80 
       8     8 2     0.446   0.472  -1.76  0.327    -2.40     -1.12 
       9     9 2     0.744   0.472   0.193 0.375    -0.543     0.929
      10    10 2     1.00    0.472   3.86  0.641     2.60      5.11 
      11    11 3     0.00131 0.472  -0.830 0.739    -2.28      0.619
      12    12 3     0.227   0.472   8.06  0.445     7.19      8.93 
      13    13 3     0.446   0.472   3.99  0.489     3.03      4.95 
      14    14 3     0.744   0.472   1.93  0.433     1.09      2.78 
      15    15 3     1.00    0.472  -0.263 1.01     -2.24      1.72 

# conditional_values works with supplied factor levels

    Code
      print(cv)
    Output
      # A tibble: 200 x 8
          .row      x2 fac      x0 .fitted   .se .lower_ci .upper_ci
         <int>   <dbl> <fct> <dbl>   <dbl> <dbl>     <dbl>     <dbl>
       1     1 0.00131 2     0.472  -2.99  0.570   -4.11      -1.88 
       2     2 0.00131 3     0.472  -0.830 0.739   -2.28       0.619
       3     3 0.0114  2     0.472  -2.97  0.537   -4.02      -1.91 
       4     4 0.0114  3     0.472  -0.229 0.653   -1.51       1.05 
       5     5 0.0215  2     0.472  -2.94  0.506   -3.93      -1.95 
       6     6 0.0215  3     0.472   0.371 0.578   -0.762      1.50 
       7     7 0.0316  2     0.472  -2.91  0.477   -3.85      -1.98 
       8     8 0.0316  3     0.472   0.970 0.518   -0.0450     1.98 
       9     9 0.0417  2     0.472  -2.89  0.450   -3.77      -2.00 
      10    10 0.0417  3     0.472   1.56  0.476    0.632      2.50 
      # i 190 more rows
