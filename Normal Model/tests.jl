## A small driver file to run the tests in parallel and combine them
a = @spawn test_Gamma("gamma_exp", 10, n_grid, 8675309000, 1., 1.)
b = @spawn test_Gamma("gamma_exp", 10, n_grid, 5165164290, 1., 1.)
c = @spawn test_Gamma("gamma_exp", 10, n_grid, 5167462266, 1., 1.)

println( fetch(a) )
println( fetch(b) )
println( fetch(c) )

