## A small driver file to run the tests in parallel and combine them
a = @spawn test_Gaussian("gaussian_exp", 10, n_grid, 8675309, 3, 0)
b = @spawn test_Gaussian("gaussian_exp", 10, n_grid, 5165174290, 3, 0)
c = @spawn test_Gaussian("gaussian_exp", 10, n_grid, 7462266, 3, 0)

println( fetch(a) )
println( fetch(b) )
println( fetch(c) )

