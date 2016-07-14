## A small driver file to run the tests in parallel and combine them
a = @spawn test_UniformLikelihood(10, n_grid, 8675309, 20, 35, 1, 33)
b = @spawn test_UniformLikelihood(10, n_grid, 8675309, 20, 35, 1, 33)
c = @spawn test_UniformLikelihood(10, n_grid, 8675309, 20, 35, 1, 33)

println( fetch(a) )
println( fetch(b) )
println( fetch(c) )

