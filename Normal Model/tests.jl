## A small driver file to run the tests in parallel and combine them
a = @spawn test_UniformLikelihood(10, n_grid, 8675309, 10., 15., 1., 14.1)
b = @spawn test_UniformLikelihood(10, n_grid, 8675309, 10., 15., 1., 14.1)
c = @spawn test_UniformLikelihood(10, n_grid, 8675309, 10., 15., 1., 14.1)

println( fetch(a) )
println( fetch(b) )
println( fetch(c) )

