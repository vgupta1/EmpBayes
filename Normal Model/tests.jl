## A small driver file to run the tests in parallel and combine them
a = @spawn test_OddEven("odd_even_exp", 10, n_grid, 8675309, 2, 2)
b = @spawn test_OddEven("odd_even_exp", 10, n_grid, 5165174290, 2, 2)
c = @spawn test_OddEven("odd_even_exp", 10, n_grid, 7462266, 2, 2)

println( fetch(a) )
println( fetch(b) )
println( fetch(c) )

