## A small driver file to run the tests in parallel and combine them
tag = "uniform_exp"
a = @spawn test_Uniform(tag, 10, n_grid, 8675309000, 1., 2.)
b = @spawn test_Uniform(tag, 10, n_grid, 5165164290, 1., 2.)
c = @spawn test_Uniform(tag, 10, n_grid, 5167462266, 1., 2.)

file_a = fetch(a)
file_b = fetch(b)
file_c = fetch(c)

println( file_a )
println( file_b )
println( file_c )

##read everyone in, throw away a line
data, header = readcsv(file_a, header=true)

data_t = readcsv(file_b, skipstart=1)
data_t[:, 1] += 10
data = vcat(data, data_t)

data_t = readcsv(file_c, skipstart=1)
data_t[:, 1] += 20
data = vcat(data, data_t)

f = open("$(tag)_parallel_results.csv", "w")
writecsv(f, header)
writecsv(f, data)
close(f)