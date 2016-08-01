## A small driver file to run the tests in parallel and combine them
tag = "gaussianExp"
numRuns = 15
a = @spawn test_Gaussian(tag, numRuns, n_grid, 1675309000, 3., 0.)
b = @spawn test_Gaussian(tag, numRuns, n_grid, 2165164290, 3., 0.)
c = @spawn test_Gaussian(tag, numRuns, n_grid, 3167462266, 3., 0.)
######

file_a = fetch(a)
file_b = fetch(b)
file_c = fetch(c)

println( file_a )
println( file_b )
println( file_c )

##read everyone in, throw away a line
data, header = readcsv(file_a, header=true)

data_t = readcsv(file_b, skipstart=1)
data_t[:, 1] += numRuns
data = vcat(data, data_t)

data_t = readcsv(file_c, skipstart=1)
data_t[:, 1] += 2numRuns
data = vcat(data, data_t)

f = open("$(tag)_parallel_results.csv", "w")
writecsv(f, header)
writecsv(f, data)
close(f)