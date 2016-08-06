## A small driver file to run the tests in parallel and combine them

n_grid = [2^i for i = 8:18]

tag = "oddEven"
numRuns = 40
a = @spawn test_OddEven(tag, numRuns, n_grid, 8675309000, 2., 2.)
b = @spawn test_OddEven(tag, numRuns, n_grid, 5164174290, 2., 2.)
c = @spawn test_OddEven(tag, numRuns, n_grid, 5167462266, 2., 2.)

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

f = open("$(file_a)_parallel_results.csv", "w")
writecsv(f, header)
writecsv(f, data)
close(f)