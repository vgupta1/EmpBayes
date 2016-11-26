## A small driver file to run the tests in parallel and combine them
#run like this
#julia -p 3 -L testHarness.jl tests.jl

n_grid = [2^i for i = 7:17]

tic()
numRuns = 50
a = @spawn test_Gaussian("gaussian", numRuns, n_grid, 8675309000, 2., 0., 2.)
b = @spawn test_Gaussian("gaussian", numRuns, n_grid, 5164174290, 2., 0., 2.)
c = @spawn test_Gaussian("gaussian", numRuns, n_grid, 5167462266, 2., 0., 2.)

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

#strip the name of file_a to make the numbers better
indx= rsearch(file_a, "_")[1]
f = open("$(file_a[1:indx])_parallel_results.csv", "w")
writecsv(f, header)
writecsv(f, data)
close(f)

println("Num Paths: \t $(numRuns) \t Time:", toc() )