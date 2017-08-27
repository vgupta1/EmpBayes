## A small driver file to run the tests in parallel and combine them
#run like this
#julia -p 3 -L testHarness.jl tests.jl
#pass arguments for things via ARGS[1] is numRun per batch

n_grid = [2^i for i = 5:17]
numRuns = parse(Int, ARGS[1])

spath = "./Results_Paper/SAA_Plot"
even_vs = float(ARGS[2])

tic()
a = @spawn test_OddEven(spath, numRuns, n_grid, 8675309000, even_vs, includeReg=false)
b = @spawn test_OddEven(spath, numRuns, n_grid, 5164174290, even_vs, includeReg=false)
c = @spawn test_OddEven(spath, numRuns, n_grid, 5167462266, even_vs, includeReg=false)
time_stamp = toc()

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
println("Filea \t", file_a)
println("spath \t", spath)

indx = search(file_a, spath)[end] + 1

f = open("$spath_$(file_a[indx:end])_full_$(3numRuns).csv", "w")
writecsv(f, header)
writecsv(f, data)
close(f)

println("Num Paths: \t $(numRuns) \t Time:", time_stamp )