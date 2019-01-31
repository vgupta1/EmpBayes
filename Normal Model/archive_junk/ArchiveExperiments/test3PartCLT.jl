## A small driver file to run the 3 part tests for CLT and combine
#run like this
#julia -p 4 -L testCLTHarness.jl tests_3partCLT.jl
#pass arguments for things via ARGS[1] is numRun per batch

spath = "./Results/3PartCLT_plot"
N_grid = collect(1:20)
const n = 2^17
numRuns = parse(Int, ARGS[1])
dist_type = ARGS[2]

tic()
a = @spawn test_3PartCLT(spath, numRuns, n, N_grid, 8675309, dist_type)
b = @spawn test_3PartCLT(spath, numRuns, n, N_grid, 5164174290, dist_type)
c = @spawn test_3PartCLT(spath, numRuns, n, N_grid, 5167462266, dist_type)
d = @spawn test_3PartCLT(spath, numRuns, n, N_grid, 123456, dist_type)

######
file_a = fetch(a)
file_b = fetch(b)
file_c = fetch(c)
file_d = fetch(d)

time_stamp = toc()


##read everyone in, throw away a line
data, header = readdlm(file_a, ',', header=true)

data_t = readdlm(file_b, ',', skipstart=1)
data_t[:, 1] += numRuns
data = vcat(data, data_t)

data_t = readdlm(file_c, ',', skipstart=1)
data_t[:, 1] += 2numRuns
data = vcat(data, data_t)

data_t = readdlm(file_d, ',', skipstart=1)
data_t[:, 1] += 3numRuns
data = vcat(data, data_t)

#strip the name of file_a to make the numbers better
println("Filea \t", file_a)
println("spath \t", spath)

indx = search(file_a, spath)[end] + 1

f = open("$(spath)_$(file_a[indx:end])_full_$(4numRuns).csv", "w")
writedlm(f,  header)
writedlm(f,  data)
close(f)

println("Num Paths: \t $(numRuns) \t Time:", time_stamp )