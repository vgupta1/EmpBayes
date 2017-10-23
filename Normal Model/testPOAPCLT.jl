## A small driver file to run the 3 part tests for CLT and combine
#run like this
#julia -p 4 -L testCLTHarness.jl tests_POAPCLT.jl
#pass arguments for things via ARGS[1] is numRun per batch

const spath = "./Results/POAPCLT_plot_presentation"
const param_path = "./Results/param_portExp_mtn1.csv"
N_grid = collect(1:15)
const n = 2^17
numRuns = parse(Int, ARGS[1])
dist_type = ARGS[2]

const n = 2^15

#VG need to update here to pass through all the way to bottom.
# Gamma_min = parse(Float64, ARGS[3])
# Gamma_max = parse(Float64, ARGS[4])

tic()
a = @spawn test_POAPCLT(spath, param_path, numRuns, n, N_grid, 8675309, dist_type)
b = @spawn test_POAPCLT(spath, param_path, numRuns, n, N_grid, 5164174290, dist_type)
c = @spawn test_POAPCLT(spath, param_path, numRuns, n, N_grid, 5167462266, dist_type)
d = @spawn test_POAPCLT(spath, param_path, numRuns, n, N_grid, 123456, dist_type)

######
file_a = fetch(a)
file_b = fetch(b)
file_c = fetch(c)
file_d = fetch(d)

time_stamp = toc()


##read everyone in, throw away a line
data, header = readcsv(file_a, header=true)

data_t = readcsv(file_b, skipstart=1)
data_t[:, 1] += numRuns
data = vcat(data, data_t)

data_t = readcsv(file_c, skipstart=1)
data_t[:, 1] += 2numRuns
data = vcat(data, data_t)

data_t = readcsv(file_d, skipstart=1)
data_t[:, 1] += 3numRuns
data = vcat(data, data_t)

#strip the name of file_a to make the numbers better
println("Filea \t", file_a)
println("spath \t", spath)

indx = search(file_a, spath)[end] + 1

f = open("$(spath)_$(file_a[indx:end])_full_$(4numRuns).csv", "w")
writecsv(f, header)
writecsv(f, data)
close(f)

println("Num Paths: \t $(numRuns) \t Time:", time_stamp )