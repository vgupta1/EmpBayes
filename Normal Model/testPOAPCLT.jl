## A small driver file to run the POAP Robustness to normality tests and combine
#run like this
#julia -p 4 -L testCLTHarness.jl testPOAPCLT.jl
#pass arguments for things via ARGS[1] is numRun per batch

spath = "./Results/POAPCLT_plot_"
param_path = "./Results/param_portExp_mtn2.csv"
S_grid = collect(1:2:26)
n = 2^17
numRuns = parse(Int, ARGS[1])
dist_type = ARGS[2]

#"true" flag indicates use constantPrecision, i.e., the robustness to Normality experiments

start_time = time_ns()
a = @spawn test_POAPCLT(spath, param_path, numRuns, n, S_grid, 8675309, dist_type, 1.0, 100., 1., onlyOracle=true)
b = @spawn test_POAPCLT(spath, param_path, numRuns, n, S_grid, 5164174290, dist_type, 1.0, 100., 1., onlyOracle=true)
c = @spawn test_POAPCLT(spath, param_path, numRuns, n, S_grid, 5167462266, dist_type, 1.0, 100., 1., onlyOracle=true)
d = @spawn test_POAPCLT(spath, param_path, numRuns, n, S_grid, 123456, dist_type, 1.0, 100., 1., onlyOracle=true)

######
file_a = fetch(a)
file_b = fetch(b)
file_c = fetch(c)
file_d = fetch(d)

time_stamp = (time_ns() - start_time) * 1e-9

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
writedlm(f,  header, ',')
writedlm(f,  data, ',')
close(f)

println("Num Paths: \t $(numRuns) \t Time:", time_stamp )