## A small driver file to run the portfolio tests and combine
#run like this
#julia -p 4 -L testHarness_Paper.jl tests_port.jl
#ARGS[1] is numRun per batch

spath = "./Results/portExp"
const param_path = "./Results/param_portExp_mtn1.csv"

n_grid = [2^i for i = 5:17]

numRuns = parse(Int, ARGS[1])

tic()
a = @spawn test_ReadData(spath, numRuns, n_grid, 8675309000, param_path, Gamma_min=1, Gamma_max=20.)
b = @spawn test_ReadData(spath, numRuns, n_grid, 5164174290, param_path, Gamma_min=1, Gamma_max=20.)
c = @spawn test_ReadData(spath, numRuns, n_grid, 123456, param_path, Gamma_min=1, Gamma_max=20.)
d = @spawn test_ReadData(spath, numRuns, n_grid, 5167462266, param_path, Gamma_min=1, Gamma_max=20.)

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