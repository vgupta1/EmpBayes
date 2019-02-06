## A small driver file to run the three part tests and combine
#run like this
#julia -p 4 -L testHarness_Paper.jl test3Partb.jl
#ARGS[1] is numRun per batch

using Distributed

const spath = "./Results/3Partb"
const param_path = "./Results/param_3Part.csv"

n_grid = [2^i for i = 5:17]
n_grid = [2^i for i = 5:8]
numRuns = parse(Int, ARGS[1])

#Won't be using the regularization results
#So just choose params to make them fast
start_time = time_ns()
a = @spawn test_ReadData(spath, numRuns, n_grid, 8675309, param_path, Gamma_min=5., Gamma_max=20., Gamma_step = 1.)
b = @spawn test_ReadData(spath, numRuns, n_grid, 5164173, param_path, Gamma_min=5., Gamma_max=20., Gamma_step = 1.)
c = @spawn test_ReadData(spath, numRuns, n_grid, 1234563, param_path, Gamma_min=5., Gamma_max=20., Gamma_step = 1.)
d = @spawn test_ReadData(spath, numRuns, n_grid, 5167463, param_path, Gamma_min=5., Gamma_max=20., Gamma_step = 1.)

######
file_a = fetch(a)
file_b = fetch(b)
file_c = fetch(c)
file_d = fetch(d)

time_stamp = (time_ns() - start_time) * 1e-9


##read everyone in, throw away a line
data, header = readdlm(file_a, ',', header=true)

data_t = readdlm(file_b, ',', skipstart=1)

data_t[:, 1] .+= numRuns
data = vcat(data, data_t)

data_t = readdlm(file_c, ',', skipstart=1)
data_t[:, 1] .+= 2numRuns
data = vcat(data, data_t)

data_t = readdlm(file_d, ',', skipstart=1)
data_t[:, 1] .+= 3numRuns
data = vcat(data, data_t)

#strip the name of file_a to make the numbers better
println("Filea \t", file_a)
println("spath \t", spath)

#indx = search(file_a, spath)[end] + 1
indx =  findfirst(spath, file_a)[end] + 1

f = open("$(spath)_$(file_a[indx:end])_full_$(4numRuns).csv", "w")
writedlm(f,  header, ',')
writedlm(f,  data, ',')
close(f)

println("Num Paths: \t $(numRuns) \t Time:", time_stamp )
