## A small driver file to run the saturation tests and combine
#run like this
#julia -p 4 -L testSaturation_Harness.jl test_stauration.jl
#ARGS[1] is numRun per batch

using Distributed

const spath = "./Results/satExp"
const param_path = "./Results/param_portExp_mtn2.csv"

n_grid = [2^i for i = 5:17]
numRuns = parse(Int, ARGS[1])
alpha = parse(Float64, ARGS[2])


n_grid = [2^i for i = 13:14]

start_time = time_ns()
a = @spawn test_saturation(spath, numRuns, n_grid, 8675309, param_path, alpha=alpha)
b = @spawn test_saturation(spath, numRuns, n_grid, 5164290, param_path, alpha=alpha)
c = @spawn test_saturation(spath, numRuns, n_grid, 1236, param_path, alpha=alpha)
d = @spawn test_saturation(spath, numRuns, n_grid, 51666, param_path, alpha=alpha)

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

indx =  findfirst(spath, file_a)[end] + 1

f = open("$(spath)_$(file_a[indx:end])_full_$(round(alpha, digits=2))_$(4numRuns).csv", "w")
writedlm(f,  header, ',')
writedlm(f,  data, ',')
close(f)

println("Num Paths: \t $(numRuns) \t Time:", time_stamp )
