## A small driver file to run the 3 part tests and combine
#run like this
#julia -p 3 -L testHarness_Paper.jl tests_3part.jl
#pass arguments for things via ARGS[1] is numRun per batch

spath = "./Results/3Part_plot"
n_grid = [2^i for i = 5:17]

numRuns = parse(Int, ARGS[1])

#extract the relevant things from ARGS
theta_l, v_l, c_l = ARGS[2:4]
theta_m, v_m, c_m = ARGS[5:7]
theta_h, v_h, c_h = ARGS[8:10]

tic()
a = @spawn test_threePart(spath, numRuns, n_grid, 8675309000, 
						float(theta_l), float(v_l), float(c_l), float(theta_m), float(v_m), float(c_m), float(theta_h), float(v_h), float(c_h), true)
b = @spawn test_threePart(spath, numRuns, n_grid, 5164174290, 
						float(theta_l), float(v_l), float(c_l), float(theta_m), float(v_m), float(c_m), float(theta_h), float(v_h), float(c_h), true)
c = @spawn test_threePart(spath, numRuns, n_grid, 5167462266, 
						float(theta_l), float(v_l), float(c_l), float(theta_m), float(v_m), float(c_m), float(theta_h), float(v_h), float(c_h), true)
d = @spawn test_threePart(spath, numRuns, n_grid, 123456, 
						float(theta_l), float(v_l), float(c_l), float(theta_m), float(v_m), float(c_m), float(theta_h), float(v_h), float(c_h), true)
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