##Tests for the Small-Data KP

include("KPNormal.jl")
using Distributions, KP

#define simulation objects that describe how costs/thetas/taus are gen for large n

#cs and taus fixed across all runs
#thetas drawn as gaussians tau0
type NormalBayesExp
	cs::Vector{Float64}
	taus::Vector{Float64}
	sigs::Vector{Float64}
	thetas::Vector{Float64}
	xs::Vector{Float64}
	ys::Vector{Float64}
	zs::Vector{Float64}
	noise::Vector{Float64}
	n_grid::Vector{Int64}
	tau0::Float64
end

function NormalBayesExp(seed, tau0, n_grid, frac_fit = .1, avg_tau = 4.)
	srand(seed)
	const n_max = maximum(n_grid)
	cs = rand(n_max) * 2./frac_fit   
	taus = 2. * avg_tau * rand(n_max)
	sigs = 1./sqrt(taus)
	thetas = zeros(Float64, n_max)
	xs = zeros(Float64, n_max)
	ys = zeros(Float64, n_max)
	zs = zeros(Float64, n_max)
	noise = zeros(Float64, n_max)
	NormalBayesExp(cs, taus, sigs, thetas, xs, ys, zs, noise, n_grid, tau0)
end

function sim!(o::NormalBayesExp)
	#cs taus and sigs remain fixed
	randn!(o.thetas) 
	o.thetas /= sqrt(o.tau0)
	o.zs[:] = randn!(o.zs) .* o.sigs + o.thetas
	o.noise[:] = randn!(o.noise) .* o.sigs
	o.xs[:] = o.zs + o.noise
	o.ys[:] = o.zs - o.noise 
end

#cs and taus fixed across all runs
#thetas drawn as gaussians tau0
type OddEvenExp
	cs::Vector{Float64}
	thetas::Vector{Float64}
	taus::Vector{Float64}
	xs::Vector{Float64}
	ys::Vector{Float64}
	zs::Vector{Float64}
	noise::Vector{Float64}
	n_grid::Vector{Int64}
end

function OddEvenExp(seed, n_grid, odd_theta, odd_tau, frac_fit = .1)
	srand(seed)
	const n_max = maximum(n_grid)
	cs = 1./frac_fit * ones(Float64, n_max)

	thetas = ones(Float64, n_max)
	thetas[1:2:n_max] = odd_theta
	taus = ones(Float64, n_max)
	taus[1:2:n_max] = odd_tau

	xs = zeros(Float64, n_max)
	ys = zeros(Float64, n_max)
	zs = zeros(Float64, n_max)
	noise = zeros(Float64, n_max)
	OddEvenExp(cs, thetas, taus, xs, ys, zs, noise, n_grid)
end

function sim!(o::OddEvenExp)
	#cs, taus, and thetas remain fixed
	const n_max = maximum(o.n_grid)
	o.zs[:] = randn!(o.zs) ./ sqrt(o.taus) + o.thetas
	o.noise[:] = randn!(o.noise) ./ sqrt(o.taus)
	o.xs[:] = o.zs + o.noise
	o.ys[:] = o.zs - o.noise 
end


#cs and taus fixed across all runs
#thetas drawn as gaussians tau0
type ExponentialExp
	cs::Vector{Float64}
	thetas::Vector{Float64}
	taus::Vector{Float64}
	xs::Vector{Float64}
	ys::Vector{Float64}
	zs::Vector{Float64}
	noise::Vector{Float64}
	n_grid::Vector{Int64}
	alpha::Float64
	beta::Float64
end

#generate thetas once according to an exponential distribution
#costs are random but fixed.  
function ExponentialExp(seed, alpha, beta, n_grid, frac_fit = .1)
	srand(seed)
	const n_max = maximum(n_grid)
	cs = rand(n_max) * 2./frac_fit   

	thetas = rand(Gamma(alpha, beta), n_max)
	taus = thetas.^2

	xs = zeros(Float64, n_max)
	ys = zeros(Float64, n_max)
	zs = zeros(Float64, n_max)
	noise = zeros(Float64, n_max)
	ExponentialExp(cs, thetas, taus, xs, ys, zs, noise, n_grid, alpha, beta)
end

function sim!(o::ExponentialExp)
	#cs, taus, and thetas remain fixed
	const n_max = maximum(o.n_grid)
	o.zs[:] = randn!(o.zs) ./ sqrt(o.taus) + o.thetas
	o.noise[:] = randn!(o.noise) ./ sqrt(o.taus)
	o.xs[:] = o.zs + o.noise
	o.ys[:] = o.zs - o.noise 
end



### The usual bayesian normal setup
### costs and taus are random, but fixed over runs
function test1_a(file_out, numRuns, n_grid, seed, tau0)
	#build the sim object
	o = NormalBayesExp(seed, tau0, n_grid)

	#run the testharness
	#output files
	f = open("$(file_out)_$(tau0)_$(seed).csv", "w")
	writecsv(f, ["Run" "n" "Method" "YVal" "thetaVal" "time" "tau0"])
	test_harness(f, numRuns, o)
	close(f)
end


### The odd even set-up.  
function test2_a(file_out, numRuns, n_grid, seed, odd_theta, odd_tau)
	#build the sim object
	o = OddEvenExp(seed, n_grid, odd_theta, odd_tau)

	#run the testharness
	#output files
	f = open("$(file_out)_$(odd_theta)_$(odd_tau)_$(seed).csv", "w")
	writecsv(f, ["Run" "n" "Method" "YVal" "thetaVal" "time" "tau0"])
	test_harness(f, numRuns, o)
	close(f)
end


### The Exponential Test
function test3_a(file_out, numRuns, n_grid, seed, alpha, beta)
	o = ExponentialExp(seed, alpha, beta, n_grid)

	f = open("$(file_out)_Exponential_$(alpha)_$(beta)_$(seed).csv", "w")
	writecsv(f, ["Run" "n" "Method" "YVal" "thetaVal" "time" "tau0"])
	test_harness(f, numRuns, o)
	close(f)

end


######################
####
#For increasing n, many simulations compare
#		#Naive x method
#		#Naive z method
#		#FullInfo Zstar
#		#Our method X
#		#Our method Avg X Y
#		#Oracle x method
#		#Oracle Z method
#		#Best randomized solution with n^k choices k = 2, 4, 6
#Compare on metrics
#		#Value wrt theta
#		#Time to compute
#		#optimal value of tau0
function test_harness(f, numRuns, o)
	for iRun = 1:numRuns
		#generate the entire path up to n_max
		sim!(o)

		for n in o.n_grid
			#Compute performance of each method
			#Naive x
			tic()
			qs = q(o.cs[1:n], o.xs[1:n])
			t = toc()
			yval = dot(o.ys[1:n], qs)/n
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "NaiveX" yval thetaval t 0.])

			#Naive Z
			tic()
			qs = q(o.cs[1:n], o.zs[1:n])
			t = toc()
			yval = dot(o.ys[1:n], qs)/n
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "NaiveZ" yval thetaval t 0.])

			#fullInfo val
			tic()
			qs = q(o.cs[1:n], o.thetas[1:n])
			t = toc()
			yval = dot(o.ys[1:n], qs)/n
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "FullInfo" yval thetaval t 0.])

			#The "Bayes" value.. only possible bc we if we we are in bayesian
			if :tau0 in fieldnames(o)
				tic()
				qs = q(o.cs[1:n], shrink(o.zs[1:n], o.taus[1:n], o.tau0))
				t = toc()
				yval = dot(o.ys[1:n], qs)/n
				thetaval = dot(o.thetas[1:n], qs)/n
				writecsv(f, [iRun n "Bayes" yval thetaval t o.tau0])
			end

			#Tau MLE
			tic()
			tauMLE, qs = q_MLE(o.cs[1:n], o.zs[1:n], o.taus[1:n])
			t = toc()
			yval = dot(o.ys[1:n], qs)/n
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "MLE" yval thetaval t tauMLE])

			#Tau MM
			tic()
			tauMM, qs = q_MM(o.cs[1:n], o.zs[1:n], o.taus[1:n])
			t = toc()
			yval = dot(o.ys[1:n], qs)/n
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "MM" yval thetaval t tauMM])

			#Our sample split method x
			tic()
			qs, vals, objs = best_q_tau(o.cs[1:n], o.xs[1:n], o.taus[1:n], o.ys[1:n])
			t = toc()
			yval = dot(o.ys[1:n], qs)/n
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "EmpBayesX" yval thetaval t vals[indmax(objs)]])

			#rescaling the sample_split method
			tau_RS = vals[indmax(objs)]/2
			qs = q(o.cs[1:n], shrink(o.zs[1:n], o.taus[1:n], tau_RS))
			yval = dot(o.ys[1:n], qs)/n
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "Rescaled" yval thetaval t tau_RS])

			#Some rescalings by wrong factors
			for i = -3:3
				i == -1 && continue  #already computed the rescaled value
				taus_RS = vals[indmax(objs)]/ 2.0^i
				qs = q(o.cs[1:n], shrink(o.zs[1:n], o.taus[1:n], tau_RS))
				yval = dot(o.ys[1:n], qs)/n
				thetaval = dot(o.thetas[1:n], qs)/n
				writecsv(f, [iRun n "Rescaled_$(round(taus_RS, 3))" yval thetaval t tau_RS])
			end

			#The ZZ method
			tic()
			qs, vals, objs = best_q_tau(o.cs[1:n], o.zs[1:n], o.taus[1:n], o.zs[1:n])
			t = toc()
			yval = dot(o.ys[1:n], qs)/n
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "ZZ" yval thetaval t vals[indmax(objs)]])

			#Oracle X method
			tic()
			qs, vals, objs = best_q_tau(o.cs[1:n], o.xs[1:n], o.taus[1:n], o.thetas[1:n])
			t = toc()
			yval = dot(o.ys[1:n], qs)/n
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "OracleX" yval thetaval t vals[indmax(objs)]])

			#Oracle Z method
			tic()
			qs, vals, objs = best_q_tau(o.cs[1:n], o.zs[1:n], o.taus[1:n], o.thetas[1:n])
			t = toc()
			yval = dot(o.ys[1:n], qs)/n
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "OracleZ" yval thetaval t vals[indmax(objs)]])

		end
		flush(f)
	end
end

n_grid = [2^i for i = 8:12]

#run some small examples to pre-compile for optimization
#test1_a("temp", 5, [100, 150], 87, 3)
#test2_a("temp2", 5,[100, 150], 87, 2, 2)
#test3_a("temp3", 5, [100, 150], 87, 1., 1.)
