##Tests for the Small-Data KP

include("KPNormal.jl")
using Distributions, KP

#Given a simulation object, generaton of data is same model to model
function simData!(o, xs, ys, zs)
	zs[:] = randn!(zs) ./ sqrt(o.taus) + o.thetas
	o.noise[:] = randn!(o.noise) ./ sqrt(o.taus)
	xs[:] = zs + o.noise
	ys[:] = zs - o.noise 
end

#Default implementation
#Leaves cs, thetas and taus fixed in o
sim!(o, x, y, z) = simData!(o, x, y, z)

#cs and taus fixed across all runs
#thetas drawn as gaussians tau0
type NormalBayesExp
	cs
	taus
	thetas
	noise
	tau0
	mu0
end

function NormalBayesExp(seed::Integer, tau0, n_max, mu0, frac_fit=.1, avg_tau=4.)
	srand(seed)
	cs = rand(n_max) * 2./frac_fit   
	taus = 2. * avg_tau * rand(n_max)
	NormalBayesExp(cs, taus, zeros(Float64, n_max), zeros(Float64, n_max), tau0, mu0)
end

#cs taus remain fixed
function sim!(o::NormalBayesExp, x, y, z)
	randn!(o.thetas) 
	o.thetas /= sqrt(o.tau0)
	o.thetas += o.mu0
	simData!(o, x, y, z)
end

type OddEvenExp
	cs
	thetas
	taus
	noise
end

#cs and taus fixed across all runs
function OddEvenExp(seed::Integer, n_max::Integer, odd_theta, odd_tau, frac_fit=.1)
	srand(seed)
	cs = 1./frac_fit * ones(Float64, n_max)

	thetas = ones(Float64, n_max)
	thetas[1:2:n_max] = odd_theta
	taus = ones(Float64, n_max)
	taus[1:2:n_max] = odd_tau

	noise = zeros(Float64, n_max)
	OddEvenExp(cs, thetas, taus, noise)
end

#cs and taus fixed across all runs
#thetas drawn as gamma(alpha, beta)
type GammaExp
	cs
	thetas
	taus
	noise
	alpha
	beta
end

#thetas follow an exponential distribution
#Costs are uniform
function GammaExp(seed::Integer, alpha, beta, n_max, frac_fit = .1, avg_tau=2)
	srand(seed)
	cs = rand(n_max) * 2./frac_fit   
	taus = 2. * avg_tau * rand(n_max)
	GammaExp(cs, zeros(Float64, n_max), taus, zeros(Float64, n_max), alpha, beta)
end

function sim!(o::GammaExp, x, y, z)
	#cs & taus are fixed
	rand!(Gamma(o.alpha, o.beta), o.thetas) 	
	simData!(o, x, y, z)	
end

#cs and taus fixed across all runs
#thetas drawn as uniform [a, b]
type UniformExp
	cs
	thetas
	taus
	noise
	a
	b
end

function UniformExp(seed::Integer, a, b, n_max, frac_fit = .1, avg_tau=2)
	srand(seed)
	cs = rand(n_max) * 2./frac_fit   
	taus = 2. * avg_tau * rand(n_max)
	UniformExp(cs, zeros(Float64, n_max), taus, zeros(Float64, n_max), a, b)
end

function sim!(o::UniformExp, x, y, z)
	#cs & taus are fixed
	rand!(Uniform(o.a, o.b), o.thetas) 	
	simData!(o, x, y, z)	
end


#cs and taus fixed across all runs
#thetas drawn as beta[a, b]
type BetaExp
	cs
	thetas
	taus
	noise
	a
	b
end

function BetaExp(seed::Integer, a, b, n_max, frac_fit = .1, avg_tau=2)
	srand(seed)
	cs = rand(n_max) * 2./frac_fit   
	taus = 2. * avg_tau * rand(n_max)
	BetaExp(cs, zeros(Float64, n_max), taus, zeros(Float64, n_max), a, b)
end

function sim!(o::BetaExp, x, y, z)
	#cs & taus are fixed
	rand!(Beta(o.a, o.b), o.thetas) 	
	simData!(o, x, y, z)	
end

#######################
### The usual bayesian normal setup
### costs and taus are random, but fixed over runs
function test_Gaussian(file_out, numRuns, n_grid, seed, tau0)
	#build the sim object
	o = NormalBayesExp(seed, tau0, maximum(n_grid), 0)

	#run the testharness
	#output files
	f = open("$(file_out)_$(tau0)_$(seed).csv", "w")
	writecsv(f, ["Run" "n" "Method" "YVal" "thetaVal" "time" "tau0"])
	test_harness(f, numRuns, o, n_grid)
	close(f)
end

### Bayesian set-up with a non-zero mean
### costs and taus are random, but fixed over runs
function test_Gaussian(file_out, numRuns, n_grid, seed, tau0, mu0)
	#build the sim object
	o = NormalBayesExp(seed, tau0, maximum(n_grid), mu0)

	#run the testharness
	#output files
	f = open("$(file_out)_$(mu0)_$(tau0)_$(seed).csv", "w")
	writecsv(f, ["Run" "n" "Method" "YVal" "thetaVal" "time" "tau0"])
	test_harness(f, numRuns, o, n_grid)
	close(f)
end

### The odd even set-up.  
function test_OddEven(file_out, numRuns, n_grid, seed, odd_theta, odd_tau)
	#build the sim object
	o = OddEvenExp(seed, maximum(n_grid), odd_theta, odd_tau)

	#run the testharness
	#output files
	f = open("$(file_out)_$(odd_theta)_$(odd_tau)_$(seed).csv", "w")
	writecsv(f, ["Run" "n" "Method" "YVal" "thetaVal" "time" "tau0"])
	test_harness(f, numRuns, o, n_grid)
	close(f)
end

### The Gamma Test
function test_Gamma(file_out, numRuns, n_grid, seed, alpha, beta)
	o = GammaExp(seed, alpha, beta, maximum(n_grid))

	f = open("$(file_out)_Gamma_$(alpha)_$(beta)_$(seed).csv", "w")
	writecsv(f, ["Run" "n" "Method" "YVal" "thetaVal" "time" "tau0"])
	test_harness(f, numRuns, o, n_grid)
	close(f)
end

### The Uniform Test
function test_Uniform(file_out, numRuns, n_grid, seed, a, b)
	o = UniformExp(seed, a, b, maximum(n_grid))

	f = open("$(file_out)_Uniform_$(a)_$(b)_$(seed).csv", "w")
	writecsv(f, ["Run" "n" "Method" "YVal" "thetaVal" "time" "tau0"])
	test_harness(f, numRuns, o, n_grid)
	close(f)
end

### The Beta Test
function test_Beta(file_out, numRuns, n_grid, seed, a, b)
	o = BetaExp(seed, a, b, maximum(n_grid))

	f = open("$(file_out)_Beta_$(a)_$(b)_$(seed).csv", "w")
	writecsv(f, ["Run" "n" "Method" "YVal" "thetaVal" "time" "tau0"])
	test_harness(f, numRuns, o, n_grid)
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
function test_harness(f, numRuns, o, n_grid)
	const n_max = maximum(n_grid)
	xs = zeros(Float64, n_max)
	ys = zeros(Float64, n_max)
	zs = zeros(Float64, n_max)

	for iRun = 1:numRuns
		#generate the entire path up to n_max
		sim!(o, xs, ys, zs)

		for n in n_grid
			#Compute performance of each method
			#Naive x
			tic()
			qs = q(o.cs[1:n], xs[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "NaiveX" yval thetaval t 0.])

			#Naive Z
			#Same as the ZZ method
			tic()
			qs = q(o.cs[1:n], zs[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "NaiveZ" yval thetaval t 0.])

			#fullInfo val
			tic()
			qs = q(o.cs[1:n], o.thetas[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "FullInfo" yval thetaval t 0.])

			#The "Bayes" value.. only possible bc we if we we are in bayesian
			if :tau0 in fieldnames(o)
				tic()
				qs = q(o.cs[1:n], shrink(zs[1:n], o.taus[1:n], o.tau0))
				t = toc()
				yval = dot(ys[1:n], qs)/n
				thetaval = dot(o.thetas[1:n], qs)/n
				writecsv(f, [iRun n "Bayes" yval thetaval t o.tau0])
			end

			#Tau MLE
			tic()
			tauMLE, qs = q_MLE(o.cs[1:n], zs[1:n], o.taus[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "MLE" yval thetaval t tauMLE])

			#Tau MM
			tic()
			tauMM, qs = q_MM(o.cs[1:n], zs[1:n], o.taus[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "MM" yval thetaval t tauMM])

			#Our sample split method x
			tic()
			qs, vals, objs = best_q_tau(o.cs[1:n], xs[1:n], o.taus[1:n], ys[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "EmpBayesX" yval thetaval t vals[indmax(objs)]])

			#rescaling the sample_split method
			tau_RS = vals[indmax(objs)]/2
			qs = q(o.cs[1:n], shrink(zs[1:n], o.taus[1:n], tau_RS))
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "Rescaled" yval thetaval t tau_RS])

			#Some rescalings by wrong factors
			for i = -3:3
				i == -1 && continue  #already computed the rescaled value
				tau_RS = vals[indmax(objs)]/ 2.0^i
				qs = q(o.cs[1:n], shrink(zs[1:n], o.taus[1:n], tau_RS))
				yval = dot(ys[1:n], qs)/n
				thetaval = dot(o.thetas[1:n], qs)/n
				writecsv(f, [iRun n "Rescaled_$(round(2.0^i, 3))" yval thetaval t tau_RS])
			end

			#Oracle X method
			tic()
			qs, vals, objs = best_q_tau(o.cs[1:n], xs[1:n], o.taus[1:n], o.thetas[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "OracleX" yval thetaval t vals[indmax(objs)]])

			#Oracle Z method
			tic()
			qs, vals, objs = best_q_tau(o.cs[1:n], zs[1:n], o.taus[1:n], o.thetas[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(o.thetas[1:n], qs)/n
			@assert abs(thetaval - maximum(objs)) <= 1e-12 "Weird Mismatch? \t $thetaval \t $(maximum(objs))"
			writecsv(f, [iRun n "OracleZ" yval thetaval t vals[indmax(objs)]])

		end
		flush(f)
	end
end

n_grid = [2^i for i = 8:17]

#run some small examples to pre-compile for optimization
test_Gaussian("temp_Gaussian2", 5, [100, 150], 87, 3, 5)
test_Gaussian("temp", 20, n_grid, 8675309, 2, 0)


test_Gaussian("temp_Gaussian", 5, [100, 150], 87, 3)
test_OddEven("temp_OddEven", 5,[100, 150], 87, 2, 2)
test_Gamma("temp_Gamma", 5, [100, 150], 87, 1., 1.)
test_Uniform("temp_Uniform", 5, [100, 150], 87, 1, 2)
test_Gaussian("temp_Gaussian2", 5, [100, 150], 87, 3, 5)
test_Beta("temp_Uniform", 5, [100, 150], 87, .5, .5)
test_Gaussian("temp_Gaussian2", 5, [100, 150], 87, 3, 5)