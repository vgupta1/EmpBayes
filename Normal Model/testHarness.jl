##Tests for the Small-Data KP
include("KPNormal.jl")
using Distributions, KP

#Default simulation of normal random variates
function simData!(o, xs, ys, zs)
	zs[:] = randn!(zs) ./ sqrt(o.vs) + o.thetas
	o.noise[:] = randn!(o.noise) ./ sqrt(o.vs)
	xs[:] = zs + o.noise
	ys[:] = zs - o.noise 
end

#Default implementation
#Leaves cs, thetas and vs fixed in o
sim!(o, x, y, z) = simData!(o, x, y, z)

#Following simulates N rvs according to the given distributions
#AS N increases, this should approximately have gaussian likelihood
type CLTExp
	cs
	vs
	thetas
	dists
	N
end

#thetas are generated as gaussian
#widths are randomly genrated between [0, width_max]
#intervals are centered on thetas
function CLTExp(seed::Integer, width_max, n_max, N; tau0=2, frac_fit = .1)
	srand(seed)
	@assert rem(N, 2) == 0
	cs = rand(n_max) * 2./frac_fit   
	widths = rand(n_max) * width_max + .1
	vs = 12N ./ widths.^2
	thetas = randn(n_max) ./ sqrt(tau0)

	dists = Array(Distributions.Uniform, n_max)
	for ix = 1:n_max
		dists[ix] = Uniform(thetas[ix] - widths[ix]/2, thetas[ix] + widths[ix]/2)
	end
	CLTExp(cs, vs, thetas, dists, N)
end

function sim!(o::CLTExp, xs, ys, zs)	
	const n_max = length(o.cs)
	zetas = zeros(o.N)
	const half_N = round(Int, o.N/2)
	for ix = 1:n_max
		rand!(o.dists[ix], zetas)
		xs[ix] = mean(zetas[1:half_N])
		ys[ix] = mean(zetas[(half_N + 1):o.N])
		zs[ix] = mean(zetas)
	end
end

## The following constitutes a counter-example to tauxy/2 policy
type UniformLikelihoodExp
	cs
	vs
	thetas
	dist_odd
	dist_even
end

#even numbers have higher mean, lower std. dev
function UniformLikelihoodExp(seed::Integer, theta_even, theta_odd, width_even, width_odd, n_max, lamp)
	srand(seed)
	@assert rem(n_max, 2) == 0

	#compute the value of c and populate
	@assert theta_even + width_even < lamp < theta_odd + width_odd
	c = 4 * width_odd / (theta_odd + width_odd - lamp)	
	println("C Value:\t", c)
	cs = c * ones(n_max)

	#check that we've picked values corresponding to the "interesting" case
	@assert theta_odd + width_odd > theta_even + width_even
	thetas = theta_odd * ones(n_max)
	thetas[2:2:n_max] = theta_even
	widths = width_odd * ones(n_max)
	widths[2:2:n_max] = width_even

	#create the distributions for later
	dist_odd = Uniform(theta_odd - width_odd, theta_odd + width_odd)
	dist_even = Uniform(theta_even - width_even, theta_even + width_even)
	UniformLikelihoodExp(cs, widths.^(-2),  thetas, dist_odd, dist_even)
end

##Specialized implementation for non-gaussian likelihood
function sim!(o::UniformLikelihoodExp, xs, ys, zs)	
	n_max = length(o.cs)
	half_n_max = round(Int, n_max/2)

	xs[1:2:n_max] = rand(o.dist_odd, half_n_max)
	xs[2:2:n_max] = rand(o.dist_even, half_n_max)
	ys[1:2:n_max] = rand(o.dist_odd, half_n_max)
	ys[2:2:n_max] = rand(o.dist_even, half_n_max)
	zs[:] = .5 * (xs + ys);
end

#cs and vs fixed across all runs
#thetas drawn as gaussians tau0
type NormalBayesExp
	cs
	vs
	thetas
	noise
	tau0
	mu0
end

function NormalBayesExp(seed::Integer, tau0, n_max, mu0, frac_fit=.1, avg_tau=4.)
	srand(seed)
	cs = rand(n_max) * 2./frac_fit   
	vs = 2. * avg_tau * rand(n_max)
	NormalBayesExp(cs, vs, zeros(Float64, n_max), zeros(Float64, n_max), tau0, mu0)
end

#cs vs remain fixed
function sim!(o::NormalBayesExp, x, y, z)
	randn!(o.thetas) 
	o.thetas /= sqrt(o.tau0)
	o.thetas += o.mu0
	simData!(o, x, y, z)
end

type OddEvenExp
	cs
	thetas
	vs
	noise
end

#cs and vs fixed across all runs
function OddEvenExp(seed::Integer, n_max::Integer, odd_theta, odd_tau, frac_fit=.1)
	srand(seed)
	cs = 1./frac_fit * ones(Float64, n_max)

	thetas = ones(Float64, n_max)
	thetas[1:2:n_max] = odd_theta
	vs = ones(Float64, n_max)
	vs[1:2:n_max] = odd_tau

	noise = zeros(Float64, n_max)
	OddEvenExp(cs, thetas, vs, noise)
end

#cs and vs fixed across all runs
#thetas drawn as gamma(alpha, beta)
type GammaExp
	cs
	thetas
	vs
	noise
	alpha
	beta
end

#thetas follow an exponential distribution
#Costs are uniform
function GammaExp(seed::Integer, alpha, beta, n_max, frac_fit = .1, avg_tau=2)
	srand(seed)
	cs = rand(n_max) * 2./frac_fit   
	vs = 2. * avg_tau * rand(n_max)
	GammaExp(cs, zeros(Float64, n_max), vs, zeros(Float64, n_max), alpha, beta)
end

function sim!(o::GammaExp, x, y, z)
	#cs & vs are fixed
	rand!(Gamma(o.alpha, o.beta), o.thetas) 	
	simData!(o, x, y, z)	
end

#cs and vs fixed across all runs
#thetas drawn as uniform [a, b]
type UniformExp
	cs
	thetas
	vs
	noise
	a
	b
end

function UniformExp(seed::Integer, a, b, n_max, frac_fit = .1, avg_tau=2)
	srand(seed)
	cs = rand(n_max) * 2./frac_fit   
	vs = 2. * avg_tau * rand(n_max)
	UniformExp(cs, zeros(Float64, n_max), vs, zeros(Float64, n_max), a, b)
end

function sim!(o::UniformExp, x, y, z)
	#cs & vs are fixed
	rand!(Uniform(o.a, o.b), o.thetas) 	
	simData!(o, x, y, z)	
end


#cs and vs fixed across all runs
#thetas drawn as beta[a, b]
type BetaExp
	cs
	thetas
	vs
	noise
	a
	b
end

function BetaExp(seed::Integer, a, b, n_max, frac_fit = .1, avg_tau=2)
	srand(seed)
	cs = rand(n_max) * 2./frac_fit   
	vs = 2. * avg_tau * rand(n_max)
	BetaExp(cs, zeros(Float64, n_max), vs, zeros(Float64, n_max), a, b)
end

function sim!(o::BetaExp, x, y, z)
	#cs & vs are fixed
	rand!(Beta(o.a, o.b), o.thetas) 	
	simData!(o, x, y, z)	
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
				qs = q(o.cs[1:n], shrink(zs[1:n], o.vs[1:n], o.tau0))
				t = toc()
				yval = dot(ys[1:n], qs)/n
				thetaval = dot(o.thetas[1:n], qs)/n
				writecsv(f, [iRun n "Bayes" yval thetaval t o.tau0])
			end

			#Tau MLE
			tic()
			tauMLE, qs = q_MLE(o.cs[1:n], zs[1:n], o.vs[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "MLE" yval thetaval t tauMLE])

			#Tau MM
			tic()
			tauMM, qs = q_MM(o.cs[1:n], zs[1:n], o.vs[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "MM" yval thetaval t tauMM])

			#Our sample split method x
			tic()
			qs, vals, objs = best_q_tau(o.cs[1:n], xs[1:n], o.vs[1:n], ys[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "EmpBayesX" yval thetaval t vals[indmax(objs)]])

			#rescaling the sample_split method
			tau_RS = vals[indmax(objs)]/2
			qs = q(o.cs[1:n], shrink(zs[1:n], o.vs[1:n], tau_RS))
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "Rescaled" yval thetaval t tau_RS])

			# #Some rescalings by wrong factors
			# for i = -3:3
			# 	i == -1 && continue  #already computed the rescaled value
			# 	tau_RS = vals[indmax(objs)]/ 2.0^i
			# 	qs = q(o.cs[1:n], shrink(zs[1:n], o.vs[1:n], tau_RS))
			# 	yval = dot(ys[1:n], qs)/n
			# 	thetaval = dot(o.thetas[1:n], qs)/n
			# 	writecsv(f, [iRun n "Rescaled_$(round(2.0^i, 3))" yval thetaval t tau_RS])
			# end

			#The idealized stein approach
			tic()
			qs, vals, objs = KP.stein_q_tau_exact(o.cs[1:n], zs[1:n], o.vs[1:n], o.thetas[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "ExactStein" yval thetaval t vals[indmax(objs)]])

			#The impuse stein approach
			#use the optimized rate, i.e. h_n = n^-1/6
			h = n^-.16666
			tic()
			qs, vals, objs = KP.stein_q_tau_impulse(o.cs[1:n], zs[1:n], o.vs[1:n], h, tau_step = .1)
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "ImpulseStein_6" yval thetaval t vals[indmax(objs)]])


			h = n^-.33333
			tic()
			qs, vals, objs = KP.stein_q_tau_impulse(o.cs[1:n], zs[1:n], o.vs[1:n], h, tau_step = .1)
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "ImpulseStein_3" yval thetaval t vals[indmax(objs)]])

			h = n^-.111111
			tic()
			qs, vals, objs = KP.stein_q_tau_impulse(o.cs[1:n], zs[1:n], o.vs[1:n], h, tau_step = .1)
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "ImpulseStein_9" yval thetaval t vals[indmax(objs)]])

			h = n^-.45
			tic()
			qs, vals, objs = KP.stein_q_tau_impulse(o.cs[1:n], zs[1:n], o.vs[1:n], h, tau_step = .1)
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "ImpulseStein_45" yval thetaval t vals[indmax(objs)]])




			# #The primal stein approach
			# #use the optimized rate, i.e. h_n = n^-1/6
			# h = n^-.16666
			# tic()
			# qs, vals, objs = KP.stein_q_tau_primal(o.cs[1:n], zs[1:n], o.vs[1:n], h)
			# t = toc()
			# yval = dot(ys[1:n], qs)/n
			# thetaval = dot(o.thetas[1:n], qs)/n
			# writecsv(f, [iRun n "PrimalStein" yval thetaval t vals[indmax(objs)]])

			#Oracle X method
			tic()
			qs, vals, objs = best_q_tau(o.cs[1:n], xs[1:n], o.vs[1:n], o.thetas[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "OracleX" yval thetaval t vals[indmax(objs)]])

			#Oracle Z method
			tic()
			qs, vals, objs = best_q_tau(o.cs[1:n], zs[1:n], o.vs[1:n], o.thetas[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(o.thetas[1:n], qs)/n
			@assert abs(thetaval - maximum(objs)) <= 1e-4 "Weird Mismatch? \t $thetaval \t $(maximum(objs))"
			writecsv(f, [iRun n "OracleZ" yval thetaval t vals[indmax(objs)]])

			#Ridge Proxy
			tic()
			qs = q_ridge(o.cs[1:n], xs[1:n], ys[1:n], o.vs[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "Ridge" yval thetaval t vals[indmax(objs)]])

			#weighted l2 regularization  mu = .1
			mu = .1
			tic()
			qs = q_l2reg(o.cs[1:n], zs[1:n], o.vs[1:n], mu)[1]
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "Regularization_.1" yval thetaval t vals[indmax(objs)]])

			mu = 1
			tic()
			qs = q_l2reg(o.cs[1:n], zs[1:n], o.vs[1:n], mu)[1]
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "Regularization_1" yval thetaval t vals[indmax(objs)]])

			mu = .01
			tic()
			qs = q_l2reg(o.cs[1:n], zs[1:n], o.vs[1:n], mu)[1]
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "Regularization_.01" yval thetaval t vals[indmax(objs)]])

		end
		flush(f)
	end
end

#######################
### The usual bayesian normal setup
### costs and vs are random, but fixed over runs
function test_Gaussian(file_out, numRuns, n_grid, seed, tau0)
	#build the sim object
	o = NormalBayesExp(seed, tau0, maximum(n_grid), 0)

	#run the testharness
	#output files
	f = open("$(file_out)_$(tau0)_$(seed).csv", "w")
	writecsv(f, ["Run" "n" "Method" "YVal" "thetaVal" "time" "tau0"])
	test_harness(f, numRuns, o, n_grid)
	close(f)
	return "$(file_out)_$(tau0)_$(seed).csv"
end

### Bayesian set-up with a non-zero mean
### costs and vs are random, but fixed over runs
function test_Gaussian(file_out, numRuns, n_grid, seed, tau0, mu0)
	#build the sim object
	o = NormalBayesExp(seed, tau0, maximum(n_grid), mu0)

	#run the testharness
	#output files
	f = open("$(file_out)_$(mu0)_$(tau0)_$(seed).csv", "w")
	writecsv(f, ["Run" "n" "Method" "YVal" "thetaVal" "time" "tau0"])
	test_harness(f, numRuns, o, n_grid)
	close(f)
	return "$(file_out)_$(mu0)_$(tau0)_$(seed).csv"
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
	return "$(file_out)_$(odd_theta)_$(odd_tau)_$(seed).csv"
end

### The Gamma Test
function test_Gamma(file_out, numRuns, n_grid, seed, alpha, beta)
	o = GammaExp(seed, alpha, beta, maximum(n_grid))

	f = open("$(file_out)_Gamma_$(alpha)_$(beta)_$(seed).csv", "w")
	writecsv(f, ["Run" "n" "Method" "YVal" "thetaVal" "time" "tau0"])
	test_harness(f, numRuns, o, n_grid)
	close(f)
	return "$(file_out)_Gamma_$(alpha)_$(beta)_$(seed).csv"
end

### The Uniform Test
function test_Uniform(file_out, numRuns, n_grid, seed, a, b)
	o = UniformExp(seed, a, b, maximum(n_grid))

	f = open("$(file_out)_Uniform_$(a)_$(b)_$(seed).csv", "w")
	writecsv(f, ["Run" "n" "Method" "YVal" "thetaVal" "time" "tau0"])
	test_harness(f, numRuns, o, n_grid)
	close(f)
	return "$(file_out)_Uniform_$(a)_$(b)_$(seed).csv"
end

### The Beta Test
function test_Beta(file_out, numRuns, n_grid, seed, a, b)
	o = BetaExp(seed, a, b, maximum(n_grid))

	f = open("$(file_out)_Beta_$(a)_$(b)_$(seed).csv", "w")
	writecsv(f, ["Run" "n" "Method" "YVal" "thetaVal" "time" "tau0"])
	test_harness(f, numRuns, o, n_grid)
	close(f)
	return "$(file_out)_Beta_$(a)_$(b)_$(seed).csv"
end

### The Beta Test
function test_UniformLikelihood(numRuns, n_grid, seed, theta_high, width_low, width_high, lamp)
	o = UniformLikelihoodExp(seed, theta_high, 1., width_high, width_low, maximum(n_grid), lamp)
	f = open("UniformLikelihood_$(theta_high)_$(width_high)_$(width_low)_$(lamp)_$(seed).csv", "w")
	writecsv(f, ["Run" "n" "Method" "YVal" "thetaVal" "time" "tau0"])
	test_harness(f, numRuns, o, n_grid)
	close(f)
	return "UniformLikelihood_$(theta_high)_$(width_high)_$(width_low)_$(lamp)_$(seed).csv"
end

function test_CLTExp(file_out, numRuns, n_grid, seed, N, width_max)
	o = CLTExp(seed, width_max, maximum(n_grid), N)
	tag = "$(file_out)_CLTExp_$(N)_$(width_max)_$(seed).csv"
	f = open(tag, "w")
	writecsv(f, ["Run" "n" "Method" "YVal" "thetaVal" "time" "tau0"])
	test_harness(f, numRuns, o, n_grid)
	close(f)
	return tag
end

#########
n_grid = [2^i for i = 8:17]

#run some small examples to pre-compile for optimization
test_CLTExp("tempCLT", 5, [100, 150], 86, 10, 2)
test_Gaussian("temp_Gaussian2", 5, [100, 150], 87, 3, 5)
test_OddEven("temp_OddEven", 5,[100, 150], 87, 2, 2)
test_Gamma("temp_Gamma", 5, [100, 150], 87, 1., 1.)
test_Uniform("temp_Uniform", 5, [100, 150], 87, 1, 2)
test_Beta("temp_Beta", 5, [100, 150], 87, .5, .5)

# #The counterexample
# n_grid = [2^i for i = 8:20]
# test_UniformLikelihood(2, [100, 150], 8675309, 10., 15., 1., 14.1)

