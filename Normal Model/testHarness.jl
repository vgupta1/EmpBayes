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

#Default experiments fix a single theta, cs, vs across all simulation
type DefaultExp
	cs
	thetas
	vs
	noise
end

#cs and vs fixed across all runs
function OddEvenExp(seed::Integer, n_max::Integer, odd_theta, odd_tau, frac_fit=.1)
	srand(seed)
	cs = 1./frac_fit * ones(n_max)

	thetas = ones(n_max)
	thetas[1:2:n_max] = odd_theta
	vs = ones(n_max)
	vs[1:2:n_max] = odd_tau

	noise = zeros(n_max)
	DefaultExp(cs, thetas, vs, noise)
end

function threePartExp(seed::Integer, n_max::Integer, theta_l, v_l, theta_h, v_h, frac_fit = .1)
	srand(seed)
	cs = 1./frac_fit * ones(n_max)
	thetas = ones(n_max)
	vs = ones(n_max)
	thetas[1:3:n_max] = theta_l
	vs[1:3:n_max] = v_l
	thetas[3:3:n_max] = theta_h
	vs[3:3:n_max] = v_h

	noise = zeros(n_max)
	DefaultExp(cs, thetas, vs, noise)
end

#Following simulates N rvs according to the given distributions
#AS N increases, this should approximately have gaussian likelihood
#dists should be mean zero.
type CLTExp
	cs
	vs
	thetas
	dists
	N
end

#thetas are initially generated as gaussians
#widths are randomly genrated between [0, width_max]
#Likelihoods are uniform, centered on theta with width
function CLTExp(seed::Integer, width_max, n_max, N; tau0=2, frac_fit = .1)
	srand(seed)
	@assert rem(N, 2) == 0
	cs = rand(n_max) * 2./frac_fit   

	widths = rand(n_max) * width_max + .1
	vs = 12 ./ widths.^2
	thetas = randn(n_max) ./ sqrt(tau0)

	dists = Array(Distributions.Uniform, n_max)
	for ix = 1:n_max
		dists[ix] = Uniform(-widths[ix]/2, widths[ix]/2)
	end
	CLTExp(cs, vs, thetas, dists, N)
end

function sim!(o::CLTExp, xs, ys, zs)	
	const n_max = length(o.cs)
	zetas = zeros(o.N)
	const half_N = round(Int, o.N/2)
	for ix = 1:n_max
		rand!(o.dists[ix], zetas)
		xs[ix] = mean(zetas[1:half_N]) * sqrt(half_N) + o.thetas[ix]
		ys[ix] = mean(zetas[(half_N + 1):o.N]) * sqrt(half_N) + o.thetas[ix]
		zs[ix] = mean(zetas) * sqrt(o.N) + o.thetas[ix]
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

function NormalBayesExp(seed::Integer, tau0, n_max; mu0 = 0., frac_fit=.1, avg_tau=tau0)
	srand(seed)
	cs = rand(n_max) * 2./frac_fit   
	vs = 2. * avg_tau * rand(n_max)
	NormalBayesExp(cs, vs, zeros(n_max), zeros(n_max), tau0, mu0)
end

#cs vs remain fixed but thetas change
function sim!(o::NormalBayesExp, x, y, z)
	randn!(o.thetas) 
	o.thetas /= sqrt(o.tau0)
	o.thetas += o.mu0
	simData!(o, x, y, z)
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

function UniformExp(seed::Integer, a, b, n_max, frac_fit = .1, avg_vs=2)
	srand(seed)
	cs = rand(n_max) * 2./frac_fit   
	vs = 2. * avg_vs * rand(n_max)
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
#		#SAA/ Naive Z
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

	#write a header
	writecsv(f, ["Run" "n" "Method" "thetaVal" "time" "tau0"])

	for iRun = 1:numRuns
		#generate the entire path up to n_max
		sim!(o, xs, ys, zs)

		for n in n_grid
			#Compute performance of each method

			#Naive Z
			#Same as the ZZ/SAA
			tic()
			qs = q(o.cs[1:n], zs[1:n])
			t = toc()
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "SAA" thetaval t 0.])

			#fullInfo val
			tic()
			qs = q(o.cs[1:n], o.thetas[1:n])
			t = toc()
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "FullInfo" thetaval t 0.])

			#The "Bayes" value.. only possible bc we if we we are in bayesian
			if :tau0 in fieldnames(o)
				tic()
				qs = q(o.cs[1:n], shrink(zs[1:n], o.vs[1:n], o.tau0))
				t = toc()
				thetaval = dot(o.thetas[1:n], qs)/n
				writecsv(f, [iRun n "Bayes" thetaval t o.tau0])
			end

			#Tau MLE
			tic()
			tauMLE, qs = q_MLE(o.cs[1:n], zs[1:n], o.vs[1:n])
			t = toc()
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "MLE" thetaval t tauMLE])

			#Tau MM
			tic()
			tauMM, qs = q_MM(o.cs[1:n], zs[1:n], o.vs[1:n])
			t = toc()
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "MM" thetaval t tauMM])

			#The half-sample heuristic
			tic()
			qs, vals, objs = best_q_tau(o.cs[1:n], xs[1:n], o.vs[1:n], ys[1:n])
			t = toc()
			
			tau_RS = vals[indmax(objs)]/2
			qs = q(o.cs[1:n], shrink(zs[1:n], o.vs[1:n], tau_RS))
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "Rescaled" thetaval t tau_RS])

			#The idealized stein approach
			tic()
			qs, vals, objs = KP.stein_q_tau_exact(o.cs[1:n], zs[1:n], o.vs[1:n], o.thetas[1:n])
			t = toc()
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "ExactStein" thetaval t vals[indmax(objs)]])

			#The stein approach with various kernels
			#Box with the optimized rate, i.e. h_n = n^-1/6
			h = n^-.16666
			tic()
			qs, vals, objs = KP.stein_q_tau_impulse(o.cs[1:n], zs[1:n], o.vs[1:n], h, KP.box, tau_step = .05)
			t = toc()
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "Box" thetaval t vals[indmax(objs)]])

			#Box with the optimized rate, i.e. h_n = n^-1/6 and scaling
			h = n^-.16666
			tic()
			qs, vals, objs = KP.stein_q_tau_impulse(o.cs[1:n], zs[1:n], o.vs[1:n], h, KP.box, tau_step = .05, scale_h=true)
			t = toc()
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "BoxS" thetaval t vals[indmax(objs)]])

			#gauss with MISE rate n^-1/5
			h = n^-.2
			tic()
			qs, vals, objs = KP.stein_q_tau_impulse(o.cs[1:n], zs[1:n], o.vs[1:n], h, KP.gauss, tau_step = .05)
			t = toc()
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "Gauss" thetaval t vals[indmax(objs)]])

			#gauss with MISE rate n^-1/5 with scaling
			h = n^-.2
			tic()
			qs, vals, objs = KP.stein_q_tau_impulse(o.cs[1:n], zs[1:n], o.vs[1:n], h, KP.gauss, tau_step = .05, scale_h = true)
			t = toc()
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "GaussS" thetaval t vals[indmax(objs)]])

			#gauss with optimized  rate n^-1/6
			h = n^-.16666
			tic()
			qs, vals, objs = KP.stein_q_tau_impulse(o.cs[1:n], zs[1:n], o.vs[1:n], h, KP.gauss, tau_step = .05)
			t = toc()
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "GaussO" thetaval t vals[indmax(objs)]])

			#gauss with optimized rate n^-1/6 with scaling
			h = n^-.16666
			tic()
			qs, vals, objs = KP.stein_q_tau_impulse(o.cs[1:n], zs[1:n], o.vs[1:n], h, KP.gauss, tau_step = .05, scale_h = true)
			t = toc()
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "GaussOS" thetaval t vals[indmax(objs)]])

			#sinc with 1/log(n)
			h = 1/log(n + 1)
			tic()
			qs, vals, objs = KP.stein_q_tau_impulse(o.cs[1:n], zs[1:n], o.vs[1:n], h, sinc, tau_step = .1)
			t = toc()
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "Sinc" thetaval t vals[indmax(objs)]])

			#sinc with 1/log(n) with scaling
			h = 1/log(n + 1)
			tic()
			qs, vals, objs = KP.stein_q_tau_impulse(o.cs[1:n], zs[1:n], o.vs[1:n], h, sinc, tau_step = .1, scale_h=true)
			t = toc()
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "SincS" thetaval t vals[indmax(objs)]])

			#Shrinkage via Cross-Val ("Ridge")
			tic()
			qs, tau_CV = q_CVShrink(o.cs[1:n], xs[1:n], ys[1:n], o.vs[1:n])
			t = toc()
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "CV_Shrink" thetaval t tau_CV])

			#weighted l2 regularization.  uses the oracle value for now
			tic()
			qs, mu = q_l2reg_CV(o.cs[1:n], zs[1:n], o.vs[1:n], o.thetas[1:n])[1:2]		
			t = toc()
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "OracleReg" thetaval t mu])

			#weighted l2 regularization.  uses Cross-val to find value and then resolves
			tic()
			qs, mu = q_l2reg_CV(o.cs[1:n], xs[1:n], o.vs[1:n], ys[1:n])[1:2]		
			t = toc()
			thetaval = dot(o.thetas[1:n], q_l2reg(o.cs[1:n], zs[1:n], o.vs[1:n], mu)[1])/n
			writecsv(f, [iRun n "RegCV" thetaval t mu])

			#Plug-in SURE estimation for L2
			tic()
			qs, tau = q_sure(o.cs[1:n], zs[1:n], o.vs[1:n])
			t = toc()
			thetaval = dot(o.thetas[1:n], qs)/n
			writecsv(f, [iRun n "SURE" thetaval t tau])

			# #The primal stein approach
			# #use the optimized rate, i.e. h_n = n^-1/6
			# h = n^-.16666
			# tic()
			# qs, vals, objs = KP.stein_q_tau_primal(o.cs[1:n], zs[1:n], o.vs[1:n], h)
			# t = toc()
			# yval = dot(ys[1:n], qs)/n
			# thetaval = dot(o.thetas[1:n], qs)/n
			# writecsv(f, [iRun n "PrimalStein" yval thetaval t vals[indmax(objs)]])

			#Oracle Value
			tic()
			qs, vals, objs = best_q_tau(o.cs[1:n], zs[1:n], o.vs[1:n], o.thetas[1:n])
			t = toc()
			thetaval = dot(o.thetas[1:n], qs)/n
			if abs(thetaval - maximum(objs)) > 1e-5
				f2 = open("Mismatch_log.csv", "w")
				writecsv(f2, [n iRun])
				writecsv(f2, xs')
				writecsv(f2, ys')
				writecsv(f2, zs')
				writecsv(f2, o.vs[1:n]')
				writecsv(f2, o.thetas[1:n]')
				writecsv(f2, o.cs[1:n]')
				close(f2)
				exit()
			end

#			@assert abs(thetaval - maximum(objs)) <= 1e-5 "Weird Mismatch? \t $thetaval \t $(maximum(objs))"
			writecsv(f, [iRun n "OracleZ" thetaval t vals[indmax(objs)]])

		end
		flush(f)
	end
end

#######################

### Bayesian set-up with a non-zero mean
### costs and vs are random, but fixed over runs
function test_Gaussian(file_out, numRuns, n_grid, seed, tau0, mu0, avg_tau)
	#build the sim object
	o = NormalBayesExp(seed, tau0, maximum(n_grid), mu0=mu0, avg_tau=avg_tau)

	#run the testharness
	#output files
	f = open("$(file_out)_$(tau0)_$(mu0)_$(avg_tau)_$(seed).csv", "w")
	test_harness(f, numRuns, o, n_grid)
	close(f)
	return "$(file_out)_$(tau0)_$(mu0)_$(avg_tau)_$(seed).csv"
end

### The odd even set-up.  
function test_OddEven(file_out, numRuns, n_grid, seed, odd_theta, odd_tau)
	#build the sim object
	o = OddEvenExp(seed, maximum(n_grid), odd_theta, odd_tau)

	#run the testharness
	#output files
	f = open("$(file_out)_$(odd_theta)_$(odd_tau)_$(seed).csv", "w")
	test_harness(f, numRuns, o, n_grid)
	close(f)
	return "$(file_out)_$(odd_theta)_$(odd_tau)_$(seed).csv"
end

### The three parts experimental set-up
function test_threePart(file_out, numRuns, n_grid, seed, theta_l, v_l, theta_h, v_h)
	o = threePartExp(seed, maximum(n_grid), theta_l, v_l, theta_h, v_h)
	f = open("$(file_out)_$(theta_l)_$(v_l)_$(theta_h)_$(v_h)_$(seed).csv", "w")
	test_harness(f, numRuns, o, n_grid)
	close(f)
	return "$(file_out)_$(theta_l)_$(v_l)_$(theta_h)_$(v_h)_$(seed).csv"
end

### The Gamma Test
function test_Gamma(file_out, numRuns, n_grid, seed, alpha, beta)
	o = GammaExp(seed, alpha, beta, maximum(n_grid))

	f = open("$(file_out)_Gamma_$(alpha)_$(beta)_$(seed).csv", "w")
	test_harness(f, numRuns, o, n_grid)
	close(f)
	return "$(file_out)_Gamma_$(alpha)_$(beta)_$(seed).csv"
end

### The Uniform Test
function test_Uniform(file_out, numRuns, n_grid, seed, a, b)
	o = UniformExp(seed, a, b, maximum(n_grid))

	f = open("$(file_out)_Uniform_$(a)_$(b)_$(seed).csv", "w")
	test_harness(f, numRuns, o, n_grid)
	close(f)
	return "$(file_out)_Uniform_$(a)_$(b)_$(seed).csv"
end

### The Beta Test
function test_Beta(file_out, numRuns, n_grid, seed, a, b)
	o = BetaExp(seed, a, b, maximum(n_grid))

	f = open("$(file_out)_Beta_$(a)_$(b)_$(seed).csv", "w")
	test_harness(f, numRuns, o, n_grid)
	close(f)
	return "$(file_out)_Beta_$(a)_$(b)_$(seed).csv"
end

### The Beta Test
function test_UniformLikelihood(numRuns, n_grid, seed, theta_high, width_low, width_high, lamp)
	o = UniformLikelihoodExp(seed, theta_high, 1., width_high, width_low, maximum(n_grid), lamp)
	f = open("UniformLikelihood_$(theta_high)_$(width_high)_$(width_low)_$(lamp)_$(seed).csv", "w")
	test_harness(f, numRuns, o, n_grid)
	close(f)
	return "UniformLikelihood_$(theta_high)_$(width_high)_$(width_low)_$(lamp)_$(seed).csv"
end

function test_CLTExp(file_out, numRuns, n_grid, seed, N, width_max)
	o = CLTExp(seed, width_max, maximum(n_grid), N)
	tag = "$(file_out)_CLTExp_$(N)_$(width_max)_$(seed).csv"
	f = open(tag, "w")
	test_harness(f, numRuns, o, n_grid)
	close(f)
	return tag
end


#######
# A by-hand test for bandwidths and reg parameters
# o is assumed pre-initialized
function test_bandwidth(file_out, numRuns, o)
	const n = length(o.thetas)
	xs = zeros(Float64, n)
	ys = zeros(Float64, n)
	zs = zeros(Float64, n)

	f = open(file_out, "w")
	writecsv(f, ["Run" "n" "Method" "bandwidth" "thetaVal" "tau0"])
	bandwidths = linspace(-1, -.001, 30)
	#regs = linspace(.05, 5, 15)

	for iRun = 1:numRuns
		sim!(o, xs, ys, zs)

		for h in bandwidths
			qs, vals, objs = KP.stein_q_tau_impulse(o.cs, zs, o.vs, n^h, KP.box, tau_step = .1)
			thetaval = dot(o.thetas, qs)/n
			writecsv(f, [iRun n "Box" h thetaval vals[indmax(objs)]])

			qs, vals, objs = KP.stein_q_tau_impulse(o.cs, zs, o.vs, n^h, KP.gauss, tau_step = .1)
			thetaval = dot(o.thetas, qs)/n
			writecsv(f, [iRun n "Gauss" h thetaval vals[indmax(objs)]])

			qs, vals, objs = KP.stein_q_tau_impulse(o.cs, zs, o.vs, n^h, sinc, tau_step = .1)
			thetaval = dot(o.thetas, qs)/n
			writecsv(f, [iRun n "Sinc" h thetaval vals[indmax(objs)]])
		end
	end
	close(f)
end

#########
n_grid = [2^i for i = 8:17]

#run some small examples to pre-compile for optimization
test_threePart("./temp/tempThreePart", 5, [100, 150], 876, .01, .01, 1.5, .1)
test_CLTExp("./temp/tempCLT", 5, [100, 150], 86, 10, 2)
# test_Gaussian("./temp/temp_Gaussian", 5, [100, 150], 87, 3, 1, 3)
# test_OddEven("./temp/temp_OddEven", 5,[100, 150], 87, 2, 2)


# test_Gamma("./temp/temp_Gamma", 5, [100, 150], 87, 1., 1.)
# test_Uniform("./temp/temp_Uniform", 5, [100, 150], 87, 1, 2)
# test_Beta("./temp/temp_Beta", 5, [100, 150], 87, .5, .5)
# test_bandwidth("./temp/tempbandwidth", 10, NormalBayesExp(8675309, 3, 100))

# #The counterexample
# n_grid = [2^i for i = 8:20]
# test_UniformLikelihood(2, [100, 150], 8675309, 10., 15., 1., 14.1)
