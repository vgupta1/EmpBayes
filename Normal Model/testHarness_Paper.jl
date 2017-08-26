##Tests used in paper for the Small-Data KP
include("KPNormal.jl")
using Distributions, KP

##VG Sanity Checks
# 1) Does the emp-bayes policy converge for a variety of priors?
# 2) For a single run, is the regularized policy producing something useful?

#Default experiments fix a single theta, cs, vs across all simulation
type DefaultExp
	cs
	thetas
	vs
end

function sim!(o, muhat)
	muhat[:] = randn!(muhat) ./ sqrt(o.vs) + o.thetas
end


######################
####
#For increasing n, many simulations compare
#		#SAA
#		#FullInfo
#		#EB MM
#		#EB MLE
#		#SURE MSE
#		#OR MSE
#		#Diract Stein
#		#Box Stein
#		#tauOR
#Compare on metrics
#		#Value wrt theta
#		#Time to compute
#		#optimal value of tau0
function test_harness(f, numRuns, o, n_grid)
	const n_max = maximum(n_grid)
	muhat = zeros(Float64, n_max)

	#write a header
	writecsv(f, ["Run" "n" "Method" "thetaVal" "time" "tau0"])

	for iRun = 1:numRuns
		#generate the entire path up to n_max
		sim!(o, muhat)

		for n in n_grid
			#Compute performance of each method

			#SAA
			tic()
			xs = x(o.cs[1:n], muhat[1:n])
			t = toc()
			thetaval = dot(o.thetas[1:n], xs)/n
			writecsv(f, [iRun n "SAA" thetaval t 0.])

			#fullInfo val
			tic()
			xs = x(o.cs[1:n], o.thetas[1:n])
			t = toc()
			thetaval = dot(o.thetas[1:n], xs)/n
			writecsv(f, [iRun n "FullInfo" thetaval t 0.])

			#The "Bayes" value.. only possible bc we if we we are in bayesian
			if :tau0 in fieldnames(o)
				tic()
				xs = x(o.cs[1:n], shrink(muhat[1:n], o.vs[1:n], o.tau0))
				t = toc()
				thetaval = dot(o.thetas[1:n], xs)/n
				writecsv(f, [iRun n "Bayes" thetaval t o.tau0])
			end

			#Tau MLE
			tic()
			tauMLE, xs = x_MLE(o.cs[1:n], muhat[1:n], o.vs[1:n])
			t = toc()
			thetaval = dot(o.thetas[1:n], xs)/n
			writecsv(f, [iRun n "EB_MLE" thetaval t tauMLE])

			#Tau MM
			tic()
			tauMM, xs = x_MM(o.cs[1:n], muhat[1:n], o.vs[1:n])
			t = toc()
			thetaval = dot(o.thetas[1:n], xs)/n
			writecsv(f, [iRun n "EB_MM" thetaval t tauMM])

			#Oracle MSE
			tic()
			xs, tau_CV = x_OR_MSE(o.cs[1:n], muhat[1:n], o.thetas[1:n], o.vs[1:n])
			t = toc()
			thetaval = dot(o.thetas[1:n], xs)/n
			writecsv(f, [iRun n "OR_MSE" thetaval t tau_CV])

			#Sure MSE
			tic()
			xs, tau_CV = x_sure_MSE(o.cs[1:n], muhat[1:n], o.vs[1:n])
			t = toc()
			thetaval = dot(o.thetas[1:n], xs)/n
			writecsv(f, [iRun n "SURE_MSE" thetaval t tau_CV])


			#Dirac Stein
			tic()
			xs, vals, objs = x_stein_exact(o.cs[1:n], muhat[1:n], o.vs[1:n], o.thetas[1:n])
			t = toc()
			thetaval = dot(o.thetas[1:n], xs)/n
			writecsv(f, [iRun n "DiracStein" thetaval t vals[indmax(objs)]])

			#The stein approach with various kernels
			#Box with the optimized rate, i.e. h_n = n^-1/6 and scaling, altKernel
			h = n^-.16666
			tic()
			xs, vals, objs = x_stein_box(o.cs[1:n], muhat[1:n], o.vs[1:n], h, tau_step = .05)
			t = toc()
			thetaval = dot(o.thetas[1:n], xs)/n
			writecsv(f, [iRun n "BoxStein" thetaval t vals[indmax(objs)]])

			#weighted l2 regularization.  uses the oracle value for now
			tic()
			xs, mu = x_l2reg_CV(o.cs[1:n], muhat[1:n], o.vs[1:n], o.thetas[1:n])[1:2]		
			t = toc()
			thetaval = dot(o.thetas[1:n], xs)/n
			writecsv(f, [iRun n "OracleReg" thetaval t mu])

			# #The primal stein approach
			# #use the optimized rate, i.e. h_n = n^-1/6
			# h = n^-.16666
			# tic()
			# xs, vals, objs = KP.x_stein_primal(o.cs[1:n], muhat[1:n], o.vs[1:n], h)
			# t = toc()
			# yval = dot(ys[1:n], xs)/n
			# thetaval = dot(o.thetas[1:n], xs)/n
			# writecsv(f, [iRun n "PrimalStein" yval thetaval t vals[indmax(objs)]])

			#Oracle Value
			tic()
			xs, vals, objs = best_x_tau(o.cs[1:n], muhat[1:n], o.vs[1:n], o.thetas[1:n])
			t = toc()
			thetaval = dot(o.thetas[1:n], xs)/n
			@assert abs(thetaval - maximum(objs)) <= 1e-5 "Weird Mismatch? \t $thetaval \t $(maximum(objs))"
			writecsv(f, [iRun n "OR" thetaval t vals[indmax(objs)]])

		end
		flush(f)
	end
end

#######################

### Bayesian set-up with a non-zero mean
### costs and vs are random, but fixed over runs
function test_Gaussian(file_out, numRuns, n_grid, seed, tau0, mu0, avg_tau, 
						frac_fit = .1)
	#build the sim object
	const n_max = maximum(n_grid)
	srand(seed)
	cs = rand(n_max) * 2. / frac_fit   
	vs = 2. * avg_tau * rand(n_max)
	thetas = randn(n_max) / sqrt(tau0) + mu0
	o = DefaultExp(cs, thetas, vs)
	
	#run the testharness
	#output files
	file_name = "$(file_out)_Gaussian_$(tau0)_$(mu0)_$(avg_tau)_$(seed).csv"
	f = open(file_name, "w")
	test_harness(f, numRuns, o, n_grid)
	close(f)
	return file_name
end

### The odd even set-up.  
function test_OddEven(file_out, numRuns, n_grid, seed, even_v; even_theta=1., frac_fit=.1)
	#build the sim object
	srand(seed)
	const n_max = maximum(n_grid)
	cs = 1./frac_fit * ones(n_max)
	thetas = zeros(n_max)
	thetas[2:2:n_max] = even_theta
	vs = ones(n_max)
	vs[2:2:n_max] = even_v
	o = DefaultExp(cs, thetas, vs)

	#run the testharness
	#output files
	file_name = "$(file_out)_OddEven_$(even_theta)_$(even_v)_$(seed).csv"
	f = open(file_name, "w")
	test_harness(f, numRuns, o, n_grid)
	close(f)
	return file_name
end

# # ### The three parts experimental set-up
# # function test_threePart(file_out, numRuns, n_grid, seed, theta_l, v_l, theta_h, v_h)
# # 	o = threePartExp(seed, maximum(n_grid), theta_l, v_l, theta_h, v_h)
# # 	f = open("$(file_out)_$(theta_l)_$(v_l)_$(theta_h)_$(v_h)_$(seed).csv", "w")
# # 	test_harness(f, numRuns, o, n_grid)
# # 	close(f)
# # 	return "$(file_out)_$(theta_l)_$(v_l)_$(theta_h)_$(v_h)_$(seed).csv"
# # end

# ### The Gamma Test
# function test_Gamma(file_out, numRuns, n_grid, seed, alpha, beta)
# 	o = GammaExp(seed, alpha, beta, maximum(n_grid))

# 	f = open("$(file_out)_Gamma_$(alpha)_$(beta)_$(seed).csv", "w")
# 	test_harness(f, numRuns, o, n_grid)
# 	close(f)
# 	return "$(file_out)_Gamma_$(alpha)_$(beta)_$(seed).csv"
# end

# ### The Uniform Test
# function test_Uniform(file_out, numRuns, n_grid, seed, a, b)
# 	o = UniformExp(seed, a, b, maximum(n_grid))

# 	f = open("$(file_out)_Uniform_$(a)_$(b)_$(seed).csv", "w")
# 	test_harness(f, numRuns, o, n_grid)
# 	close(f)
# 	return "$(file_out)_Uniform_$(a)_$(b)_$(seed).csv"
# end

# ### The Beta Test
# function test_Beta(file_out, numRuns, n_grid, seed, a, b)
# 	o = BetaExp(seed, a, b, maximum(n_grid))

# 	f = open("$(file_out)_Beta_$(a)_$(b)_$(seed).csv", "w")
# 	test_harness(f, numRuns, o, n_grid)
# 	close(f)
# 	return "$(file_out)_Beta_$(a)_$(b)_$(seed).csv"
# end

# ### The Beta Test
# function test_UniformLikelihood(numRuns, n_grid, seed, theta_high, width_low, width_high, lamp)
# 	o = UniformLikelihoodExp(seed, theta_high, 1., width_high, width_low, maximum(n_grid), lamp)
# 	f = open("UniformLikelihood_$(theta_high)_$(width_high)_$(width_low)_$(lamp)_$(seed).csv", "w")
# 	test_harness(f, numRuns, o, n_grid)
# 	close(f)
# 	return "UniformLikelihood_$(theta_high)_$(width_high)_$(width_low)_$(lamp)_$(seed).csv"
# end

# function test_CLTExp(file_out, numRuns, n_grid, seed, N, width_max)
# 	o = CLTExp(seed, width_max, maximum(n_grid), N)
# 	tag = "$(file_out)CLTExp_$(N)_$(width_max)_$(seed).csv"
# 	f = open(tag, "w")
# 	test_harness(f, numRuns, o, n_grid)
# 	close(f)
# 	return tag
# end


# #cs and vs fixed across all runs
# function OddEvenExp(seed::Integer, n_max::Integer, odd_theta, odd_tau, frac_fit=.1)
# 	srand(seed)
# 	cs = 1./frac_fit * ones(n_max)

# 	thetas = ones(n_max)
# 	thetas[1:2:n_max] = odd_theta
# 	vs = ones(n_max)
# 	vs[1:2:n_max] = odd_tau

# 	noise = zeros(n_max)
# 	DefaultExp(cs, thetas, vs, noise)
# end

# function threePartExp(seed::Integer, n_max::Integer, theta_l, v_l, theta_h, v_h, frac_fit = .1)
# 	srand(seed)
# 	cs = 1./frac_fit * ones(n_max)
# 	thetas = ones(n_max)
# 	vs = ones(n_max)
# 	thetas[1:3:n_max] = theta_l
# 	vs[1:3:n_max] = v_l
# 	thetas[3:3:n_max] = theta_h
# 	vs[3:3:n_max] = v_h

# 	noise = zeros(n_max)
# 	DefaultExp(cs, thetas, vs, noise)
# end

# #Following simulates N rvs according to the given distributions
# #AS N increases, this should approximately have gaussian likelihood
# #dists should be mean zero.
# type CLTExp
# 	cs
# 	vs
# 	thetas
# 	dists  #should be a mean zero r.v. 
# 	N
# end

#thetas /vs follow the three pt experiment
#width calculated accordingly
#Obsservations are uniform, centered on theta with width
# function CLTExp(seed::Integer, width_max, n_max, N; frac_fit = .1)
# 	srand(seed)
# 	@assert rem(N, 2) == 0
# 	cs = 1./frac_fit * ones(n_max)
# 	thetas = ones(n_max)
# 	vs = ones(n_max)
# 	thetas[1:3:n_max] = .01
# 	vs[1:3:n_max] = .01
# 	thetas[3:3:n_max] = 1.5
# 	vs[3:3:n_max] = .1

# 	widths = sqrt(12 ./ vs)

# 	# widths = rand(n_max) * width_max + .1
# 	# vs = 12 ./ widths.^2
# 	# thetas = randn(n_max) ./ sqrt(tau0)

# 	dists = Array(Distributions.Uniform, n_max)
# 	for ix = 1:n_max
# 		dists[ix] = Uniform(-widths[ix]/2, widths[ix]/2 )
# 	end
# 	CLTExp(cs, vs, thetas, dists, N)
# end

# function sim!(o::CLTExp, xs, ys, muhat)	
# 	const n_max = length(o.cs)
# 	zetas = zeros(o.N)
# 	const half_N = round(Int, o.N/2)
# 	for ix = 1:n_max
# 		rand!(o.dists[ix], zetas)
# 		xs[ix] = mean(zetas[1:half_N]) * sqrt(half_N) + o.thetas[ix]
# 		ys[ix] = mean(zetas[(half_N + 1):o.N]) * sqrt(half_N) + o.thetas[ix]
# 		muhat[ix] = mean(zetas) * sqrt(o.N) + o.thetas[ix]
# 	end
# end

# function NormalBayesExp(seed::Integer, tau0, n_max; mu0 = 0., frac_fit=.1, avg_tau=tau0)
# 	srand(seed)
# 	cs = rand(n_max) * 2./frac_fit   
# 	vs = 2. * avg_tau * rand(n_max)
# 	NormalBayesExp(cs, vs, zeros(n_max), zeros(n_max), tau0, mu0)
# end

# #cs vs remain fixed but thetas change
# function sim!(o::NormalBayesExp, x, y, z)
# 	randn!(o.thetas) 
# 	o.thetas /= sqrt(o.tau0)
# 	o.thetas += o.mu0
# 	simData!(o, x, y, z)
# end





#######
# A by-hand test for bandwidths and reg parameters
# o is assumed pre-initialized
# function test_bandwidth(file_out, numRuns, o)
# 	const n = length(o.thetas)
# 	xs = zeros(Float64, n)
# 	ys = zeros(Float64, n)
# 	muhat = zeros(Float64, n)

# 	f = open(file_out, "w")
# 	writecsv(f, ["Run" "n" "Method" "bandwidth" "thetaVal" "tau0"])
# 	bandwidths = linspace(-1, -.001, 30)
# 	#regs = linspace(.05, 5, 15)

# 	for iRun = 1:numRuns
# 		sim!(o, xs, ys, muhat)

# 		for h in bandwidths
# 			xs, vals, objs = KP.x_stein_box(o.cs, muhat, o.vs, n^h, KP.box, tau_step = .1)
# 			thetaval = dot(o.thetas, xs)/n
# 			writecsv(f, [iRun n "Box" h thetaval vals[indmax(objs)]])

# 			xs, vals, objs = KP.x_stein_box(o.cs, muhat, o.vs, n^h, KP.gauss, tau_step = .1)
# 			thetaval = dot(o.thetas, xs)/n
# 			writecsv(f, [iRun n "Gauss" h thetaval vals[indmax(objs)]])

# 			xs, vals, objs = KP.x_stein_box(o.cs, muhat, o.vs, n^h, sinc, tau_step = .1)
# 			thetaval = dot(o.thetas, xs)/n
# 			writecsv(f, [iRun n "Sinc" h thetaval vals[indmax(objs)]])
# 		end
# 	end
# 	close(f)
# end

#########
n_grid = [2^i for i = 8:17]
test_Gaussian("./temp/temp_Gaussian", 5, [100, 150], 87, 3, 1, 3)

#run some small examples to pre-compile for optimization
# test_threePart("./temp/tempThreePart", 5, [100, 150], 876, .01, .01, 1.5, .1)
# test_CLTExp("./temp/tempCLT", 5, [100, 150], 86, 10, 2)
# test_OddEven("./temp/temp_OddEven", 5,[100, 150], 87, 2, 2)


# test_Gamma("./temp/temp_Gamma", 5, [100, 150], 87, 1., 1.)
# test_Uniform("./temp/temp_Uniform", 5, [100, 150], 87, 1, 2)
# test_Beta("./temp/temp_Beta", 5, [100, 150], 87, .5, .5)
# test_bandwidth("./temp/tempbandwidth", 10, NormalBayesExp(8675309, 3, 100))

# #The counterexample
# n_grid = [2^i for i = 8:20]
# test_UniformLikelihood(2, [100, 150], 8675309, 10., 15., 1., 14.1)
