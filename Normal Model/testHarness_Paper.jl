##Tests used in paper for the Small-Data KP
include("KPNormal.jl")
using Distributions, KP

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
function test_harness(f, numRuns, o, n_grid; includeReg=true)
	const n_max = maximum(n_grid)
	muhat = zeros(Float64, n_max)
	noise = zeros(Float64, n_max)

	#write a header
	writecsv(f, ["Run" "n" "Method" "thetaVal" "time" "tau0"])

	for iRun = 1:numRuns
		#generate the entire path up to n_max
		sim!(o, muhat)
		noise[:] = randn!(noise) ./ sqrt(o.vs)

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

			#Box with the optimized rate, i.e. h_n = n^-1/6 and scaling, altKernel
			h = n^-.16666
			tic()
			xs, vals, objs = x_stein_box(o.cs[1:n], muhat[1:n], o.vs[1:n], h, tau_step = .05)
			t = toc()
			thetaval = dot(o.thetas[1:n], xs)/n
			writecsv(f, [iRun n "BoxStein" thetaval t vals[indmax(objs)]])

			#Oracle Value
			tic()
			xs, vals, objs = best_x_tau(o.cs[1:n], muhat[1:n], o.vs[1:n], o.thetas[1:n])
			t = toc()
			thetaval = dot(o.thetas[1:n], xs)/n
			@assert abs(thetaval - maximum(objs)) <= 1e-5 "Weird Mismatch? \t $thetaval \t $(maximum(objs))"
			writecsv(f, [iRun n "OR" thetaval t vals[indmax(objs)]])

			if includeReg
				#Oracle Regularization
				tic()
				xs, Gamma_grid, objs = KP.x_l2reg_CV_warm(o.cs[1:n], muhat[1:n], o.vs[1:n], o.thetas[1:n], Gamma_max = 20.)
				t = toc()
				Gammahat = Gamma_grid[indmax(objs)]
				thetaval = dot(o.thetas[1:n], xs)/n
				writecsv(f, [iRun n "OracleReg" thetaval t Gammahat])

				#Our Stein Approach to Regularization
				tic()
				xs, Gamma_grid, objs = KP.x_stein_reg(o.cs[1:n], muhat[1:n], o.vs[1:n], Gamma_max = 20.)
				t = toc()
				Gammahat = Gamma_grid[indmax(objs)]
				thetaval = dot(o.thetas[1:n], xs)/n
				writecsv(f, [iRun n "SteinReg" thetaval t Gammahat])

				#Stein Appraoch with Bounds
				tic()
				xs, Gamma_grid, objs = KP.x_stein_reg(o.cs[1:n], muhat[1:n], o.vs[1:n], Gamma_min = 10., Gamma_max = 20.)
				t = toc()
				Gammahat = Gamma_grid[indmax(objs)]
				thetaval = dot(o.thetas[1:n], xs)/n
				writecsv(f, [iRun n "SteinRegBnded" thetaval t Gammahat])

				#RO heuristic for Gamma
				#eps = .1				
				tic()
				xs, lam = KP.x_l2reg(o.cs[1:n], muhat[1:n], o.vs[1:n], 2.563103)
				t = toc()
				thetaval = dot(o.thetas[1:n], xs)/n
				writecsv(f, [iRun n "RO_Eps_.1" thetaval t 2.563103])

				#eps = .05				
				tic()
				xs, lam = KP.x_l2reg(o.cs[1:n], muhat[1:n], o.vs[1:n], 3.289707)
				t = toc()
				thetaval = dot(o.thetas[1:n], xs)/n
				writecsv(f, [iRun n "RO_Eps_.05" thetaval t 3.289707])

				#eps = .01				
				tic()
				xs, lam = KP.x_l2reg(o.cs[1:n], muhat[1:n], o.vs[1:n], 4.652696)
				t = toc()
				thetaval = dot(o.thetas[1:n], xs)/n
				writecsv(f, [iRun n "RO_Eps_.01" thetaval t 4.652696])

				#Leave one out validation (LOO)
				tic()
				xs, Gamma_grid, objs = KP.x_LOO_reg(o.cs[1:n], muhat[1:n] + noise[1:n], muhat[1:n] - noise[1:n], o.vs[1:n])
				t = toc()
				thetaval = dot(o.thetas[1:n], xs)/n
				writecsv(f, [iRun n "LOO" thetaval t Gamma_grid[indmax(objs)]])
			end

			# #The primal stein approach
			# #use the optimized rate, i.e. h_n = n^-1/6
			# h = n^-.16666
			# tic()
			# xs, vals, objs = KP.x_stein_primal(o.cs[1:n], muhat[1:n], o.vs[1:n], h)
			# t = toc()
			# yval = dot(ys[1:n], xs)/n
			# thetaval = dot(o.thetas[1:n], xs)/n
			# writecsv(f, [iRun n "PrimalStein" yval thetaval t vals[indmax(objs)]])


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
function test_OddEven(file_out, numRuns, n_grid, seed, even_v; 
						even_theta=1., frac_fit=.1, includeReg=true)
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
	test_harness(f, numRuns, o, n_grid, includeReg=includeReg)
	close(f)
	return file_name
end

### The three parts experimental set-up
function test_threePart(file_out, numRuns, n_grid, seed, 
						theta_l, v_l, c_l, theta_m, v_m, c_m, theta_h, v_h, c_h, includeReg)
	srand(seed)
	const n_max = maximum(n_grid)
	cs = ones(n_max)
	thetas = ones(n_max)
	vs = ones(n_max)

	thetas[1:3:n_max] = theta_l
	vs[1:3:n_max] = v_l
	cs[1:3:n_max] = c_l

	thetas[2:3:n_max] = theta_m
	vs[2:3:n_max] = v_m
	cs[2:3:n_max] = c_m

	thetas[3:3:n_max] = theta_h
	vs[3:3:n_max] = v_h
	cs[3:3:n_max] = c_h

	o = DefaultExp(cs, thetas, vs)
	file_name = "$(file_out)_3part_$(theta_l)_$(v_l)_$(c_l)_$(theta_m)_$(v_m)_$(c_m)_$(theta_h)_$(v_h)_$(c_h)_$(seed).csv"
	f = open(file_name, "w")
	test_harness(f, numRuns, o, n_grid; includeReg=includeReg)
	close(f)
	return file_name
end


###  Reads in a theta/cs/vs specification
function test_ReadData(file_out, numRuns, n_grid, seed, param_path; 
						includeReg = true)
	srand(seed)
	const n_max = maximum(n_grid)
	dat, header = readcsv(param_path, header=true)

	#confirm that n_max works
	@assert n_max <= size(dat, 1) "Param file too short for n_max"
	cs = dat[1:n_max, 3]
	thetas = dat[1:n_max, 1]
	vs = dat[1:n_max, 2]

	#rescale cs so budget is sensible. 
	cs /= quantile(cs, .2)

	o = DefaultExp(cs, thetas, vs)
	file_name = "$(file_out)_$(seed).csv"
	f = open(file_name, "w")
	test_harness(f, numRuns, o, n_grid; includeReg=includeReg)
	close(f)
	return file_name	
end



### A by-hand testing multiple Gamma values and their convergence
# to the true oracle
function test_MultiGamma(file_out, numRuns, n_grid, seed, gamma_Grid;
						even_theta = 1., even_v = 2, frac_fit = .1)
	srand(seed)
	const n_max = maximum(n_grid)

	#build the sim object with the odd_even setup
	srand(seed)
	const n_max = maximum(n_grid)
	cs = 1./frac_fit * ones(n_max)
	thetas = zeros(n_max)
	thetas[2:2:n_max] = even_theta
	vs = ones(n_max)
	vs[2:2:n_max] = even_v
	o = DefaultExp(cs, thetas, vs)
	muhat = zeros(n_max)

	#output file
	file_name = "$(file_out)_multi_gamma_$(seed).csv"
	f = open(file_name, "w")
	#write a header
	writecsv(f, ["Run" "n" "Method" "val" "Gamma"])

	for iRun in 1:numRuns
		sim!(o, muhat)
		for n in n_grid
			for Gamma in gamma_Grid
				xhat, lam = x_l2reg(o.cs[1:n], muhat[1:n], o.vs[1:n], Gamma)
				const valOR = dot(o.thetas[1:n], xhat)/n
				const val = dot(muhat[1:n], xhat)/n - KP.reg_bias(o.cs[1:n], o.vs[1:n], muhat[1:n], lam, Gamma)

				writecsv(f, [iRun n "Oracle" valOR Gamma])			
				writecsv(f, [iRun n "Estimate" val Gamma])
			end
		end
		flush(f)
	end
	close(f)
	return file_name
end



########################################################
#########
n_grid = [2^i for i = 5:17]
#small run for pre-compilation
#test_Gaussian("./temp/temp_Gaussian", 5, [100, 150], 87, 3, 1, 3)
#test_OddEven("./temp/temp_OddEvenReg", 5, [100, 150], 8675309000, 2.1, includeReg=true)
test_ReadData("./temp/temp_PortExp", 5, [100, 150], 8675309, "./Results/param_portExp_Linear_.5.csv")
