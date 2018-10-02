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
	muhat[:] = randn!(muhat) ./ sqrt.(o.vs) .+ o.thetas
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
function test_harness(f, numRuns, o, n_grid; includeReg=true, Gamma_min=1., Gamma_max=20.)
	const n_max = maximum(n_grid)
	muhat = zeros(Float64, n_max)
	noise = zeros(Float64, n_max)
	xs  = zeros(Float64, n_max)
	lam_t = 0.

	#write a header
	writecsv(f, ["Run" "n" "Method" "thetaVal" "time" "tau0"])

	for iRun = 1:numRuns
		#generate the entire path up to n_max
		sim!(o, muhat)
		noise[:] = randn!(noise) ./ sqrt.(o.vs)

		for n in n_grid
			#Take Views on Everything to speed up garbage collection
			x_t = view(xs, 1:n)
			vs  = view(o.vs, 1:n)
			cs  = view(o.cs, 1:n)
			thetas = view(o.thetas, 1:n)
			muhat_t = view(muhat, 1:n)
			noise_t = view(noise, 1:n)

			#Compute performance of each method

			#SAA
			tic()
			x_t[:] = x(cs, muhat_t)
			t = toc()
			thetaval = dot(thetas, x_t)/n
			writecsv(f, [iRun n "SAA" thetaval t 0.])

			#fullInfo val
			tic()
			x_t[:] = x(cs, thetas)
			t = toc()
			thetaval = dot(thetas, x_t)/n
			writecsv(f, [iRun n "FullInfo" thetaval t 0.])

			# #The "Bayes" value.. only possible bc we if we we are in bayesian
			# if :tau0 in fieldnames(o)
			# 	tic()
			# 	xs = x(cs, shrink(muhat_t, vs, o.tau0))
			# 	t = toc()
			# 	thetaval = dot(thetas, xs)/n
			# 	writecsv(f, [iRun n "Bayes" thetaval t o.tau0])
			# end

			#Tau MLE
			tic()
			tauMLE, x_t[:] = x_MLE(cs, muhat_t, vs)
			t = toc()
			thetaval = dot(thetas, x_t)/n
			writecsv(f, [iRun n "EB_MLE" thetaval t tauMLE])

			#Tau MM
			tic()
			tauMM, x_t[:] = x_MM(cs, muhat_t, vs)
			t = toc()
			thetaval = dot(thetas, x_t)/n
			writecsv(f, [iRun n "EB_MM" thetaval t tauMM])

			#Oracle MSE
			tic()
			x_t[:], tau_CV = x_OR_MSE(cs, muhat_t, thetas, vs)
			t = toc()
			thetaval = dot(thetas, x_t)/n
			writecsv(f, [iRun n "OR_MSE" thetaval t tau_CV])

			#Sure MSE
			tic()
			x_t[:], tau_CV = x_sure_MSE(cs, muhat_t, vs)
			t = toc()
			thetaval = dot(thetas, x_t)/n
			writecsv(f, [iRun n "SURE_MSE" thetaval t tau_CV])

			#Dirac Stein
			tic()
			x_t[:], vals, objs = x_stein_exact(cs, muhat_t, vs, thetas)
			t = toc()
			thetaval = dot(thetas, x_t)/n
			writecsv(f, [iRun n "DiracStein" thetaval t vals[indmax(objs)]])

			#Box with the optimized rate, i.e. h_n = n^-1/6 and scaling, altKernel
			h = n^-.16666
			tic()
			x_t[:], vals, objs = x_stein_box(cs, muhat_t, vs, h, tau_step = .05)
			t = toc()
			thetaval = dot(thetas, x_t)/n
			writecsv(f, [iRun n "BoxStein" thetaval t vals[indmax(objs)]])

			#Oracle Value
			tic()
			x_t[:], vals, objs = best_x_tau(cs, muhat_t, vs, thetas)
			t = toc()
			thetaval = dot(thetas, x_t)/n
			writecsv(f, [iRun n "OR" thetaval t vals[indmax(objs)]])

			if includeReg
				#Oracle Regularization
				tic()
				x_t[:], Gamma_grid, objs = KP.x_l2reg_CV(cs, muhat_t, vs, thetas, 
															Gamma_min=Gamma_min, Gamma_max=Gamma_max)
				t = toc()
				Gammahat = Gamma_grid[indmax(objs)]
				thetaval = dot(thetas, x_t)/n
				writecsv(f, [iRun n "OracleReg" thetaval t Gammahat])

				#From the same values, extract oracles for Other pairs
				#For Gamma_min = 5
				ind_min = findfirst(Gamma_grid .>= 5.0)
				Gammahat = Gamma_grid[ind_min:end][indmax(objs[ind_min:end])]
				lam_t = KP.x_l2reg2!(cs, muhat_t, vs, Gammahat, x_t)
				thetaval = dot(thetas, x_t)/n
				writecsv(f, [iRun n "OracleReg_5" thetaval t Gammahat])

				#For Gamma_min = 10
				#VG Needs fixing to match new signatures if you want to uncomment
				# ind_min = findfirst(Gamma_grid .>= 10.0)
				# Gammahat = Gamma_grid[ind_min:end][indmax(objs[ind_min:end])]
				# xs = KP.x_l2reg(cs, muhat_t, vs, Gammahat)[1]
				# thetaval = dot(thetas, xs)/n
				# writecsv(f, [iRun n "OracleReg_10" thetaval t Gammahat])

				#Our Stein Approach to Regularization
				tic()
				x_t[:], Gamma_grid, objs = KP.x_stein_reg(cs, muhat_t, vs, 
											Gamma_min=Gamma_min, Gamma_max=Gamma_max)
				t = toc()
				Gammahat = Gamma_grid[indmax(objs)]
				thetaval = dot(thetas, x_t)/n
				writecsv(f, [iRun n "SteinReg" thetaval t Gammahat])

				#Again, from same values, extract values for gamma_min
				ind_min = findfirst(Gamma_grid .>= 5.0)
				Gammahat = Gamma_grid[ind_min:end][indmax(objs[ind_min:end])]
				lam_t = KP.x_l2reg2!(cs, muhat_t, vs, Gammahat, x_t)
				thetaval = dot(thetas, x_t)/n
				writecsv(f, [iRun n "SteinReg_5" thetaval t Gammahat])

				#VG these need fixing to match new signatures if you want to uncomment
				# ind_min = findfirst(Gamma_grid .>= 10.0)
				# Gammahat = Gamma_grid[ind_min:end][indmax(objs[ind_min:end])]
				# xs = KP.x_l2reg(cs, muhat_t, vs, Gammahat)[1]
				# thetaval = dot(thetas, xs)/n
				# writecsv(f, [iRun n "SteinReg_10" thetaval t Gammahat])

				#NEW RO method with new threshold 
				tic()
				thresh = sqrt(2*log(1/.1))
				x_t[:] = KP.x_robFW(cs, muhat_t, vs, thresh, TOL=1e-4)
				t = toc()
				thetaval = dot(thetas, x_t)/n
				writecsv(f, [iRun n "FWRO_Eps_.1" thetaval t thresh])

				tic()
				thresh = sqrt(2*log(1/.05))
				x_t[:] = KP.x_robFW(cs, muhat_t, vs, thresh, TOL=1e-4)
				t = toc()
				thetaval = dot(thetas, x_t)/n
				writecsv(f, [iRun n "FWRO_Eps_.05" thetaval t thresh])

				tic()
				thresh = sqrt(2*log(1/.01))
				x_t[:] = KP.x_robFW(cs, muhat_t, vs, thresh, TOL=1e-4)
				t = toc()
				thetaval = dot(thetas, x_t)/n
				writecsv(f, [iRun n "FWRO_Eps_.01" thetaval t thresh])


				# ####
				# #This section was used for the intial paper submission Feb2017
				# # The thresholds have been updated.  If want to rerun, consider using the FW implementation too.  
				# #RO heuristic for Gamma
				# #eps = .05				
				# tic()
				# xs, lam = KP.x_rob(cs, muhat_t, vs, 1.6448536269514717)
				# t = toc()
				# thetaval = dot(thetas, xs)/n
				# writecsv(f, [iRun n "RO_Eps_.05" thetaval t 1.6448536269514717])

				# #eps = .01				
				# tic()
				# xs, lam = KP.x_rob(cs, muhat_t, vs, 2.326347874040845)
				# t = toc()
				# thetaval = dot(thetas, xs)/n
				# writecsv(f, [iRun n "RO_Eps_.01" thetaval t 2.326347874040845])

				#Leave one out validation (LOO)
				tic()
				x_t[:], Gamma_grid, objs = KP.x_LOO_reg(cs, muhat_t + noise_t, muhat_t - noise_t, vs, 
													Gamma_min=Gamma_min, Gamma_max=Gamma_max)
				t = toc()
				thetaval = dot(thetas, x_t)/n
				GammaLOO = Gamma_grid[indmax(objs)]
				writecsv(f, [iRun n "LOO" thetaval t GammaLOO])

				#Use same values to extract for other gamma_min
				ind_min = findfirst(Gamma_grid .>= 5.0)
				GammaLOO = Gamma_grid[ind_min:end][indmax(objs[ind_min:end])]
				lam_t = KP.x_l2reg2!(cs, muhat_t, vs, GammaLOO, x_t)
				thetaval = dot(thetas, x_t)/n
				writecsv(f, [iRun n "LOO_5" thetaval t GammaLOO])

				# ind_min = findfirst(Gamma_grid .>= 10.0)
				# GammaLOO = Gamma_grid[ind_min:end][indmax(objs[ind_min:end])]
				# xs = KP.x_l2reg(cs, muhat_t, vs, GammaLOO)[1]
				# thetaval = dot(thetas, xs)/n
				# writecsv(f, [iRun n "LOO_10" thetaval t GammaLOO])
			end

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
						theta_l, v_l, c_l, theta_m, v_m, c_m, theta_h, v_h, c_h, 
						includeReg, Gamma_min, Gamma_max)
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
	test_harness(f, numRuns, o, n_grid; includeReg=includeReg, Gamma_min=Gamma_min, Gamma_max=Gamma_max)
	close(f)
	return file_name
end


###  Reads in a theta/cs/vs specification
function test_ReadData(file_out, numRuns, n_grid, seed, param_path; 
						includeReg = true, Gamma_min = 1., Gamma_max = 20)
	srand(seed)
	const n_max = maximum(n_grid)
	dat, header = readcsv(param_path, header=true)

	#confirm that n_max works
	@assert n_max <= size(dat, 1) "Param file too short for n_max"
	cs = dat[1:n_max, 3]
	thetas = dat[1:n_max, 1]
	vs = dat[1:n_max, 2]

	#rescale cs so budget is sensible. 
	#VG To generate the LOO Plots, you need to undo this scaling.
	cs /= quantile(cs, .2)

	o = DefaultExp(cs, thetas, vs)
	file_name = "$(file_out)_$(seed).csv"
	f = open(file_name, "w")
	test_harness(f, numRuns, o, n_grid; includeReg=includeReg, Gamma_min=Gamma_min, Gamma_max=Gamma_max)
	close(f)
	return file_name	
end

#uses the three-part set-up to create density plots
function create_density_plot(file_name; n=2^17, seed=8675309)
	srand(seed)
	f = open(file_name, "w")

	cs = ones(n)
	thetas = ones(n)
	vs = ones(n)

	thetas[1:3:n] = theta_l
	vs[1:3:n] = v_l
	cs[1:3:n] = c_l

	thetas[2:3:n] = theta_m
	vs[2:3:n] = v_m
	cs[2:3:n] = c_m

	thetas[3:3:n] = theta_h
	vs[3:3:n] = v_h
	cs[3:3:n] = c_h
end

function tauDependence3PartPlot(file_name, n, seed)
	srand(seed)
	f = open("$(file_name)_$(n)_$(seed).csv", "w")
	#present the true value and relative to full-info
	writecsv(f, ["n" "tauVal" "thetaVal" "RelFullInfo"])

	#sim some data
	theta_l, v_l, c_l = 0.0, 0.1, 20.001
	theta_m, v_m, c_m = 1.0, 1.0, 20.001
	theta_h, v_h, c_h = 0.3, 4.0, 20.001

	cs = ones(n)
	thetas = ones(n)
	vs = ones(n)

	thetas[1:3:n] = theta_l
	vs[1:3:n] = v_l
	cs[1:3:n] = c_l

	thetas[2:3:n] = theta_m
	vs[2:3:n] = v_m
	cs[2:3:n] = c_m

	thetas[3:3:n] = theta_h
	vs[3:3:n] = v_h
	cs[3:3:n] = c_h

	muhat = randn(n) ./ sqrt.(vs) + thetas

	#compute full-info for scaling
	xstar, lamstar = KP.x_dual(cs, thetas)
	zstar = dot(xstar, thetas)/n

	for tau in linspace(0, 5, 500)
		rs = shrink(muhat, vs, tau)
		xs, lam = KP.x_dual(cs, rs)
		thetaVal = dot(xs, thetas)/n
		writecsv(f, [n tau thetaVal thetaVal/zstar ])
	end
	close(f)
end





########################################################
#########
n_grid = [2^i for i = 5:8]
#small run for pre-compilation
#test_Gaussian("./temp/temp_Gaussian", 5, [100, 150], 87, 3, 1, 3)
#test_OddEven("./temp/temp_OddEvenReg", 5, [100, 150], 8675309000, 2.1, includeReg=true)
test_ReadData("./temp/temp_PortExp", 5, [100, 150], 8675309, "./Results/param_portExp_mtn1.csv")

@time test_ReadData("./temp/temp_PortExp", 5, [32], 8675309, "./Results/param_portExp_mtn1.csv")


