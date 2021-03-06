##Tests used in paper for the Small-Data KP
using Distributions, Random, DelimitedFiles, LinearAlgebra
include("KPNormal.jl")
using .KP

#Default experiments fix a single theta, cs, vs across all simulation
mutable struct DefaultExp
	cs
	thetas
	vs
end

function sim!(o, muhat)
	muhat[:] = randn!(muhat) ./ sqrt.(o.vs) .+ o.thetas
end

#simulates random normal stuff, but via a dataset
#dat is assumed n, S
function sim!(o, muhat, dat)
	S = size(dat, 2)
	n = size(dat, 1)
	for i = 1:S
		dat[:, i] .= randn(n) ./ sqrt.(o.vs/S)  .+ o.thetas #notice scaling by S
	end
	muhat[:] = mean(dat, dims=2)
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
function test_harness(f, numRuns, o, n_grid; Gamma_min=1., Gamma_max=20., Gamma_step=.01)
	n_max = maximum(n_grid)
	S = 10
	muhat = zeros(Float64, n_max)
	xs  = zeros(Float64, n_max)
	dat = zeros(Float64, n_max, S)
	lam_t = 0.

	#write a header
	writedlm(f,  ["Run" "n" "Method" "thetaVal" "time" "tau0"], ',')

	for iRun = 1:numRuns
		#generate the entire path up to n_max
		#sim!(o, muhat)
		sim!(o, muhat, dat)

		for n in n_grid
			#Take Views on Everything to speed up garbage collection
			x_t = view(xs, 1:n)
			vs  = view(o.vs, 1:n)
			cs  = view(o.cs, 1:n)
			thetas = view(o.thetas, 1:n)
			muhat_t = view(muhat, 1:n)
			dat_t = view(dat, 1:n, 1:S)

			#Compute performance of each method

			#SAA
			t = @elapsed x_t[:] = x(cs, muhat_t)
			thetaval = dot(thetas, x_t)/n
			writedlm(f,  [iRun n "SAA" thetaval t 0.], ',')

			#fullInfo val
			t = @elapsed x_t[:] = x(cs, thetas)
			thetaval = dot(thetas, x_t)/n
			writedlm(f,  [iRun n "FullInfo" thetaval t 0.], ',')

			# #The "Bayes" value.. only possible bc we if we we are in bayesian
			# if :tau0 in fieldnames(o)
			# 	tic()
			# 	xs = x(cs, shrink(muhat_t, vs, o.tau0))
			# 	t = toc()
			# 	thetaval = dot(thetas, xs)/n
			# 	writedlm(f,  [iRun n "Bayes" thetaval t o.tau0], ',')
			# end

			#Tau MLE
			t = 
			  @elapsed tauMLE, x_t[:] = x_MLE(cs, muhat_t, vs)
			
			thetaval = dot(thetas, x_t)/n
			writedlm(f,  [iRun n "EB_MLE" thetaval t tauMLE], ',')

			#Tau MM
			t = 
			  @elapsed tauMM, x_t[:] = x_MM(cs, muhat_t, vs)
			thetaval = dot(thetas, x_t)/n
			writedlm(f,  [iRun n "EB_MM" thetaval t tauMM], ',')

			#Oracle MSE
			t = 
			  @elapsed x_t[:], tau_CV = x_OR_MSE(cs, muhat_t, thetas, vs)
			
			thetaval = dot(thetas, x_t)/n
			writedlm(f,  [iRun n "OR_MSE" thetaval t tau_CV], ',')

			#Sure MSE
			t = 
			  @elapsed x_t[:], tau_CV = x_sure_MSE(cs, muhat_t, vs)
			
			thetaval = dot(thetas, x_t)/n
			writedlm(f,  [iRun n "SURE_MSE" thetaval t tau_CV], ',')

			#Dirac Stein
			t = 
			  @elapsed x_t[:], vals, objs = x_stein_exact(cs, muhat_t, vs, thetas)

			thetaval = dot(thetas, x_t)/n
			writedlm(f,  [iRun n "DiracStein" thetaval t vals[argmax(objs)]], ',')

			#Box with the optimized rate, i.e. h_n = n^-1/6 and scaling, altKernel
			h = n^-.16666
			t = 
			  @elapsed x_t[:], vals, objs = x_stein_box(cs, muhat_t, vs, h, tau_step = .05)
			
			thetaval = dot(thetas, x_t)/n
			writedlm(f,  [iRun n "BoxStein" thetaval t vals[argmax(objs)]], ',')

			#Oracle Value
			t = 
			  @elapsed x_t[:], vals, objs = best_x_tau(cs, muhat_t, vs, thetas)
			
			thetaval = dot(thetas, x_t)/n
			writedlm(f,  [iRun n "OR" thetaval t vals[argmax(objs)]], ',')

			#Hold-Out validation
			t = 
			  @elapsed x_t[:], tau_grid, objs = KP.x_kFoldCV(cs, vs/S, dat_t, 2; tau_grid=collect(0.:.05:5))
			
			thetaval = dot(thetas, x_t)/n
			writedlm(f,  [iRun n "HO_EB" thetaval t tau_grid[argmax(objs)]], ',')

			#5-Fold validation
			t = 
			  @elapsed x_t[:], tau_grid, objs = KP.x_kFoldCV(cs, vs/S, dat_t, 5; tau_grid=collect(0.:.05:5))
			
			thetaval = dot(thetas, x_t)/n
			writedlm(f,  [iRun n "K5_EB" thetaval t tau_grid[argmax(objs)]], ',')

			#LOO validation
			t = 
			  @elapsed x_t[:], tau_grid, objs = KP.x_kFoldCV(cs, vs/S, dat_t, S; tau_grid=collect(0.:.05:5))
			
			thetaval = dot(thetas, x_t)/n
			writedlm(f,  [iRun n "LOO_EB" thetaval t tau_grid[argmax(objs)]], ',')


			###Regularization Methods
			#Oracle Regularization
			t = 
			  @elapsed x_t[:], Gamma_grid, objs = KP.x_l2reg_CV(cs, muhat_t, vs, thetas, 
														Gamma_min=Gamma_min, Gamma_max=Gamma_max, Gamma_step=Gamma_step)
			
			Gammahat = Gamma_grid[argmax(objs)]
			thetaval = dot(thetas, x_t)/n
			writedlm(f,  [iRun n "OracleReg" thetaval t Gammahat], ',')

			# #From the same values, extract oracles for Other pairs
			# #For Gamma_min = 5
			ind_min = findfirst(Gamma_grid .>= 5.0)
			Gammahat = Gamma_grid[ind_min:end][argmax(objs[ind_min:end])]
			lam_t = KP.x_l2reg2!(cs, muhat_t, vs, Gammahat, x_t)
			thetaval = dot(thetas, x_t)/n
			writedlm(f,  [iRun n "OracleReg_5" thetaval t Gammahat], ',')

			# #Our Stein Approach to Regularization
			t = 
			  @elapsed x_t[:], Gamma_grid, objs = KP.x_stein_reg(cs, muhat_t, vs, 
										Gamma_min=Gamma_min, Gamma_max=Gamma_max, Gamma_step=Gamma_step)
			
			Gammahat = Gamma_grid[argmax(objs)]
			thetaval = dot(thetas, x_t)/n
			writedlm(f,  [iRun n "SteinReg" thetaval t Gammahat], ',')

			# #Again, from same values, extract values for gamma_min
			ind_min = findfirst(Gamma_grid .>= 5.0)
			Gammahat = Gamma_grid[ind_min:end][argmax(objs[ind_min:end])]
			lam_t = KP.x_l2reg2!(cs, muhat_t, vs, Gammahat, x_t)
			thetaval = dot(thetas, x_t)/n
			writedlm(f,  [iRun n "SteinReg_5" thetaval t Gammahat], ',')

			# #NEW RO method with new threshold 		
			thresh = sqrt(2*log(1/.1))
			t = 
			  @elapsed x_t[:] = KP.x_robFW(cs, muhat_t, vs, thresh, TOL=1e-4)
			
			thetaval = dot(thetas, x_t)/n
			writedlm(f,  [iRun n "FWRO_Eps_.1" thetaval t thresh], ',')

			thresh = sqrt(2*log(1/.05))
			t = 
			  @elapsed x_t[:] = KP.x_robFW(cs, muhat_t, vs, thresh, TOL=1e-4)
			
			thetaval = dot(thetas, x_t)/n
			writedlm(f,  [iRun n "FWRO_Eps_.05" thetaval t thresh], ',')

			thresh = sqrt(2*log(1/.01))
			t = 
			  @elapsed x_t[:] = KP.x_robFW(cs, muhat_t, vs, thresh, TOL=1e-4)
			
			thetaval = dot(thetas, x_t)/n
			writedlm(f,  [iRun n "FWRO_Eps_.01" thetaval t thresh], ',')


			# ####
			# #This section was used for the intial paper submission Feb2017
			# # The thresholds have been updated.  If want to rerun, consider using the FW implementation too.  
			# #RO heuristic for Gamma
			# #eps = .05				
			# tic()
			# xs, lam = KP.x_rob(cs, muhat_t, vs, 1.6448536269514717)
			# t = toc()
			# thetaval = dot(thetas, xs)/n
			# writedlm(f,  [iRun n "RO_Eps_.05" thetaval t 1.6448536269514717], ',')

			# #eps = .01				
			# tic()
			# xs, lam = KP.x_rob(cs, muhat_t, vs, 2.326347874040845)
			# t = toc()
			# thetaval = dot(thetas, xs)/n
			# writedlm(f,  [iRun n "RO_Eps_.01" thetaval t 2.326347874040845], ',')

			#VG Replaced with proper KFold
			#Leave one out validation (LOO)
			# tic()
			# x_t[:], Gamma_grid, objs = KP.x_LOO_reg(cs, muhat_t + noise_t, muhat_t - noise_t, vs, 
			# 									Gamma_min=Gamma_min, Gamma_max=Gamma_max, Gamma_step=Gamma_step)
			# t = toc()
			# thetaval = dot(thetas, x_t)/n
			# GammaLOO = Gamma_grid[argmax(objs)]
			# writedlm(f,  [iRun n "LOO" thetaval t GammaLOO], ',')

			# #Use same values to extract for other gamma_min
			# ind_min = findfirst(Gamma_grid .>= 5.0)
			# GammaLOO = Gamma_grid[ind_min:end][argmax(objs[ind_min:end])]
			# lam_t = KP.x_l2reg2!(cs, muhat_t, vs, GammaLOO, x_t)
			# thetaval = dot(thetas, x_t)/n
			# writedlm(f,  [iRun n "LOO_5" thetaval t GammaLOO], ',')
			#######
			#########	
			#######

			#Hold Out Validation choice
			t = 
			  @elapsed x_t[:], Gamma_grid, objs = KP.x_l2reg_kFoldCV(cs, vs/S, dat_t, 2; 
                    Gamma_step=Gamma_step, Gamma_min=Gamma_min, Gamma_max=Gamma_max)
			
			Gammahat = Gamma_grid[argmax(objs)]
			thetaval = dot(thetas, x_t)/n
			writedlm(f,  [iRun n "HO_Reg" thetaval t Gammahat], ',')

			# #5-Fold Validation 
			t = 
			  @elapsed x_t[:], Gamma_grid, objs = KP.x_l2reg_kFoldCV(cs, vs/S, dat_t, 5; 
                    Gamma_step=Gamma_step, Gamma_min=Gamma_min, Gamma_max=Gamma_max)
			
			Gammahat = Gamma_grid[argmax(objs)]
			thetaval = dot(thetas, x_t)/n
			writedlm(f,  [iRun n "K5_Reg" thetaval t Gammahat], ',')

			#LOO Validation 
			t = 
			  @elapsed x_t[:], Gamma_grid, objs = KP.x_l2reg_kFoldCV(cs, vs/S, dat_t, S; 
	                Gamma_step=Gamma_step, Gamma_min=Gamma_min, Gamma_max=Gamma_max)
			Gammahat = Gamma_grid[argmax(objs)]
			thetaval = dot(thetas, x_t)/n
			writedlm(f,  [iRun n "LOO_Reg" thetaval t Gammahat], ',')
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
	n_max = maximum(n_grid)
	Random.seed!(seed)
	cs = rand(n_max) .* 2 ./ frac_fit   
	vs = 2. * avg_tau .* rand(n_max)
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
						even_theta=1., frac_fit=.1)
	#build the sim object
	Random.seed!(seed)
	n_max = maximum(n_grid)
	cs = 1/frac_fit .* ones(n_max)
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

### The three parts experimental set-up
function test_threePart(file_out, numRuns, n_grid, seed, 
						theta_l, v_l, c_l, theta_m, v_m, c_m, theta_h, v_h, c_h, 
						Gamma_min, Gamma_max)
	Random.seed!(seed)
	n_max = maximum(n_grid)
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
	test_harness(f, numRuns, o, n_grid; Gamma_min=Gamma_min, Gamma_max=Gamma_max)
	close(f)
	return file_name
end


###  Reads in a theta/cs/vs specification
function test_ReadData(file_out, numRuns, n_grid, seed, param_path; 
						Gamma_min = 1., Gamma_max = 20, Gamma_step=.01)
	Random.seed!(seed)
	n_max = maximum(n_grid)
	dat, header = readdlm(param_path, ',', header=true)

	#confirm that n_max works
	@assert n_max <= size(dat, 1) "Param file too short for n_max"
	cs = dat[1:n_max, 3]
	thetas = dat[1:n_max, 1]
	vs = dat[1:n_max, 2]

	o = DefaultExp(cs, thetas, vs)
	file_name = "$(file_out)_$(seed).csv"
	f = open(file_name, "w")
	test_harness(f, numRuns, o, n_grid; Gamma_min=Gamma_min, Gamma_max=Gamma_max, Gamma_step=Gamma_step)
	close(f)
	return file_name	
end

#uses the three-part set-up to create density plots
function create_density_plot(file_name; n=2^17, seed=8675309)
	Random.seed!(seed)
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
	Random.seed!(seed)
	f = open("$(file_name)_$(n)_$(seed).csv", "w")
	#present the true value and relative to full-info
	writedlm(f,  ["n" "tauVal" "thetaVal" "RelFullInfo"], ',')

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
		writedlm(f,  [n tau thetaVal thetaVal/zstar ], ',')
	end
	close(f)
end





########################################################
#########
#Small run to precompiles stuff
#test_ReadData("./temp/temp_PortExp", 5, [100, 150], 8675309, "./Results/param_portExp_mtn2.csv")


