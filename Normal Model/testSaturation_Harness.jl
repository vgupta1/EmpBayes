##Saturation Tests
using Distributions, Random, DelimitedFiles, LinearAlgebra
include("KPNormal.jl")
using .KP

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
#For increasing n, many simulations compare methods
#on metrics
#		#Value wrt theta
#		#Time to compute
#		#optimal value of tau0
#Key is that the cs are scaled so that 
# the constraint is effectively 
#      c^T x \leq min(alpha*n, 1) --> 1/ n max(1/alpha ,n) *  c^T x \leq 1
#Since opt value is O(1), we do not scale opt value by n for stability.
function test_harness(f, numRuns, o, n_grid, Gamma_min, Gamma_max, Gamma_step, alpha)
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
		sim!(o, muhat, dat)

		for n in n_grid
			#Take Views on Everything to speed up garbage collection
			x_t = view(xs, 1:n)
			vs  = view(o.vs, 1:n)
			thetas = view(o.thetas, 1:n)
			muhat_t = view(muhat, 1:n)
			dat_t = view(dat, 1:n, 1:S)

			#rescale the c's so that they work properly n what follows.
			cs  = o.cs[1:n] * max(1/alpha, n)


			#Compute performance of each method
			#SAA
			t = @elapsed x_t[:] = x(cs, muhat_t)
			thetaval = dot(thetas, x_t)
			writedlm(f,  [iRun n "SAA" thetaval t 0.], ',')

			#fullInfo val
			t = @elapsed x_t[:] = x(cs, thetas)
			thetaval = dot(thetas, x_t)
			writedlm(f,  [iRun n "FullInfo" thetaval t 0.], ',')

			#Tau MLE
			t = 
			  @elapsed tauMLE, x_t[:] = x_MLE(cs, muhat_t, vs)
			
			thetaval = dot(thetas, x_t)
			writedlm(f,  [iRun n "EB_MLE" thetaval t tauMLE], ',')

			#Tau MM
			t = 
			  @elapsed tauMM, x_t[:] = x_MM(cs, muhat_t, vs)
			thetaval = dot(thetas, x_t)
			writedlm(f,  [iRun n "EB_MM" thetaval t tauMM], ',')

			#Oracle MSE
			t = 
			  @elapsed x_t[:], tau_CV = x_OR_MSE(cs, muhat_t, thetas, vs)
			
			thetaval = dot(thetas, x_t)
			writedlm(f,  [iRun n "OR_MSE" thetaval t tau_CV], ',')

			#Sure MSE
			t = 
			  @elapsed x_t[:], tau_CV = x_sure_MSE(cs, muhat_t, vs)
			
			thetaval = dot(thetas, x_t)
			writedlm(f,  [iRun n "SURE_MSE" thetaval t tau_CV], ',')

			#Dirac Stein
			t = 
			  @elapsed x_t[:], vals, objs = x_stein_exact(cs, muhat_t, vs, thetas)

			thetaval = dot(thetas, x_t)
			writedlm(f,  [iRun n "DiracStein" thetaval t vals[argmax(objs)]], ',')

			#Box with the optimized rate, i.e. h_n = n^-1/6 and scaling, altKernel
			h = n^-.16666
			t = 
			  @elapsed x_t[:], vals, objs = x_stein_box(cs, muhat_t, vs, h, tau_step = .05)
			
			thetaval = dot(thetas, x_t)
			writedlm(f,  [iRun n "BoxStein" thetaval t vals[argmax(objs)]], ',')

			#Oracle Value
			t = 
			  @elapsed x_t[:], vals, objs = best_x_tau(cs, muhat_t, vs, thetas)
			
			thetaval = dot(thetas, x_t)
			writedlm(f,  [iRun n "OR" thetaval t vals[argmax(objs)]], ',')

			#Hold-Out validation
			t = 
			  @elapsed x_t[:], tau_grid, objs = KP.x_kFoldCV(cs, vs/S, dat_t, 2; tau_grid=collect(0.:.05:5))
			
			thetaval = dot(thetas, x_t)
			writedlm(f,  [iRun n "HO_EB" thetaval t tau_grid[argmax(objs)]], ',')

			#5-Fold validation
			t = 
			  @elapsed x_t[:], tau_grid, objs = KP.x_kFoldCV(cs, vs/S, dat_t, 5; tau_grid=collect(0.:.05:5))
			
			thetaval = dot(thetas, x_t)
			writedlm(f,  [iRun n "K5_EB" thetaval t tau_grid[argmax(objs)]], ',')

			#LOO validation
			t = 
			  @elapsed x_t[:], tau_grid, objs = KP.x_kFoldCV(cs, vs/S, dat_t, S; tau_grid=collect(0.:.05:5))
			
			thetaval = dot(thetas, x_t)
			writedlm(f,  [iRun n "LOO_EB" thetaval t tau_grid[argmax(objs)]], ',')


			###Regularization Methods
			#Oracle Regularization
			t = 
			  @elapsed x_t[:], Gamma_grid, objs = KP.x_l2reg_CV(cs, muhat_t, vs, thetas, 
														Gamma_min=Gamma_min, Gamma_max=Gamma_max, Gamma_step=Gamma_step)
			
			Gammahat = Gamma_grid[argmax(objs)]
			thetaval = dot(thetas, x_t)
			writedlm(f,  [iRun n "OracleReg" thetaval t Gammahat], ',')

			# #Our Stein Approach to Regularization
			t = 
			  @elapsed x_t[:], Gamma_grid, objs = KP.x_stein_reg(cs, muhat_t, vs, 
										Gamma_min=Gamma_min, Gamma_max=Gamma_max, Gamma_step=Gamma_step)
			
			Gammahat = Gamma_grid[argmax(objs)]
			thetaval = dot(thetas, x_t)
			writedlm(f,  [iRun n "SteinReg" thetaval t Gammahat], ',')

			# #NEW RO method with new threshold 		
			thresh = sqrt(2*log(1/.1))
			t = 
			  @elapsed x_t[:] = KP.x_robFW(cs, muhat_t, vs, thresh, TOL=1e-4)
			
			thetaval = dot(thetas, x_t)
			writedlm(f,  [iRun n "FWRO_Eps_.1" thetaval t thresh], ',')

			thresh = sqrt(2*log(1/.05))
			t = 
			  @elapsed x_t[:] = KP.x_robFW(cs, muhat_t, vs, thresh, TOL=1e-4)
			
			thetaval = dot(thetas, x_t)
			writedlm(f,  [iRun n "FWRO_Eps_.05" thetaval t thresh], ',')

			thresh = sqrt(2*log(1/.01))
			t = 
			  @elapsed x_t[:] = KP.x_robFW(cs, muhat_t, vs, thresh, TOL=1e-4)
			
			thetaval = dot(thetas, x_t)
			writedlm(f,  [iRun n "FWRO_Eps_.01" thetaval t thresh], ',')

			#Hold Out Validation choice
			t = 
			  @elapsed x_t[:], Gamma_grid, objs = KP.x_l2reg_kFoldCV(cs, vs/S, dat_t, 2; 
                    Gamma_step=Gamma_step, Gamma_min=Gamma_min, Gamma_max=Gamma_max)
			
			Gammahat = Gamma_grid[argmax(objs)]
			thetaval = dot(thetas, x_t)
			writedlm(f,  [iRun n "HO_Reg" thetaval t Gammahat], ',')

			# #5-Fold Validation 
			t = 
			  @elapsed x_t[:], Gamma_grid, objs = KP.x_l2reg_kFoldCV(cs, vs/S, dat_t, 5; 
                    Gamma_step=Gamma_step, Gamma_min=Gamma_min, Gamma_max=Gamma_max)
			
			Gammahat = Gamma_grid[argmax(objs)]
			thetaval = dot(thetas, x_t)
			writedlm(f,  [iRun n "K5_Reg" thetaval t Gammahat], ',')

			#LOO Validation 
			t = 
			  @elapsed x_t[:], Gamma_grid, objs = KP.x_l2reg_kFoldCV(cs, vs/S, dat_t, S; 
	                Gamma_step=Gamma_step, Gamma_min=Gamma_min, Gamma_max=Gamma_max)
			Gammahat = Gamma_grid[argmax(objs)]
			thetaval = dot(thetas, x_t)
			writedlm(f,  [iRun n "LOO_Reg" thetaval t Gammahat], ',')
		end
		flush(f)
	end
end

#######################
###  Reads in a theta/cs/vs specification
function test_saturation(file_out, numRuns, n_grid, seed, param_path; 
						Gamma_min = 1., Gamma_max = 100., Gamma_step=.5, alpha=1.)
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
	test_harness(f, numRuns, o, n_grid, Gamma_min, Gamma_max, Gamma_step, alpha)
	close(f)
	return file_name	
end

