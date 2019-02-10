using Distributions, Random, DelimitedFiles, LinearAlgebra

include("KPNormal.jl")
using .KP
########################################################
###
# Robustness to Non_Normality
##########
mutable struct CLTExp
	cs
	thetas
	vs
	data  #contains the raw runs
	dist
end

#Simulates by drawing S obs according to dist, and scaling by sqrt S
# useConstantPrecision ensure that muhats have same precision for all S
# else it scales like sqrt(S)
# For efficiency, assumes data has been sized correctly!!
function sim!(o::CLTExp, muhat, useConstantPrecision)
	S = size(o.data, 2)
	n = size(o.data, 1)
	shifts = mean(o.dist)
	scales = 1/ std(o.dist) ./ sqrt.(o.vs)
	if useConstantPrecision
		scales[:] .*= sqrt(S)
	end 

	for i = 1:S
		#first generate variables with centered, scaled
		o.data[:, i] .= (rand(o.dist, n) .- mean(o.dist)) .* scales .+ o.thetas
	end
	muhat[:] = vec(mean(o.data, dims=2))
end

######################
####
#For increasing N, fixed n, many simulations compare
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
##### useConstantPrecision = true gives the CLT test.  false yields large sample test. 

function test_CLTharness(f, numRuns, o, S_grid, Gamma_min, Gamma_max, Gamma_step, useConstantPrecision)
	n = length(o.cs)
	muhat = zeros(Float64, n)
	xs = zeros(Float64, n)
	lam_t = 0.

	#write a header
	writedlm(f,  ["Run" "N" "Method" "thetaVal" "time" "tau0"], ',')

	for iRun = 1:numRuns
		for S in S_grid
			#reset the object and simulate
			o.data = zeros(n, S)
			sim!(o, muhat, useConstantPrecision)  #CLT experiment uses constant precision

			#SAA
			t = 
			  @elapsed xs[:] = KP.x(o.cs, muhat)
			thetaval = dot(o.thetas, xs)/n
			writedlm(f,  [iRun S "SAA" thetaval t 0.], ',')

			#fullInfo val
			t = 
			  @elapsed xs[:] = x(o.cs, o.thetas)
			thetaval = dot(o.thetas, xs)/n
			writedlm(f,  [iRun S "FullInfo" thetaval t 0.], ',')

			#Tau MLE
			t = 
			  @elapsed tauMLE, xs[:] = x_MLE(o.cs, muhat, o.vs)
			thetaval = dot(o.thetas, xs)/n
			writedlm(f,  [iRun S "EB_MLE" thetaval t tauMLE], ',')

			#Tau MM
			t = 
			  @elapsed tauMM, xs[:] = x_MM(o.cs, muhat, o.vs)
			thetaval = dot(o.thetas, xs)/n
			writedlm(f,  [iRun S "EB_MM" thetaval t tauMM], ',')

			#Oracle MSE
			t = 
			  @elapsed xs[:], tau_CV = x_OR_MSE(o.cs, muhat, o.thetas, o.vs)
			thetaval = dot(o.thetas, xs)/n
			writedlm(f,  [iRun S "OR_MSE" thetaval t tau_CV], ',')

			#Sure MSE
			t = 
			  @elapsed xs[:], tau_CV = x_sure_MSE(o.cs, muhat, o.vs)
			thetaval = dot(o.thetas, xs)/n
			writedlm(f,  [iRun S "SURE_MSE" thetaval t tau_CV], ',')

			#Box with the optimized rate, i.e. h_n = n^-1/6
			h = n^-.16666
			t = 
			  @elapsed xs[:], vals, objs = x_stein_box(o.cs, muhat, o.vs, h, tau_step = .05)
			thetaval = dot(o.thetas, xs)/n
			writedlm(f,  [iRun S "BoxStein" thetaval t vals[argmax(objs)]], ',')

			#Oracle Value
			t = 
			  @elapsed xs[:], vals, objs = best_x_tau(o.cs, muhat, o.vs, o.thetas)
			thetaval = dot(o.thetas, xs)/n
			writedlm(f,  [iRun S "OR" thetaval t vals[argmax(objs)]], ',')

			#Hold-Out validation
			if S >= 2
				t = 
				  @elapsed xs[:], tau_grid, objs = KP.x_kFoldCV(o.cs, o.vs/S, o.data, 2; tau_grid=collect(0.:.05:5))
				
				thetaval = dot(o.thetas, xs)/n
				writedlm(f,  [iRun S "HO_EB" thetaval t tau_grid[argmax(objs)]], ',')

			end

			#5-Fold validation
			if S >= 5
				t = 
				  @elapsed xs[:], tau_grid, objs = KP.x_kFoldCV(o.cs, o.vs/S, o.data, 5; tau_grid=collect(0.:.05:5))
				
				thetaval = dot(o.thetas, xs)/n
				writedlm(f,  [iRun S "K5_EB" thetaval t tau_grid[argmax(objs)]], ',')
			end

			#LOO validation
			if S >= 2
				t = 
				  @elapsed xs[:], tau_grid, objs = KP.x_kFoldCV(o.cs, o.vs/S, o.data, S; tau_grid=collect(0.:.05:5))
				
				thetaval = dot(o.thetas, xs)/n
				writedlm(f,  [iRun S "LOO_EB" thetaval t tau_grid[argmax(objs)]], ',')
			end

			####
			# Regularization 
			####
			#Oracle Regularization
			t = 
			  @elapsed xs[:], Gamma_grid, objs = KP.x_l2reg_CV(o.cs, muhat, o.vs, o.thetas, 
															Gamma_min=Gamma_min, Gamma_max=Gamma_max, Gamma_step=Gamma_step)


			Gammahat = Gamma_grid[argmax(objs)]
			thetaval = dot(o.thetas, xs)/n
			writedlm(f,  [iRun S "OracleReg" thetaval t Gammahat], ',')

			## Reuse same values to look at different Gamma_min
			#For Gamma_min =5
			ind_min = findfirst(Gamma_grid .>= 5.0)
			Gammahat = Gamma_grid[ind_min:end][argmax(objs[ind_min:end])]
			#xs = KP.x_l2reg(o.cs, muhat, o.vs, Gammahat)[1]
			t = 
			  @elapsed KP.x_l2reg2!(o.cs, muhat, o.vs, Gammahat, xs)
			thetaval = dot(o.thetas, xs)/n
			writedlm(f,  [iRun S "OracleReg_5" thetaval t Gammahat], ',')

			#Our Stein Approach to Regularization
			t = 
			  @elapsed xs[:], Gamma_grid, objs = KP.x_stein_reg(o.cs, muhat, o.vs, 
											Gamma_min=Gamma_min, Gamma_max=Gamma_max, Gamma_step=Gamma_step)

			Gammahat = Gamma_grid[argmax(objs)]
			thetaval = dot(o.thetas, xs)/n
			writedlm(f,  [iRun S "SteinReg" thetaval t Gammahat], ',')

			#Reuse values for shortened gamma
			ind_min = findfirst(Gamma_grid .>= 5.0)
			Gammahat = Gamma_grid[ind_min:end][argmax(objs[ind_min:end])]
			t = 
			  @elapsed KP.x_l2reg2!(o.cs, muhat, o.vs, Gammahat, xs)
			thetaval = dot(o.thetas, xs)/n
			writedlm(f,  [iRun S "SteinReg_5" thetaval t Gammahat], ',')


			### Initial Paper submission:  The old code for robust
			### Sow uses a better threshold and algorithm
			#RO heuristic for Gamma
			#eps = .05				
			# tic()
			# xs, lam = KP.x_rob(o.cs, muhat, o.vs, 1.6448536269514717)
			# t = toc()
			# thetaval = dot(o.thetas, xs)/n
			# writedlm(f,  [iRun S "RO_Eps_.05" thetaval t 1.6448536269514717], ',')

			# #eps = .01				
			# tic()
			# xs, lam = KP.x_l2reg(o.cs, muhat, o.vs, 2.326347874040845)
			# t = toc()
			# thetaval = dot(o.thetas, xs)/n
			# writedlm(f,  [iRun S "RO_Eps_.01" thetaval t 2.326347874040845], ',')

			thresh = sqrt(2*log(1/.1))
			t = 
			  @elapsed xs[:] = KP.x_robFW(o.cs, muhat, o.vs, thresh, TOL=1e-4)
			thetaval = dot(o.thetas, xs)/n
			writedlm(f,  [iRun S "FWRO_Eps_.1" thetaval t thresh], ',')

			thresh = sqrt(2*log(1/.05))
			t = 
			  @elapsed xs[:] = KP.x_robFW(o.cs, muhat, o.vs, thresh, TOL=1e-4)
			thetaval = dot(o.thetas, xs)/n
			writedlm(f,  [iRun S "FWRO_Eps_.05" thetaval t thresh], ',')

			thresh = sqrt(2*log(1/.01))
			t = 
			  @elapsed xs[:] = KP.x_robFW(o.cs, muhat, o.vs, thresh, TOL=1e-4)
			thetaval = dot(o.thetas, xs)/n
			writedlm(f,  [iRun S "FWRO_Eps_.01" thetaval t thresh], ',')


			#####
			## Add K Fold values IF S is adequately large
			#####
			#Hold Out Validation choice
			if S >= 2
				t = 
				  @elapsed xs[:], Gamma_grid, objs = KP.x_l2reg_kFoldCV(o.cs, o.vs/S, o.data, 2; 
	                    Gamma_step=Gamma_step, Gamma_min=Gamma_min, Gamma_max=Gamma_max)
				
				Gammahat = Gamma_grid[argmax(objs)]
				thetaval = dot(o.thetas, xs)/n
				writedlm(f,  [iRun S "HO_Reg" thetaval t Gammahat], ',')
			end

			# #5-Fold Validation 
			if S >= 5
				t = 
				  @elapsed xs[:], Gamma_grid, objs = KP.x_l2reg_kFoldCV(o.cs, o.vs/S, o.data, 5; 
	                    Gamma_step=Gamma_step, Gamma_min=Gamma_min, Gamma_max=Gamma_max)
				
				Gammahat = Gamma_grid[argmax(objs)]
				thetaval = dot(o.thetas, xs)/n
				writedlm(f,  [iRun S "K5_Reg" thetaval t Gammahat], ',')
			end

			#LOO Validation 
			if S >= 2
				t = 
				  @elapsed xs[:], Gamma_grid, objs = KP.x_l2reg_kFoldCV(o.cs, o.vs/S, o.data, S; 
		                Gamma_step=Gamma_step, Gamma_min=Gamma_min, Gamma_max=Gamma_max)
				Gammahat = Gamma_grid[argmax(objs)]
				thetaval = dot(o.thetas, xs)/n
				writedlm(f,  [iRun S "LOO_Reg" thetaval t Gammahat], ',')
			end

			# #Leave one out validation (LOO)
			# tic()
			# xs[:], Gamma_grid, objs = KP.x_LOO_reg(o.cs, muhat + noise, muhat - noise, o.vs, 
			# 									Gamma_min=Gamma_min, Gamma_max=Gamma_max, Gamma_step=Gamma_step)

			# t = toc()
			# thetaval = dot(o.thetas, xs)/n
			# writedlm(f,  [iRun S "LOO" thetaval t Gamma_grid[argmax(objs)]], ',')

			# #Use same values to extract for other gamma_min
			# ind_min = findfirst(Gamma_grid .>= 5.0)
			# GammaLOO = Gamma_grid[ind_min:end][argmax(objs[ind_min:end])]
			# #xs = KP.x_l2reg(o.cs, muhat, o.vs, GammaLOO)[1]
			# KP.x_l2reg2!(o.cs, muhat, o.vs, GammaLOO, xs)
			# thetaval = dot(o.thetas, xs)/n
			# writedlm(f,  [iRun S "LOO_5" thetaval t GammaLOO], ',')
		end
		flush(f)
	end
end

#helper
function get_dist(dist_type)
	if dist_type == "uniform"
		dist = Uniform()
	elseif dist_type == "bernoulli"
		dist = Bernoulli(.5)
	elseif dist_type == "exponential"
		dist = Exponential()
	elseif dist_type == "t"
		dist = TDist(3)
	elseif dist_type == "pareto"
		dist = Pareto(3)
	else
		throw("Distribution Type Not Recognized: $(dist_type)")
	end
	return dist
end

function test_POAPCLT(file_out, param_path, numRuns, n, S_grid, seed, dist_type, Gammamin, Gammamax, Gamma_step)
	Random.seed!(seed)
	dat, header = readdlm(param_path, ',', header=true)

	@assert n <= size(dat, 1) "Param file too short for n"
	cs = vec(dat[1:n, 3])
	thetas = vec(dat[1:n, 1])
	vs = vec(dat[1:n, 2])

	dist = get_dist(dist_type)
	o = CLTExp(cs, thetas, vs, S_grid[1], dist)

	file_name = "$(file_out)_$(dist_type)_$(seed).csv"
	f = open(file_name, "w")
	test_CLTharness(f, numRuns, o, S_grid, Gammamin, Gammamax, Gamma_step)
	close(f)
	return file_name
end



########################################################
#########
#Small run for pre0compilation
#test_POAPCLT("Results/temp_POAPCLT", "Results/param_portExp_mtn2.csv", 2, 100, [2 3], 8675309, "exponential")


