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
	vs   #always represents precisions of muhat
	dist

	data  #contains the raw runs, always size (n, S)
	vs_orig  #original precisions pre-scaling by S
	useConstantPrecision #ensures muhat has same precision for all S
end

function CLTExp(cs, thetas, vs_, dist, S, useConstantPrecision) 
	vs = copy(vs_)
	if !useConstantPrecision
		vs[:] .*= S
	end
	CLTExp(cs, thetas, vs, dist, zeros(length(cs), S), 
		copy(vs), useConstantPrecision)
end

#Simulates by drawing S obs according to dist, and scaling by sqrt S
#Assumes object pre-configured!
function sim!(o::CLTExp, muhat)
	S = size(o.data, 2)
	n = size(o.data, 1)
	shifts = mean(o.dist)
	scales = 1/ std(o.dist) ./ sqrt.(o.vs_orig)
	if o.useConstantPrecision
		scales[:] .*= sqrt(S)
	end 

	for i = 1:S
		#first generate variables with centered, scaled
		o.data[:, i] .= (rand(o.dist, n) .- mean(o.dist)) .* scales .+ o.thetas
	end
	muhat[:] = vec(mean(o.data, dims=2))
end

#sets object up to be used for size S
function resize!(o::CLTExp, S)
	o.data = zeros(length(o.cs), S)
	o.vs = copy(o.vs_orig)
	if !o.useConstantPrecision
		o.vs[:] .*= S
	end
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
function test_CLTharness(f, numRuns, o, S_grid, Gamma_min, Gamma_max, Gamma_step)
	n = length(o.cs)
	muhat = zeros(Float64, n)
	xs = zeros(Float64, n)
	lam_t = 0.

	#write a header
	writedlm(f,  ["Run" "N" "Method" "thetaVal" "time" "tau0"], ',')

	for iRun = 1:numRuns
		for S in S_grid
			#reset the object and simulate
			resize!(o, S)
			#o.data = zeros(n, S)
			sim!(o, muhat)

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

			#Oracle
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
			#Oracle
			t = 
			  @elapsed xs[:], Gamma_grid, objs = KP.x_l2reg_CV(o.cs, muhat, o.vs, o.thetas, 
															Gamma_min=Gamma_min, Gamma_max=Gamma_max, Gamma_step=Gamma_step)


			Gammahat = Gamma_grid[argmax(objs)]
			thetaval = dot(o.thetas, xs)/n
			writedlm(f,  [iRun S "OracleReg" thetaval t Gammahat], ',')

			#Our Stein Approach to Regularization
			t = 
			  @elapsed xs[:], Gamma_grid, objs = KP.x_stein_reg(o.cs, muhat, o.vs, 
											Gamma_min=Gamma_min, Gamma_max=Gamma_max, Gamma_step=Gamma_step)

			Gammahat = Gamma_grid[argmax(objs)]
			thetaval = dot(o.thetas, xs)/n
			writedlm(f,  [iRun S "SteinReg" thetaval t Gammahat], ',')

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
	o = CLTExp(cs, thetas, vs, dist, S_grid[1], true) #true forces constant precision

	file_name = "$(file_out)_$(dist_type)_$(seed).csv"
	f = open(file_name, "w")
	test_CLTharness(f, numRuns, o, S_grid, Gammamin, Gammamax, Gamma_step)  
	close(f)
	return file_name
end


function test_POAPLargeSample(file_out, param_path, numRuns, n, S_grid, seed, dist_type, Gammamin, Gammamax, Gamma_step)
	Random.seed!(seed)
	dat, header = readdlm(param_path, ',', header=true)

	@assert n <= size(dat, 1) "Param file too short for n"
	cs = vec(dat[1:n, 3])
	thetas = vec(dat[1:n, 1])
	vs = vec(dat[1:n, 2])

	dist = get_dist(dist_type)
	o = CLTExp(cs, thetas, vs, dist, S_grid[1], false) #false forces constant precision

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


