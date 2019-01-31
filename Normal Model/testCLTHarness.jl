using Distributions, Random, DelimitedFiles

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
	N
	dist
end

#VG Probably an opportunity to speed this up via a @. operator
#Simulates by drawing N uniform obs, and scaling by sqrt N
function sim!(o::CLTExp, muhat)
	fill!(muhat, 0.)
	#Shift and scale so that after loop, 
	#muhat is standard noise
	#mean of uniform is .5, and precision = 12
	shift = mean(o.dist)
	scale = 1/ std(o.dist) / sqrt(o.N)
	n = length(o.cs)
	for k = 1:o.N
		muhat[:] += (rand(o.dist, n) .- shift) .* scale
	end
	muhat[:] = muhat[:] ./ sqrt.(o.vs) + o.thetas
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

##VG a smarter solution would keep the same observations as you incease N
function test_CLTharness(f, numRuns, o, N_grid; Gamma_min = 1.0, Gamma_max=10, Gamma_step = .1)

	n = length(o.cs)
	muhat = zeros(Float64, n)
	noise = zeros(Float64, n)
	xs = zeros(Float64, n)
	lam_t = 0.

	#write a header
	writedlm(f,  ["Run" "N" "Method" "thetaVal" "time" "tau0"], ',')

	for iRun = 1:numRuns
		for N in N_grid
			#reset the object and simulate
			o.N = N
			sim!(o, muhat)

			#VG consider taking a view on muhat, vs, thtas, cs for performance

			#VG Consider updating this to do cross-val properly
			noise[:] = randn!(noise) ./ sqrt.(o.vs)	

			#SAA
			tic()
			xs[:] = KP.x(o.cs, muhat)
			t = toc()
			thetaval = dot(o.thetas, xs)/n
			writedlm(f,  [iRun N "SAA" thetaval t 0.], ',')

			#fullInfo val
			tic()
			xs[:] = x(o.cs, o.thetas)
			t = toc()
			thetaval = dot(o.thetas, xs)/n
			writedlm(f,  [iRun N "FullInfo" thetaval t 0.], ',')

			#Tau MLE
			tic()
			tauMLE, xs[:] = x_MLE(o.cs, muhat, o.vs)
			t = toc()
			thetaval = dot(o.thetas, xs)/n
			writedlm(f,  [iRun N "EB_MLE" thetaval t tauMLE], ',')

			#Tau MM
			tic()
			tauMM, xs[:] = x_MM(o.cs, muhat, o.vs)
			t = toc()
			thetaval = dot(o.thetas, xs)/n
			writedlm(f,  [iRun N "EB_MM" thetaval t tauMM], ',')

			#Oracle MSE
			tic()
			xs[:], tau_CV = x_OR_MSE(o.cs, muhat, o.thetas, o.vs)
			t = toc()
			thetaval = dot(o.thetas, xs)/n
			writedlm(f,  [iRun N "OR_MSE" thetaval t tau_CV], ',')

			#Sure MSE
			tic()
			xs[:], tau_CV = x_sure_MSE(o.cs, muhat, o.vs)
			t = toc()
			thetaval = dot(o.thetas, xs)/n
			writedlm(f,  [iRun N "SURE_MSE" thetaval t tau_CV], ',')

			#Box with the optimized rate, i.e. h_n = n^-1/6
			h = n^-.16666
			tic()
			xs[:], vals, objs = x_stein_box(o.cs, muhat, o.vs, h, tau_step = .05)
			t = toc()
			thetaval = dot(o.thetas, xs)/n
			writedlm(f,  [iRun N "BoxStein" thetaval t vals[argmax(objs)]], ',')

			#Oracle Value
			tic()
			xs[:], vals, objs = best_x_tau(o.cs, muhat, o.vs, o.thetas)
			t = toc()
			thetaval = dot(o.thetas, xs)/n
			writedlm(f,  [iRun N "OR" thetaval t vals[argmax(objs)]], ',')

			#Oracle Regularization
			tic()
			# xs, Gamma_grid, objs = KP.x_l2reg_CV_warm(o.cs, muhat, o.vs, o.thetas, 
			# 											Gamma_min=1., Gamma_max = 20.)
			xs[:], Gamma_grid, objs = KP.x_l2reg_CV(o.cs, muhat, o.vs, o.thetas, 
															Gamma_min=Gamma_min, Gamma_max=Gamma_max, Gamma_step=Gamma_step)


			t = toc()
			Gammahat = Gamma_grid[argmax(objs)]
			thetaval = dot(o.thetas, xs)/n
			writedlm(f,  [iRun N "OracleReg" thetaval t Gammahat], ',')

			## Reuse same values to look at different Gamma_min
			#For Gamma_min =5
			ind_min = findfirst(Gamma_grid .>= 5.0)
			Gammahat = Gamma_grid[ind_min:end][argmax(objs[ind_min:end])]
			#xs = KP.x_l2reg(o.cs, muhat, o.vs, Gammahat)[1]
			KP.x_l2reg2!(o.cs, muhat, o.vs, Gammahat, xs)
			thetaval = dot(o.thetas, xs)/n
			writedlm(f,  [iRun N "OracleReg_5" thetaval t Gammahat], ',')

			#Our Stein Approach to Regularization
			tic()
#			xs, Gamma_grid, objs = KP.x_stein_reg(o.cs, muhat, o.vs, Gamma_min=1, Gamma_max = 20.)
			xs[:], Gamma_grid, objs = KP.x_stein_reg(o.cs, muhat, o.vs, 
											Gamma_min=Gamma_min, Gamma_max=Gamma_max, Gamma_step=Gamma_step)

			t = toc()
			Gammahat = Gamma_grid[argmax(objs)]
			thetaval = dot(o.thetas, xs)/n
			writedlm(f,  [iRun N "SteinReg" thetaval t Gammahat], ',')

			#Reuse values for shortened gamma
			ind_min = findfirst(Gamma_grid .>= 5.0)
			Gammahat = Gamma_grid[ind_min:end][argmax(objs[ind_min:end])]
			#xs = KP.x_l2reg(o.cs, muhat, o.vs, Gammahat)[1]
			KP.x_l2reg2!(o.cs, muhat, o.vs, Gammahat, xs)
			thetaval = dot(o.thetas, xs)/n
			writedlm(f,  [iRun N "SteinReg_5" thetaval t Gammahat], ',')


			### Initial Paper submission:  The old code for robust
			### Now uses a better threshold and algorithm
			#RO heuristic for Gamma
			#eps = .05				
			# tic()
			# xs, lam = KP.x_rob(o.cs, muhat, o.vs, 1.6448536269514717)
			# t = toc()
			# thetaval = dot(o.thetas, xs)/n
			# writedlm(f,  [iRun N "RO_Eps_.05" thetaval t 1.6448536269514717], ',')

			# #eps = .01				
			# tic()
			# xs, lam = KP.x_l2reg(o.cs, muhat, o.vs, 2.326347874040845)
			# t = toc()
			# thetaval = dot(o.thetas, xs)/n
			# writedlm(f,  [iRun N "RO_Eps_.01" thetaval t 2.326347874040845], ',')

			tic()
			thresh = sqrt(2*log(1/.1))
			xs[:] = KP.x_robFW(o.cs, muhat, o.vs, thresh, TOL=1e-4)
			t = toc()
			thetaval = dot(o.thetas, xs)/n
			writedlm(f,  [iRun n "FWRO_Eps_.1" thetaval t thresh], ',')

			tic()
			thresh = sqrt(2*log(1/.05))
			xs[:] = KP.x_robFW(o.cs, muhat, o.vs, thresh, TOL=1e-4)
			t = toc()
			thetaval = dot(o.thetas, xs)/n
			writedlm(f,  [iRun n "FWRO_Eps_.05" thetaval t thresh], ',')

			tic()
			thresh = sqrt(2*log(1/.01))
			xs[:] = KP.x_robFW(o.cs, muhat, o.vs, thresh, TOL=1e-4)
			t = toc()
			thetaval = dot(o.thetas, xs)/n
			writedlm(f,  [iRun n "FWRO_Eps_.01" thetaval t thresh], ',')

			#Leave one out validation (LOO)
			tic()
			xs[:], Gamma_grid, objs = KP.x_LOO_reg(o.cs, muhat + noise, muhat - noise, o.vs, 
												Gamma_min=Gamma_min, Gamma_max=Gamma_max, Gamma_step=Gamma_step)

			t = toc()
			thetaval = dot(o.thetas, xs)/n
			writedlm(f,  [iRun N "LOO" thetaval t Gamma_grid[argmax(objs)]], ',')

			#Use same values to extract for other gamma_min
			ind_min = findfirst(Gamma_grid .>= 5.0)
			GammaLOO = Gamma_grid[ind_min:end][argmax(objs[ind_min:end])]
			#xs = KP.x_l2reg(o.cs, muhat, o.vs, GammaLOO)[1]
			KP.x_l2reg2!(o.cs, muhat, o.vs, GammaLOO, xs)
			thetaval = dot(o.thetas, xs)/n
			writedlm(f,  [iRun N "LOO_5" thetaval t GammaLOO], ',')

	
		end
		flush(f)
	end
end

function threePartCLTExp(n, N, dist)
	cs = ones(n)
	thetas = ones(n)
	vs = ones(n)

	theta_l, v_l, c_l = 0.0, 0.1, 20.001
	theta_m, v_m, c_m = 1.0, 1.0, 20.001
	theta_h, v_h, c_h = 0.3, 4.0, 20.001

	thetas[1:3:n] = theta_l
	vs[1:3:n] = v_l
	cs[1:3:n] = c_l

	thetas[2:3:n] = theta_m
	vs[2:3:n] = v_m
	cs[2:3:n] = c_m

	thetas[3:3:n] = theta_h
	vs[3:3:n] = v_h
	cs[3:3:n] = c_h

	return CLTExp(cs, thetas, vs, N, dist)
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

function test_3PartCLT(file_out, numRuns, n, N_grid, seed, dist_type)
	Random.seed!(seed)
	dist = get_dist(dist_type)
	o = threePartCLTExp(n, N_grid[1], dist)
	file_name = "$(file_out)_3partCLT_$(n)_$(seed).csv"
	f = open(file_name, "w")
	test_CLTharness(f, numRuns, o, N_grid)
	close(f)
	return file_name
end

function test_POAPCLT(file_out, param_path, numRuns, n, N_grid, seed, dist_type)
	Random.seed!(seed)
	dat, header = readdlm(param_path, ',', header=true)

	@assert n <= size(dat, 1) "Param file too short for n"
	cs = vec(dat[1:n, 3])
	thetas = vec(dat[1:n, 1])
	vs = vec(dat[1:n, 2])

	dist = get_dist(dist_type)
	o = CLTExp(cs, thetas, vs, N_grid[1], dist)

	file_name = "$(file_out)_$(dist_type)_$(seed).csv"
	f = open(file_name, "w")
	test_CLTharness(f, numRuns, o, N_grid)
	close(f)
	return file_name
end



########################################################
#########
#Small run for pre0compilation
#test_POAPCLT("Results/temp_POAPCLT", "Results/param_portExp_mtn2.csv", 2, 100, [2 3], 8675309, "exponential")


