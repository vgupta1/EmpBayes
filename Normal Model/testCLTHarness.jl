include("KPNormal.jl")
using KP, Distributions

########################################################
###
# Robustness to Non_Normality
##########
type CLTExp
	cs
	thetas
	vs
	N
	dist
end

#Simulates by drawing N uniform obs, and scaling by sqrt N
function sim!(o::CLTExp, muhat)
	fill!(muhat, 0.)
	#Shift and scale so that after loop, 
	#muhat is standard noise
	#mean of uniform is .5, and precision = 12
	const shift = mean(o.dist)#.5
	const scale = 1/ std(o.dist) / sqrt(o.N)  #sqrt(12/o.N)
	const n = length(o.cs)
	for k = 1:o.N
		muhat[:] += (rand(o.dist, n) - shift) * scale
	end
	muhat[:] = muhat[:] ./ sqrt(o.vs) + o.thetas
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
function test_CLTharness(f, numRuns, o, N_grid; includeReg=false)
	const n = length(o.cs)
	muhat = zeros(Float64, n)
	noise = zeros(Float64, n)

	#write a header
	writecsv(f, ["Run" "N" "Method" "thetaVal" "time" "tau0"])

	for iRun = 1:numRuns
		for N in N_grid
			#reset the object and simulate
			o.N = N
			sim!(o, muhat)
			noise[:] = randn!(noise) ./ sqrt(o.vs)


			#SAA
			tic()
			xs = KP.x(o.cs, muhat)
			t = toc()
			thetaval = dot(o.thetas, xs)/n
			writecsv(f, [iRun N "SAA" thetaval t 0.])

			#fullInfo val
			tic()
			xs = x(o.cs, o.thetas)
			t = toc()
			thetaval = dot(o.thetas, xs)/n
			writecsv(f, [iRun N "FullInfo" thetaval t 0.])

			#Tau MLE
			tic()
			tauMLE, xs = x_MLE(o.cs, muhat, o.vs)
			t = toc()
			thetaval = dot(o.thetas, xs)/n
			writecsv(f, [iRun N "EB_MLE" thetaval t tauMLE])

			#Tau MM
			tic()
			tauMM, xs = x_MM(o.cs, muhat, o.vs)
			t = toc()
			thetaval = dot(o.thetas, xs)/n
			writecsv(f, [iRun N "EB_MM" thetaval t tauMM])

			#Oracle MSE
			tic()
			xs, tau_CV = x_OR_MSE(o.cs, muhat, o.thetas, o.vs)
			t = toc()
			thetaval = dot(o.thetas, xs)/n
			writecsv(f, [iRun N "OR_MSE" thetaval t tau_CV])

			#Sure MSE
			tic()
			xs, tau_CV = x_sure_MSE(o.cs, muhat, o.vs)
			t = toc()
			thetaval = dot(o.thetas, xs)/n
			writecsv(f, [iRun N "SURE_MSE" thetaval t tau_CV])

			#Dirac Stein
			tic()
			xs, vals, objs = x_stein_exact(o.cs, muhat, o.vs, o.thetas)
			t = toc()
			thetaval = dot(o.thetas, xs)/n
			writecsv(f, [iRun N "DiracStein" thetaval t vals[indmax(objs)]])

			#Box with the optimized rate, i.e. h_n = n^-1/6
			h = n^-.16666
			tic()
			xs, vals, objs = x_stein_box(o.cs, muhat, o.vs, h, tau_step = .05)
			t = toc()
			thetaval = dot(o.thetas, xs)/n
			writecsv(f, [iRun N "BoxStein" thetaval t vals[indmax(objs)]])

			#Oracle Value
			tic()
			xs, vals, objs = best_x_tau(o.cs, muhat, o.vs, o.thetas)
			t = toc()
			thetaval = dot(o.thetas, xs)/n
			@assert abs(thetaval - maximum(objs)) <= 1e-5 "Weird Mismatch? \t $thetaval \t $(maximum(objs))"
			writecsv(f, [iRun N "OR" thetaval t vals[indmax(objs)]])

			if includeReg
				#Oracle Regularization
				tic()
				xs, Gamma_grid, objs = KP.x_l2reg_CV_warm(o.cs, muhat, o.vs, o.thetas, Gamma_max = 20.)
				t = toc()
				Gammahat = Gamma_grid[indmax(objs)]
				thetaval = dot(o.thetas, xs)/n
				writecsv(f, [iRun N "OracleReg" thetaval t Gammahat])

				#Our Stein Approach to Regularization
				tic()
				xs, Gamma_grid, objs = KP.x_stein_reg(o.cs, muhat, o.vs, Gamma_max = 20.)
				t = toc()
				Gammahat = Gamma_grid[indmax(objs)]
				thetaval = dot(o.thetas, xs)/n
				writecsv(f, [iRun N "SteinReg" thetaval t Gammahat])

				#Stein Appraoch with Bounds
				tic()
				xs, Gamma_grid, objs = KP.x_stein_reg(o.cs, muhat, o.vs, Gamma_min = 10., Gamma_max = 20.)
				t = toc()
				Gammahat = Gamma_grid[indmax(objs)]
				thetaval = dot(o.thetas, xs)/n
				writecsv(f, [iRun N "SteinRegBnded" thetaval t Gammahat])

				#RO heuristic for Gamma
				#eps = .1				
				tic()
				xs, lam = KP.x_rob(o.cs, muhat, o.vs, 1.2815515655446006)
				t = toc()
				thetaval = dot(o.thetas, xs)/n
				writecsv(f, [iRun N "RO_Eps_.1" thetaval t 1.2815515655446006])

				#eps = .05				
				tic()
				xs, lam = KP.x_rob(o.cs, muhat, o.vs, 1.6448536269514717)
				t = toc()
				thetaval = dot(o.thetas, xs)/n
				writecsv(f, [iRun N "RO_Eps_.05" thetaval t 1.6448536269514717])

				#eps = .01				
				tic()
				xs, lam = KP.x_l2reg(o.cs, muhat, o.vs, 2.326347874040845)
				t = toc()
				thetaval = dot(o.thetas, xs)/n
				writecsv(f, [iRun N "RO_Eps_.01" thetaval t 2.326347874040845])

				#Leave one out validation (LOO)
				tic()
				xs, Gamma_grid, objs = KP.x_LOO_reg(o.cs, muhat + noise, muhat - noise, o.vs)
				t = toc()
				thetaval = dot(o.thetas, xs)/n
				writecsv(f, [iRun N "LOO" thetaval t Gamma_grid[indmax(objs)]])



			end

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
	else
		throw("Distribution Type Not Recognized: $(dist_type)")
	end
	return dist
end

function test_3PartCLT(file_out, numRuns, n, N_grid, seed, dist_type)
	srand(seed)
	dist = get_dist(dist_type)
	o = threePartCLTExp(n, N_grid[1], dist)
	file_name = "$(file_out)_3partCLT_$(n)_$(seed).csv"
	f = open(file_name, "w")
	test_CLTharness(f, numRuns, o, N_grid, includeReg=true)
	close(f)
	return file_name
end

function test_POAPCLT(file_out, param_path, numRuns, n, N_grid, seed, dist_type)
	srand(seed)
	dat, header = readcsv(param_path, header=true)

	@assert n <= size(dat, 1) "Param file too short for n"
	cs = dat[1:n, 3]
	thetas = dat[1:n, 1]
	vs = dat[1:n, 2]
	cs /= quantile(cs, .2)  #rescale cs so budget is sensible. 

	dist = get_dist(dist_type)
	o = CLTExp(cs, thetas, vs, N_grid[1], dist)

	file_name = "$(file_out)_$(seed).csv"
	f = open(file_name, "w")
	test_CLTharness(f, numRuns, o, N_grid, includeReg=true)
	close(f)
	return file_name
end



########################################################
#########
N_grid = collect(1:10)
#Small run for pre0compilation
#test_3PartCLT("Results/temp_3PartCLT", 5, 100, [1 2], 8675309, "uniform")
test_POAPCLT("Results/temp_POAPCLT", "Results/param_portExp_mtn1.csv", 2, 100, [2 3], 8675309, "exponential")


