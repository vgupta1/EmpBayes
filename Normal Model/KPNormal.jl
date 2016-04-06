module KP
using Distributions
using Roots

shrink(Xs, taus, tau0) = (taus ./ (taus + tau0)) .* Xs

#Scales costs by n internally
function q(cs, rs)
    const n = length(rs)
    ratios = rs ./ cs
    cs = cs/n #makes a copy
    indx = sortperm(ratios, rev=true)
    out = zeros(Float64, length(indx))
    weight = 0.
    for ix = 1:length(indx)
    	#check if the item has positive value, else skip
    	if rs[indx[ix]] <=0. 
    		continue
    	end

        #check if next item fits and it has positive value
        if weight + cs[indx[ix]] <= 1
            out[indx[ix]] = 1.
            weight += cs[indx[ix]]
        else
            #take fractional part
            out[[indx[ix]]] = (1-weight)/cs[indx[ix]]
            break
        end
    end
    out
end

#used to in computing the entire tau-grid
function sbars(xs, cs, taus, f_indx)
    out = zeros(Float64, length(xs))
    for i = 1:length(out)
        out[i] = taus[i]*taus[f_indx] * (cs[i]*xs[f_indx] - cs[f_indx]*xs[i])
        out[i] /= cs[f_indx]*taus[i]*xs[i] - cs[i]*taus[f_indx]*xs[f_indx]
    end   
    out
end

#min_indx above the threshold
#returns -1 if there are no values above thresh
function indmin_thresh(vals, thresh)
	const n = length(vals)
	min_indx = -1
	min_val = Inf
	for ix = 1:n
		if thresh < vals[ix] < min_val
			min_indx = ix
			min_val = vals[min_indx]
		end
	end
	min_indx
end

## searches over bfs to figure out the best value of tau0
function best_q_tau(cs, xs, taus, ys)
	const TOL = 1e-8
	const n = length(xs)
	ratios = xs ./ cs
	indx = sortperm(ratios, rev=true)
	cs = cs/n  #makes a copy

	vals = Float64[]
	objs = Float64[]

	#first determine tau0 = 0 solution
	qs = zeros(Float64, n)
	total_weight = 0.
	frac_indx = -1
    for ix = 1:n
    	#check if the item has positive value, else skip
    	if xs[indx[ix]] <=0. 
    		continue
    	end

        #check if next item fits entirely
        if total_weight + cs[indx[ix]] <= 1
            qs[indx[ix]] = 1.
            total_weight += cs[indx[ix]]
        else
            #take fractional part
            qs[[indx[ix]]] = (1.-total_weight)/cs[indx[ix]]
            frac_indx = indx[ix]
            total_weight = 1.
            break
        end
    end

    ##check for a trivial solution
    if frac_indx == -1 || total_weight < 1
    	println("Trivial Problem Set-up:  All items fit.")
    	println("total_weight:\t $total_weight")
    	println("qs:\n", qs)
    	return qs, [0.], [dot(qs, ys)/n]
    end

    #Now search across all possible tau0, constructing grid as you go
	tau0 = 0.  #always represents current value
    svals = sbars(xs, cs, taus, frac_indx)
    max_obj = -1
    max_tau0 = 0.

   	iter = 0
    while tau0 > -1.
    	iter += 1
    	@assert iter < 1e7 "Max iterations exceeded"
    	# println("Tau0: $tau0")
    	# println("qs:\n", qs)
		# println("Svals:\n", svals)
		# println("Frac_Indx:\t $frac_indx")


	    #assess solution and record value
	    obj = dot(ys, qs)/n
		push!(objs, obj)
    	push!(vals, tau0)
    	if obj > max_obj
    		max_obj = obj
    		max_tau0 = tau0
    	end    

    	##Find the entering indx if it exists
    	indx_p = indmin_thresh(svals, tau0)
    	# println("Entering Index:\t", indx_p)
    	if indx_p == -1
    		tau0 = -1
    		continue
    	end
    
    	#pivot depends on whether element is already in basis or not.
    	if qs[indx_p] >= 1-TOL
    		#Case 1: Exiting Variable becomes zero
    		if cs[indx_p]/cs[frac_indx] <= 1 - qs[frac_indx]
    			qs[indx_p] = 0
    			qs[frac_indx] += cs[indx_p]/cs[frac_indx]
    			@assert qs[frac_indx] < 1  "Degeneracy?"

    			#frac_indx and s_vals stay the same
    			tau0 = svals[indx_p]
    		else #Case 2: Exiting Variable becomes fractional
    			qs[indx_p] -= (1.- qs[frac_indx]) * cs[frac_indx]/cs[indx_p]
    			@assert qs[indx_p] > 0 "Dengeracy?"

    			qs[frac_indx] = 1
    			tau0 = svals[indx_p]
    			frac_indx = indx_p
    			svals = sbars(xs, cs, taus, frac_indx)			
    		end
    	else  #qs[indx_p] == 0
	    	frac_weight = qs[frac_indx] * cs[frac_indx]
    		#Case 1: Entering Variable becomes 1
	    	if frac_weight > cs[indx_p]  
	    		qs[indx_p] = 1.
	    		qs[frac_indx] -= cs[indx_p]/cs[frac_indx]
	    		@assert qs[frac_indx] > 0 "Dengeracy?"
	    		
	    		#frac_indx and svals are unchanged
	    		tau0 = svals[indx_p]
	    	else #Case 2: Entering Variable becomes fractional
				qs[indx_p] = frac_weight / cs[indx_p]
				qs[frac_indx] = 0.
				tau0 = svals[indx_p]	
				frac_indx = indx_p
				svals = sbars(xs, cs, taus, frac_indx)
			end
		end
    end

    ##assume that solving one last time is short relative to costs
    ##have to unscale cs for this to work
    max_qs = KP.q(cs * n, shrink(xs, taus, max_tau0))
    return max_qs, vals, objs 
end

#Compute the method of moment estimator
function q_MM(cs, zs, taus)
	tau_mm = 1/mean(zs.^2 - 1./taus)
	tau_mm = max(0, tau_mm)
	tau_mm, q(cs, shrink(zs, taus, tau_mm))
end

function q_MLE(cs, zs, taus, max_bnd = 1e2)
	#solve the MLE using a rootfinder solver for the variance and then invert
	vs = 1./ taus
	deriv_loglik(v) = mean(zs.^2 ./ (v + vs).^2 - 1./(v + vs))
	v0 = fzero(deriv_loglik, 0, max_bnd)
	@assert v0 < max_bnd "max_bnd $(max_bnd) was insufficient"
	1/v0, q(cs, shrink(zs, taus, 1/v0))
end

ideal_val(thetas, cs) = dot(q(cs, thetas), thetas) /length(thetas)

## returns the best of n^2 draws
function rand_alg(cs, xs, ys, taus, numDraws=length(xs)^2)
	const n = length(xs)
	qalg = zeros(Float64, n)
	sigs = 1./sqrt(taus)
	max_val = 0.

	for i=1:numDraws
		rs = xs + randn(n) .* sigs 
		q_temp = q(cs, rs)
		if (ys .* q_temp)/n > qalg
			qalg = q_temp
			max_val = (ys .* q_temp)/n
		end
	end 
	qalg 
end

###### Testing Code ####
# VG think about moving this to a standalone file
# Collapse these tests into a single harness with a set-up protocol type

####
#For increasing n, many simulations compare
#		#Naive x method
#		#Naive z method
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
#Costs and taus are fixed across all runs for now
#Uses a gaussian setup with tau0
function test1(file_out, numRuns, n_grid, seed, tau0)
	srand(seed)

	#output files
	f = open("$(file_out)_$(tau0)_$(seed).csv", "w")
	writecsv(f, ["Run" "n" "Method" "YVal" "thetaVal" "time" "tau0"])

	#keep costs and taus fixed across runs
	n_max = maximum(n_grid)
	cs = rand(n_max) * 20   #about 10% of items fit
	taus = 8. * rand(n_max)
	sigs = 1./sqrt(taus)

	#pre-allocate space for efficiency
	thetas = zeros(Float64, n_max)
	xs = zeros(Float64, n_max)
	ys = zeros(Float64, n_max)
	zs = zeros(Float64, n_max)
	noise = zeros(Float64, n_max)
	for iRun = 1:numRuns
		#generate the entire path
		thetas[:] = randn!(thetas) / sqrt(tau0)
		zs[:] = randn!(zs) .* sigs + thetas
		noise[:] = randn!(noise) .* sigs
		xs[:] = zs + noise
		ys[:] = zs - noise 

		for n in n_grid
			#Compute performance of each method
			#Naive x
			tic()
			qs = q(cs[1:n], xs[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(thetas[1:n], qs)/n
			writecsv(f, [iRun n "NaiveX" yval thetaval t 0.])

			#Naive Z
			tic()
			qs = q(cs[1:n], zs[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(thetas[1:n], qs)/n
			writecsv(f, [iRun n "NaiveZ" yval thetaval t 0.])

			#fullInfo val
			tic()
			qs = q(cs[1:n], thetas[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(thetas[1:n], qs)/n
			writecsv(f, [iRun n "FullInfo" yval thetaval t 0.])

			#The "Bayes" value.. only possible bc we know the setup is bayesian
			tic()
			qs = q(cs[1:n], shrink(zs[1:n], taus[1:n], tau0))
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(thetas[1:n], qs)/n
			writecsv(f, [iRun n "Bayes" yval thetaval t tau0])

			#Tau MLE
			tic()
			tauMLE, qs = q_MLE(cs[1:n], zs[1:n], taus[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(thetas[1:n], qs)/n
			writecsv(f, [iRun n "MLE" yval thetaval t tauMLE])

			#Tau MM
			tic()
			tauMM, qs = q_MM(cs[1:n], zs[1:n], taus[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(thetas[1:n], qs)/n
			writecsv(f, [iRun n "MM" yval thetaval t tauMM])

			#Our sample split method x
			tic()
			qs, vals, objs = best_q_tau(cs[1:n], xs[1:n], taus[1:n], ys[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(thetas[1:n], qs)/n
			writecsv(f, [iRun n "EmpBayesX" yval thetaval t vals[indmax(objs)]])

			#rescaling the sample_split method
			tau_RS = vals[indmax(objs)]/2
			qs = q(cs[1:n], shrink(zs[1:n], taus[1:n], tau_RS))
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(thetas[1:n], qs)/n
			writecsv(f, [iRun n "Rescaled" yval thetaval t tau_RS])

			#rescaling by a wrong factor
			taus_RS = vals[indmax(objs)]/2
			qs = q(cs[1:n], shrink(zs[1:n], taus[1:n], tau_RS))
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(thetas[1:n], qs)/n
			writecsv(f, [iRun n "Rescaled_1" yval thetaval t tau_RS])

			#rescaling by a wrong factor
			taus_RS = vals[indmax(objs)]/4
			qs = q(cs[1:n], shrink(zs[1:n], taus[1:n], tau_RS))
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(thetas[1:n], qs)/n
			writecsv(f, [iRun n "Rescaled_4" yval thetaval t tau_RS])

			#The ZZ method
			tic()
			qs, vals, objs = best_q_tau(cs[1:n], zs[1:n], taus[1:n], zs[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(thetas[1:n], qs)/n
			writecsv(f, [iRun n "ZZ" yval thetaval t vals[indmax(objs)]])

			#Oracle X method
			tic()
			qs, vals, objs = best_q_tau(cs[1:n], xs[1:n], taus[1:n], thetas[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(thetas[1:n], qs)/n
			writecsv(f, [iRun n "OracleX" yval thetaval t vals[indmax(objs)]])

			#Oracle Z method
			tic()
			qs, vals, objs = best_q_tau(cs[1:n], zs[1:n], taus[1:n], thetas[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(thetas[1:n], qs)/n
			writecsv(f, [iRun n "OracleZ" yval thetaval t vals[indmax(objs)]])

		end
		flush(f)
	end
	close(f)
end

####
# costs are same for every item,
# odd numbered items have lower theta, but low tau
# even numbered items have theta = 1, taus = 1 
function test2(file_out, numRuns, n_grid, seed, low_theta, low_tau)
	srand(seed)

	#output files
	f = open("$(file_out)_$(low_theta)_$(low_tau)_$(seed).csv", "w")
	writecsv(f, ["Run" "n" "Method" "YVal" "thetaVal" "time" "tau0"])

	n_max = maximum(n_grid)
	cs = 10 * ones(Float64, n_max)  #10% of items fit
	taus = ones(Float64, n_max)
	taus[1:2:n_max] = low_tau
	sigs = 1./sqrt(taus)
	thetas = ones(Float64, n_max)
	thetas[1:2:n_max] = low_theta

	#pre-allocate space for efficiency
	zs = zeros(Float64, n_max)
	xs = zeros(Float64, n_max)
	ys = zeros(Float64, n_max)
	noise = zeros(Float64, n_max)
	for iRun = 1:numRuns
		#generate the entire path
		zs[:] = randn!(zs) .* sigs + thetas
		noise[:] = randn!(noise) .* sigs
		xs[:] = zs + noise
		ys[:] = zs - noise 

		for n in n_grid
			#Compute performance of each method
			#Naive x
			tic()
			qs = q(cs[1:n], xs[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(thetas[1:n], qs)/n
			writecsv(f, [iRun n "NaiveX" yval thetaval t 0.])

			#Naive Z
			tic()
			qs = q(cs[1:n], zs[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(thetas[1:n], qs)/n
			writecsv(f, [iRun n "NaiveZ" yval thetaval t 0.])

			#fullInfo val
			tic()
			qs = q(cs[1:n], thetas[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(thetas[1:n], qs)/n
			writecsv(f, [iRun n "FullInfo" yval thetaval t 0.])

			#Tau MLE
			tic()
			tauMLE, qs = q_MLE(cs[1:n], zs[1:n], taus[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(thetas[1:n], qs)/n
			writecsv(f, [iRun n "MLE" yval thetaval t tauMLE])

			#Tau MM
			tic()
			tauMM, qs = q_MM(cs[1:n], zs[1:n], taus[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(thetas[1:n], qs)/n
			writecsv(f, [iRun n "MM" yval thetaval t tauMM])

			#Our sample split method x
			tic()
			qs, vals, objs = best_q_tau(cs[1:n], xs[1:n], taus[1:n], ys[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(thetas[1:n], qs)/n
			writecsv(f, [iRun n "EmpBayesX" yval thetaval t vals[indmax(objs)]])


			#rescaling the sample_split method
			tau_RS = vals[indmax(objs)]/2
			qs = q(cs[1:n], shrink(zs[1:n], taus[1:n], tau_RS))
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(thetas[1:n], qs)/n
			writecsv(f, [iRun n "Rescaled" yval thetaval t tau_RS])

			#Oracle X method
			tic()
			qs, vals, objs = best_q_tau(cs[1:n], xs[1:n], taus[1:n], thetas[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(thetas[1:n], qs)/n
			writecsv(f, [iRun n "OracleX" yval thetaval t vals[indmax(objs)]])

			#Oracle Z method
			tic()
			qs, vals, objs = best_q_tau(cs[1:n], zs[1:n], taus[1:n], thetas[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(thetas[1:n], qs)/n
			writecsv(f, [iRun n "OracleZ" yval thetaval t vals[indmax(objs)]])

		end
		flush(f)
	end
	close(f)
end

####
# costs are same for every item,
# computed at different scalings
function test3(file_out, numRuns, n_grid, seed, tau0, tau_scale_grid)
	srand(seed)

	#output files
	f = open("$(file_out)_tau_scale_$(seed).csv", "w")
	writecsv(f, ["Run" "n" "Scale" "YVal" "thetaVal" "time" "tau0"])

	n_max = maximum(n_grid)
	cs = rand(n_max) * 20   #about 10% of items fit
	taus = 8. * rand(n_max)
	sigs = 1./sqrt(taus)

	#pre-allocate space for efficiency
	thetas = zeros(Float64, n_max)
	zs = zeros(Float64, n_max)
	xs = zeros(Float64, n_max)
	ys = zeros(Float64, n_max)
	noise = zeros(Float64, n_max)
	for iRun = 1:numRuns
		#generate the entire path
		thetas[:] = randn!(thetas) / sqrt(tau0)
		zs[:] = randn!(zs) .* sigs + thetas
		noise[:] = randn!(noise) .* sigs
		xs[:] = zs + noise
		ys[:] = zs - noise 

		for n in n_grid
			#Compute performance of each method
			#Our sample split method x
			tic()
			qs, vals, objs = best_q_tau(cs[1:n], xs[1:n], taus[1:n], ys[1:n])
			t = toc()
			tau_hat = vals[indmax(objs)]

			for scale in tau_scale_grid
				qs = q(cs[1:n], shrink(zs[1:n], taus[1:n], tau_hat*scale))
				yval = dot(ys[1:n], qs)/n
				thetaval = dot(thetas[1:n], qs)/n
				writecsv(f, [iRun n scale yval thetaval t vals[indmax(objs)]])
			end
		end
		flush(f)
	end
	close(f)
end



####
# costs are same for every item,
# computed at different scalings
function test4(file_out, numRuns, n_grid, seed, low_theta, low_tau, tau_scale_grid)
	srand(seed)

	#output files
	f = open("$(file_out)_tau_scale_$(low_theta)_$(low_tau)_$(seed).csv", "w")
	writecsv(f, ["Run" "n" "Scale" "YVal" "thetaVal" "time" "tau0"])

	n_max = maximum(n_grid)
	cs = 10 * ones(Float64, n_max)  #10% of items fit
	taus = ones(Float64, n_max)
	taus[1:2:n_max] = low_tau
	sigs = 1./sqrt(taus)
	thetas = ones(Float64, n_max)
	thetas[1:2:n_max] = low_theta

	#pre-allocate space for efficiency
	zs = zeros(Float64, n_max)
	xs = zeros(Float64, n_max)
	ys = zeros(Float64, n_max)
	noise = zeros(Float64, n_max)
	for iRun = 1:numRuns
		#generate the entire path
		zs[:] = randn!(zs) .* sigs + thetas
		noise[:] = randn!(noise) .* sigs
		xs[:] = zs + noise
		ys[:] = zs - noise 

		for n in n_grid
			#Compute performance of each method
			#Our sample split method x
			tic()
			qs, vals, objs = best_q_tau(cs[1:n], xs[1:n], taus[1:n], ys[1:n])
			t = toc()
			tau_hat = vals[indmax(objs)]

			for scale in tau_scale_grid
				qs = q(cs[1:n], shrink(zs[1:n], taus[1:n], tau_hat*scale))
				yval = dot(ys[1:n], qs)/n
				thetaval = dot(thetas[1:n], qs)/n
				writecsv(f, [iRun n scale yval thetaval t vals[indmax(objs)]])
			end
		end
		flush(f)
	end
	close(f)
end

####
# costs are random uniform
# True generating model is random Exponential
# Testing to see how scale/2 performs
function test5(file_out, numRuns, n_grid, seed, alpha, beta)
	srand(seed)

	#output files
	f = open("$(file_out)_Exp_$(alpha)_$(beta)_$(seed).csv", "w")
	writecsv(f, ["Run" "n" "Method" "YVal" "thetaVal" "time" "tau0"])

	n_max = maximum(n_grid)
	cs = 20 * rand(Float64, n_max)  #10% of items fit

	##draw thetas from a gamma(alpha, beta)
	thetas = rand(Gamma(alpha, beta), n_max)
	taus = thetas.^2
	sigs = 1./thetas

	#pre-allocate space for efficiency
	zs = zeros(Float64, n_max)
	xs = zeros(Float64, n_max)
	ys = zeros(Float64, n_max)
	noise = zeros(Float64, n_max)
	for iRun = 1:numRuns
		#Generate the entire path
		for j = 1:n_max
			zs[j] = rand(Exponential(thetas[j]))
		end
		noise[:] = randn!(noise) .* sigs
		xs[:] = zs + noise
		ys[:] = zs - noise 

		for n in n_grid
			#Compute performance of each method
			#Naive Z
			tic()
			qs = q(cs[1:n], zs[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(thetas[1:n], qs)/n
			writecsv(f, [iRun n "NaiveZ" yval thetaval t 0.])

			#fullInfo val
			tic()
			qs = q(cs[1:n], thetas[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(thetas[1:n], qs)/n
			writecsv(f, [iRun n "FullInfo" yval thetaval t 0.])

			#VG Correct this to the exponential version
			# #The "Bayes" value.. only possible bc we know the setup is bayesian
			# tic()
			# qs = q(cs[1:n], shrink(zs[1:n], taus[1:n], tau0))
			# t = toc()
			# yval = dot(ys[1:n], qs)/n
			# thetaval = dot(thetas[1:n], qs)/n
			# writecsv(f, [iRun n "Bayes" yval thetaval t tau0])

			#Tau MLE
			tic()
			tauMLE, qs = q_MLE(cs[1:n], zs[1:n], taus[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(thetas[1:n], qs)/n
			writecsv(f, [iRun n "MLE" yval thetaval t tauMLE])

			#Tau MM
			tic()
			tauMM, qs = q_MM(cs[1:n], zs[1:n], taus[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(thetas[1:n], qs)/n
			writecsv(f, [iRun n "MM" yval thetaval t tauMM])

			#Our sample split method x
			tic()
			qs, vals, objs = best_q_tau(cs[1:n], xs[1:n], taus[1:n], ys[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(thetas[1:n], qs)/n
			writecsv(f, [iRun n "EmpBayesX" yval thetaval t vals[indmax(objs)]])

			#rescaling the sample_split method
			tau_RS = vals[indmax(objs)]/2
			qs = q(cs[1:n], shrink(zs[1:n], taus[1:n], tau_RS))
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(thetas[1:n], qs)/n
			writecsv(f, [iRun n "Rescaled" yval thetaval t tau_RS])

			#Oracle X method
			tic()
			qs, vals, objs = best_q_tau(cs[1:n], xs[1:n], taus[1:n], thetas[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(thetas[1:n], qs)/n
			writecsv(f, [iRun n "OracleX" yval thetaval t vals[indmax(objs)]])

			#Oracle Z method
			tic()
			qs, vals, objs = best_q_tau(cs[1:n], zs[1:n], taus[1:n], thetas[1:n])
			t = toc()
			yval = dot(ys[1:n], qs)/n
			thetaval = dot(thetas[1:n], qs)/n
			writecsv(f, [iRun n "OracleZ" yval thetaval t vals[indmax(objs)]])
		end
		flush(f)
	end
	close(f)
end




# #### Deprecated ######
# # compute the distribution of rewards for various taus
# function perf_taus(thetas, taus, cs, numSims, tau0_grid)
# 	const n = length(thetas)
# 	out = zeros(Float64, numSims, length(tau0_grid))
# 	for iSim = 1:numSims
# 		Xs = randn(n) ./ sqrt(taus) + thetas
# 		for (ix, t) in enumerate(tau0_grid)
# 			rs = shrink(Xs, taus, t)
# 			out[iSim, ix] = dot(q(cs, rs), thetas) / n
# 		end
# 	end
# 	out
# end

# ##Used to assess ULLN

# #Computes the set of taus where non-ngative xs cross.  
# function crossings(Xs, taus)
# 	const n = length(Xs)
# 	out = Float64[]
# 	for i = 1:n
# 		for j = i+1:n
# 			s = (taus[i]*taus[j]*Xs[j] - taus[i]*taus[j]*Xs[i]) 
# 			s /= (taus[i]*Xs[i] - taus[j]*Xs[j])
# 			if s > 0
# 				push!(out, s)
# 			end
# 		end
# 	end
# 	sort!(out)
# end

# function qbar(thetas, taus, cs, numSims, tau0::Float64)
# 	const n = length(thetas)
# 	out = zeros(Float64, length(thetas))
# 	for i = 1:numSims
# 		Xs = randn(n) ./ sqrt(taus) + thetas
# 		rs = shrink(Xs, taus, tau0)
# 		out += q(cs, rs) / numSims
# 	end
# 	out
# end

# function ULLN(thetas, taus, cs, Xs, numSims)
# 	const n = length(thetas)
# 	s_grid = crossings(Xs, taus)

# 	#computes qbar obj at each of the above crossign points
# 	qbar_obj = zeros(Float64, length(s_grid))
# 	for (ix, s) in enumerate(s_grid)
# 		qb = qbar(thetas, taus, cs, numSims, s)
# 		qbar_obj[ix] = dot(qb, thetas)/n
# 	end

# 	#compute qiobj at all the crossings
# 	qi_obj = zeros(Float64, length(s_grid))
# 	for (ix, s) in enumerate(s_grid)
# 		qi = q(cs, shrink(Xs, taus, s))
# 		qi_obj[ix] = dot(qi, thetas)/n
# 	end

# 	#diffs and upperbound
# 	diffs = abs(qbar_obj - qi_obj)
# 	j = indmax(diffs)
# 	qi_obj[j], qbar_obj[j], s_grid[j], s_grid[j+1] -s_grid[j]
# end

function q_linreg(cs, rs, ts)
	a, b = linreg(rs, ts) 
	rhat = a + b*rs
	q(cs, rhat)
end

function perf_reg(thetas, cs, numSims)
	const n = length(thetas)
	out = zeros(Float64, n, 3)
	for iSim = 1:numSims
		Xs = randn(n) ./ sqrt(taus) + thetas
		Ys = randn(n) ./ sqrt(taus) + thetas
		q = q_linreg(cs, Xs, Ys)
		out[iSim, 1] = dot(q, thetas) / n
	end
	out
end



end #ends module
