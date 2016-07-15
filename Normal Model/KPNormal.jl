module KP
using Distributions
using Roots

export q, best_q_tau, q_MM, q_MLE, ideal_val, q_dual, lam, shrink, q_ridge

shrink(xs, vs, tau0) = (vs ./ (vs + tau0)) .* xs

#cs_unscaled should have entries on order unity.
#The dual return fails with value -1 if the particular run is degenerate.  
function q_dual(cs_unscaled, rs)
    const n = length(rs)
    ratios = rs ./ cs_unscaled
    cs = cs_unscaled/n #makes a copy
    indx = sortperm(ratios, rev=true)
    out = zeros(Float64, length(indx))
    weight = 0.
    dual = -1
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
            dual = ratios[indx[ix]]
            break
        end
    end
    out, dual
end

q(cs_unscaled, rs) = q_dual(cs_unscaled, rs)[1]
function lam(zs, tau, vs, cs_unscaled)
    rs = shrink(zs, vs, tau)
    q_dual(cs_unscaled, rs)[2]
end

function dual_obj(zs, tau, vs, cs_unscaled)
    l = lam(zs, tau, vs, cs_unscaled)
    l + mean(max(0, shrink(zs, vs, tau) .- l * cs_unscaled))    
end

#used to in computing the entire tau-grid
function sbars(xs, cs, vs, f_indx)
    out = zeros(Float64, length(xs))
    for i = 1:length(out)
        out[i] = vs[i]*vs[f_indx] * (cs[i]*xs[f_indx] - cs[f_indx]*xs[i])
        out[i] /= cs[f_indx]*vs[i]*xs[i] - cs[i]*vs[f_indx]*xs[f_indx]
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
#  cs_unsc has entries order unity.  cs has entires order 1/n
## returns max_qs, tau_vals, objs
## objs[ix] is theobjective at tau_vals[ix] + eps
function best_q_tau(cs_unscaled, xs, vs, ys)
	const TOL = 1e-10
	const n = length(xs)
	ratios = xs ./ cs_unscaled
	indx = sortperm(ratios, rev=true)
	cs = cs_unscaled/n  #makes a copy

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
    svals = sbars(xs, cs, vs, frac_indx)
    max_obj = -1
    max_tau0 = 0.

    #At the beginining of every iteration qs is solution for tau0 + eps
   	iter = 0
    while tau0 > -1.
    	iter += 1
    	@assert iter < 1e7 "Max iterations exceeded"

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
    			@assert qs[frac_indx] < 1  "Degeneracy? $(qs[frac_indx])"

    			#frac_indx and s_vals stay the same
    			tau0 = svals[indx_p]
    		else #Case 2: Exiting Variable becomes fractional
    			qs[indx_p] -= (1.- qs[frac_indx]) * cs[frac_indx]/cs[indx_p]
    			@assert qs[indx_p] > 0 "Dengeracy?"

    			qs[frac_indx] = 1
    			tau0 = svals[indx_p]
    			frac_indx = indx_p
    			svals = sbars(xs, cs, vs, frac_indx)			
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
				svals = sbars(xs, cs, vs, frac_indx)
			end
		end
    end

    ##assume that solving one last time is short relative to costs
    ##have to unscale cs for this to work
    max_qs = KP.q(cs_unscaled, shrink(xs, vs, max_tau0+TOL))

    return max_qs, vals, objs 
end

#Compute the method of moment estimator
function q_MM(cs, zs, vs)
	tau_mm = 1/mean(zs.^2 - 1./vs)
	tau_mm = max(0, tau_mm)
	tau_mm, q(cs, shrink(zs, vs, tau_mm))
end

##The MLE doesn't always solve properly
##Instead of throwing, just return a -1 and the qs corresponding to tau0 = 0
function q_MLE(cs, zs, vs, max_bnd = 1e2)
	#solve the MLE using a rootfinder solver for the variance and then invert
	vs = 1./ vs
	deriv_loglik(v) = mean(zs.^2 ./ (v + vs).^2 - 1./(v + vs))

    #Check the bracket
    if deriv_loglik(0) * deriv_loglik(max_bnd) > 0
        return -1, q(cs, zs)
    else
        v0 = fzero(deriv_loglik, 0, max_bnd)
        #check if you ran up against a bound
        if v0 >= max_bnd
            return -1, q(cs, zs)
        end
    return 1/v0, q(cs, shrink(zs, vs, 1/v0))
    end
end

ideal_val(thetas, cs) = dot(q(cs, thetas), thetas) /length(thetas)

## returns the best of n^2 draws
function rand_alg(cs, xs, ys, vs, numDraws=length(xs)^2)
	const n = length(xs)
	qalg = zeros(Float64, n)
	sigs = 1./sqrt(vs)
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

#a heuristic based on ridge regression
function q_ridge(cs, xs, ys, vs; max_iter = 8)
    #compute the tau that minimizes the out-of-sample l2
    function deriv_f(tau0)
        rs = shrink(xs, vs, tau0)
        mean((ys - rs) .* rs./(vs + tau0))
    end
    #bracket tau0
    max_bnd = 1
    sgn = sign(deriv_f(0))
    iter = 0
    for iter = 1:max_iter
        if sgn * deriv_f(max_bnd) < 0
            break
        else
            max_bnd *= 2
        end
    end
    #if iteration limit reached, just return zeros
    if iter == max_iter
        return zeros(length(xs))
    end

    tau_ridge = fzero(deriv_f, 0, max_bnd)
    zs = .5 * (xs + ys)
    q(cs, shrink(zs, vs, tau_ridge))
end

function q_linreg(cs, rs, ts)
	a, b = linreg(rs, ts) 
	rhat = a + b*rs
	q(cs, rhat)
end

function perf_reg(thetas, cs, numSims)
	const n = length(thetas)
	out = zeros(Float64, n, 3)
	for iSim = 1:numSims
		Xs = randn(n) ./ sqrt(vs) + thetas
		Ys = randn(n) ./ sqrt(vs) + thetas
		q = q_linreg(cs, Xs, Ys)
		out[iSim, 1] = dot(q, thetas) / n
	end
	out
end

end #ends module
