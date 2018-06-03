module KP
using Distributions, Roots, Optim

export x, best_x_tau, x_MM, x_MLE, 
      x_dual, lam, shrink, x_l2reg, 
      x_l2reg_CV, x_sure_MSE, x_stein_exact, 
      x_OR_MSE, x_stein_box, x_l2reg_warm

### The main workhorse
#cs_unsc (unscaled) should have entries on order unity.
#The dual return fails with value -1 if the particular run is degenerate.  
function x_dual(cs_unsc, rs)
    const n = length(rs)
    ratios = rs ./ cs_unsc
    cs = cs_unsc/n #makes a copy
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
            out[[indx[ix]]] = (1 - weight)/cs[indx[ix]]
            dual = ratios[indx[ix]]
            break
        end
    end
    out, dual
end

### Some generic helpers 
function shrink(muhat, vs, tau) 
    vmin = minimum(vs)
    return (vmin + tau)/vmin * (vs ./ (vs + tau)) .* muhat
end

x(cs_unsc, rs) = x_dual(cs_unsc, rs)[1]

function lam(muhat, tau, vs, cs_unsc)
    rs = shrink(muhat, vs, tau)
    x_dual(cs_unsc, rs)[2]
end

function dual_obj(muhat, tau, vs, cs_unsc)
    l = lam(muhat, tau, vs, cs_unsc)
    l + mean(max(0, shrink(muhat, vs, tau) .- l * cs_unsc))    
end


## Computing the whole grid 
# helper used to identify next crossing point
# not to be exposed
function sbars(muhat, cs, vs, f_indx)
    out = zeros(Float64, length(muhat))
    for i = 1:length(out)
        out[i] = vs[i]*vs[f_indx] * (cs[i]*muhat[f_indx] - cs[f_indx]*muhat[i])
        out[i] /= cs[f_indx]*vs[i]*muhat[i] - cs[i]*vs[f_indx]*muhat[f_indx]
    end   
    out
end

#Min_indx above a threshold
#returns -1 if there are no values above thresh
function indmin_thresh(vals, thresh)
	min_indx = -1
	min_val = Inf
	for ix = 1:length(vals)
		if thresh < vals[ix] < min_val
			min_indx = ix
			min_val = vals[ix]
		end
	end
	min_indx
end

# searches over bfs to figure out the best value of tau0 
# cs_unsc has entries order unity.  cs has entires order 1/n
# returns max_xs, tau_vals, objs
# objs[ix] is the objective at tau_vals[ix] + eps
# For oracle, use muhat, thetas
# For Hold_out, use muhat1, muhat2
function best_x_tau(cs_unsc, muhat, vs, thetas)
	const TOL = 1e-10
	const n = length(muhat)
	ratios = muhat ./ cs_unsc
	indx = sortperm(ratios, rev=true)
	cs = cs_unsc/n  #makes a copy

	vals = Float64[]
	objs = Float64[]

	#first determine tau0 = 0 solution
	xs = zeros(n)
	total_weight = 0.
	frac_indx = -1
    for ix = 1:n
    	#check if the item has positive value, else skip
    	if muhat[indx[ix]] <=0. 
    		continue
    	end

        #check if next item fits entirely
        if total_weight + cs[indx[ix]] <= 1
            xs[indx[ix]] = 1.
            total_weight += cs[indx[ix]]
        else
            #take fractional part
            xs[[indx[ix]]] = (1.-total_weight)/cs[indx[ix]]
            frac_indx = indx[ix]
            total_weight = 1.
            break
        end
    end

    ##check for a trivial solution
    if frac_indx == -1 || total_weight < 1
    	println("Trivial Problem Set-up:  All items fit.")
    	println("total_weight:\t $total_weight")
    	println("xs:\n", xs)
    	return xs, [0.], [dot(xs, thetas)/n]
    end

    #Now search across all possible tau0, constructing grid as you go
	tau0 = 0.  #always represents current value
    svals = sbars(muhat, cs, vs, frac_indx)
    max_obj = -1
    max_tau0 = 0.

    #At the beginining of every iteration xs is solution for tau0 + eps
   	iter = 0
    while tau0 > -1.
    	iter += 1
    	@assert iter < 1e7 "Max iterations exceeded"

	    #assess solution and record value
	    obj = dot(thetas, xs)/n
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
    	if xs[indx_p] >= 1-TOL
    		#Case 1: Exiting Variable becomes zero
    		if cs[indx_p]/cs[frac_indx] <= 1 - xs[frac_indx]
    			xs[indx_p] = 0
    			xs[frac_indx] += cs[indx_p]/cs[frac_indx]
    			@assert xs[frac_indx] < 1  "Degeneracy? $(xs[frac_indx])"

    			#frac_indx and s_vals stay the same
    			tau0 = svals[indx_p]
    		else #Case 2: Exiting Variable becomes fractional
    			xs[indx_p] -= (1.- xs[frac_indx]) * cs[frac_indx]/cs[indx_p]
    			@assert xs[indx_p] > 0 "Dengeracy?"

    			xs[frac_indx] = 1
    			tau0 = svals[indx_p]
    			frac_indx = indx_p
    			svals = sbars(muhat, cs, vs, frac_indx)			
    		end
    	else  #xs[indx_p] == 0
	    	frac_weight = xs[frac_indx] * cs[frac_indx]
    		#Case 1: Entering Variable becomes 1
	    	if frac_weight > cs[indx_p]  
	    		xs[indx_p] = 1.
	    		xs[frac_indx] -= cs[indx_p]/cs[frac_indx]
	    		@assert xs[frac_indx] > 0 "Dengeracy?"
	    		
	    		#frac_indx and svals are unchanged
	    		tau0 = svals[indx_p]
	    	else #Case 2: Entering Variable becomes fractional
				xs[indx_p] = frac_weight / cs[indx_p]
				xs[frac_indx] = 0.
				tau0 = svals[indx_p]	
				frac_indx = indx_p
				svals = sbars(muhat, cs, vs, frac_indx)
			end
		end
    end

    ##assume that solving one last time is short relative to costs
    indmax = findfirst(vals .>= max_tau0)
    if indmax == length(vals)
        max_xs = x(cs_unsc, shrink(muhat, vs, 1.1*max_tau0))
    else
        max_xs = x(cs_unsc, shrink(muhat, vs, .5*max_tau0 + .5*vals[indmax + 1]))
    end

    return max_xs, vals, objs 
end

#Naive grid search for best value of $\tau$
#cs_unsc has entries order unity.  cs has entires order 1/n
#returns max_xs, tau_vals, objs
#objs[ix] is the objective at tau_vals[ix] + eps
#For oracle, use muhat, thetas
#For Hold_out, use muhat1, muhat2
function best_x_tau2(cs_unsc, muhat, vs, thetas; tau_step = .01, tau_max = 5.)
    const n = length(vs)
    tau_grid = collect(0.:tau_step:tau_max)
    objs = zeros(length(tau_grid))
    obj_best = -Inf
    tau_best = -1.
    for ix = 1:length(tau_grid)
        xs, lam = x_dual(cs_unsc, shrink(muhat, vs, tau_grid[ix]))
        #compute the value corresponding to evaluating the dirac
        objs[ix] = dot(thetas, xs)/n

        if objs[ix] > obj_best
            tau_best = tau_grid[ix]
            obj_best = objs[ix]
        end
    end
    #solving once more is fast compared to rest of algorithm
    return x(cs_unsc, shrink(muhat, vs, tau_best)), tau_grid, objs
end

### Empirical Bayes Type Estimators
##EB Method of Moments
function x_MM(cs, muhat, vs)
	tau_mm = 1/mean(muhat.^2 - 1./vs)
	tau_mm = max(0, tau_mm)
	tau_mm, x(cs, shrink(muhat, vs, tau_mm))
end

##EB MLE 
#Note, mathematically MLE might not solve.  
#In this case, returns -1 for tau and the xs corresponding to tau0 = 0
function x_MLE(cs, muhat, vs, max_bnd = 1e2)
	#solve the MLE using a rootfinder solver for the variance and then invert
	vars = 1./ vs
	deriv_loglik(v) = mean(muhat.^2 ./ (v + vars).^2 - 1./(v + vars))

    #Check the bracket
    if deriv_loglik(0) * deriv_loglik(max_bnd) > 0
        return -1, x(cs, muhat)
    else
        v0 = fzero(deriv_loglik, 0, max_bnd)
        #check if you ran up against a bound
        if v0 >= max_bnd
            return -1, x(cs, muhat)
        end
    return 1/v0, x(cs, shrink(muhat, vs, 1/v0))
    end
end

### Other Estimate then Optimize Heuristics

#Selects tau0 to minimize out-of-sample MSE
#out of sample mse estimated via hold-out validation
#This approximates ridge regression
#Run on muhat1, muhat2 to get hold-out
function x_HO_MSE(cs, muhat1, muhat2, vs; max_iter = 5)   
    function deriv_f(tau0)
        rs = shrink(muhat1, vs, tau0)
        mean((muhat2 - rs) .* rs./(vs + tau0))
    end
    #bracketting tau0
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
    #if iteration limit reached, just return zero
    tau_star = 0.
    if iter != max_iter
        tau_star = fzero(deriv_f, 0, max_bnd)
    end
    muhat = .5 * (muhat1 + muhat2)
    return x(cs, shrink(muhat1, vs, tau_star)), tau_star
end

##Selects tau0 to minimize oracle MSE
#We can use the old shrinkage to compute tau.
function x_OR_MSE(cs, muhat, thetas, vs; max_iter = 10)   
    const n = length(cs)
    function deriv_f(tau0)
        #use old scaling
        rs = muhat .* vs ./ (vs + tau0)
        mean( (thetas - rs) .* rs ./ (vs + tau0))
    end   
    #bracketting tau0
    #this code is a bit lazy.
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
    #if iteration limit reached, just return zero
    tau_star = 0.
    if iter != max_iter
        tau_star = fzero(deriv_f, 0, max_bnd)
    end
    return x(cs, shrink(muhat, vs, tau_star)), tau_star
end

#Selects tau0 to minimize SURE of MSE
function x_sure_MSE(cs, muhat, vs)
    f_deriv(tau) = dot(vs.*(muhat.^2 * tau - 1) - tau, 1./(vs + tau).^3)

    #bracket for ub
    lb, ub = 0., 1.
    iter = 1
    MAX_ITER = 20
    while f_deriv(ub) < 0 && iter < MAX_ITER
        lb = ub
        ub *= 2
        iter += 1
    end
    if iter == MAX_ITER
        return x(cs, muhat), -1 #Returnthe SAA with an error flat
    end

    @assert iter < MAX_ITER "max iterations reached"
    tau_star = fzero(f_deriv, lb, ub)
    rs = shrink(muhat, vs, tau_star)

    return x_dual(cs, rs)[1], tau_star
end

### Regularization based Methods
#helper function.  not to be exposed
#cs are unscaled (each element order unity)
function g_j(Gamma, lam, cj, muhat_j, vj, sqrt_vmin)
    const t = muhat_j - lam * cj
    out = max(t, 0) - max(t - Gamma * sqrt_vmin/vj, 0)
    out *= vj/Gamma/sqrt_vmin 
    return out
end

#cs are unscaled (each element order unity)
#returns xs, lam
function x_l2reg(cs, muhat, vs, Gamma)
    const n = length(cs)

    #Find lambda by making knapsack tight
    sqrt_vmin = sqrt(minimum(vs))
    function f(lam)
        out = 0.
        for jx = 1:n
            out += cs[jx] * g_j(Gamma, lam, cs[jx], muhat[jx], vs[jx], sqrt_vmin)
        end
        out / n
    end

    #Check if we can trivially fit all items
    if f(0) <= 1
        return [g_j(Gamma, 0., cs[jx], muhat[jx], vs[jx], sqrt_vmin) for jx = 1:n], 0.
    else
        #bracketting
        lb, ub = 0., 1.
        iter = 1
        MAX_ITER = 20
        while f(ub) > 1
            @assert iter < MAX_ITER "Maximum iterations reached in regularized dual"
            lb = ub
            ub *= 2
            iter += 1
        end
        lamstar = fzero(lam -> f(lam) - 1, [lb, ub])
        return [g_j(Gamma, lamstar, cs[jx], muhat[jx], vs[jx], sqrt_vmin) for jx =1:n], lamstar
    end
end

#Solves the Regularization problem from a warm-start
#If heuristic fails, reverts to solution above
function x_l2reg_warm(cs, muhat, vs, Gamma; 
                        lambda_0 = -1., Gamma_0 = Inf)
    const n = length(cs)
    #Find lambda by making knapsack tight
    sqrt_vmin = sqrt(minimum(vs))
    function f(lam)
        out = 0.
        for jx = 1:n
            out += cs[jx] * g_j(Gamma, lam, cs[jx], muhat[jx], vs[jx], sqrt_vmin)
        end
        out / n
    end

    #Check if we can trivially fit all items
    if f(0) <= 1
        return [g_j(Gamma, 0., cs[jx], muhat[jx], vs[jx], sqrt_vmin) for jx = 1:n], 0.
    end

    #Use derivative to approximate a sol if you have warmstart
    if lambda_0 > 0
        deriv = 0.
        for jx = 1:n
            deriv += cs[jx]^2 * vs[jx] / Gamma / sqrt_vmin
        end
        deriv /= n        
        lambda_approx = max(lambda_0 + (1. - f(lambda_0))/deriv, 0)
    else
        lambda_approx = 0
    end

    #Now look for a good bracket.  
    #Should have f(lb) > 1 >= f(ub)
    lb, ub = 0., 1.
    if (Gamma_0 < Gamma) && (lambda_0 > 0)
        ub = lambda_0  #since decreasing
        if abs(f(ub) - 1.) <= 1e-10 #Just in case you nail it
            return [g_j(Gamma, ub, cs[jx], muhat[jx], vs[jx], sqrt_vmin) for jx =1:n], ub
        end            
        if f(lambda_approx) > 1
            lb = lambda_approx
        else
            ub = min(lambda_approx, ub)
        end
    else #either Gamma_0 is bigger or wasn't set
        lb = max(0, lambda_0) #since decreasing
        if abs(f(lb) - 1.) <= 1e-10  #Just in case you nail it
            return [g_j(Gamma, lb, cs[jx], muhat[jx], vs[jx], sqrt_vmin) for jx =1:n], lb
        end

        if f(lambda_approx) <= 1
            ub = lambda_approx
        else
            lb = max(lb, lambda_approx)
            ub = 2lb + 1
            iter = 1
            MAX_ITER = 20
            while f(ub) > 1
                @assert iter < MAX_ITER "Maximum iterations reached in regularized dual"
                lb = ub
                ub *= 2
                iter += 1
            end
        end
    end
    if (f(ub) > 1 || f(lb) <= 1)
        println("[$lb, $ub]")
        println("Initial $(lambda_0) $Gamma_0")
        println("Gamma $Gamma")
        println("f(ub) $(f(ub)) $(f(ub) <=1)")
        println("f(lb) $(f(lb)) $(f(lb) > 1)")
        println("lambda_approx $lambda_approx")
    end

    @assert f(lb) > 1  "lb is incorrect in bracketing [$lb, $ub]"
    @assert f(ub) <= 1  "ub is incorrect in bracketing [$lb, $ub]"

    lamstar = fzero(lam -> f(lam) - 1, [lb, ub])
    return [g_j(Gamma, lamstar, cs[jx], muhat[jx], vs[jx], sqrt_vmin) for jx =1:n], lamstar
end


#Selects Gamma via hold-out validation
#Test against muhat/thetas for oracle, use muhat1/muhat2 for holdout
#output xhat, Gamma_grid, objs
function x_l2reg_CV(cs, muhat, vs, thetas; Gamma_step = .01, Gamma_min = .1, 
                                            Gamma_max = 10)
    const n = length(muhat)
    Gamma_grid = collect(Gamma_min:Gamma_step:Gamma_max)

    #an exhaustive search
    best_val = -1.
    xhat = zeros(n)
    x_t = zeros(n)
    Gamma_hat = -1;
    objs = zeros(length(Gamma_grid))

    for (ix, Gamma) in enumerate(Gamma_grid)
        x_t[:] = x_l2reg(cs, muhat, vs, Gamma)[1]
        objs[ix] = dot(thetas, x_t)/n
        if objs[ix] > best_val
            best_val = objs[ix]
            best_Gamma = Gamma
            xhat[:] = x_t
        end
    end     
    return xhat, Gamma_grid, objs
end


## Uses warm-start information to solve and clever bounds on Gamma OR
function x_l2reg_CV_warm(cs, muhat, vs, thetas; 
                    Gamma_step = .01, Gamma_min = .1, Gamma_max = 10)
    const n = length(muhat)
    Gamma_grid = collect(Gamma_min:Gamma_step:Gamma_max)

    #an exhaustive search
    #updated as we find better values
    best_val, Gamma_hat, lam = -1., -1, -1
    #preallocated for speed
    xhat = zeros(n)
    x_t = zeros(n)
    objs = zeros(Gamma_grid)

    #used in trying to end early    
    sq_norm_u_vinv = mean(thetas.^2 .* vs)
    thresh, len_grid = -Inf, length(Gamma_grid)

    for (ix, Gamma) in enumerate(Gamma_grid)
        if ix == 1
            x_t[:], lam = x_l2reg_warm(cs, muhat, vs, Gamma)
        else
            x_t[:], lam = x_l2reg_warm(cs, muhat, vs, Gamma, lambda_0 = lam, Gamma_0=Gamma_grid[ix - 1])
        end
       #check if we can terminate early
        if mean(x_t[:].^2 ./ vs) < thresh
            len_grid = ix
            break  #all future values of x will be too small
        end

        objs[ix] = dot(thetas, x_t)/n
        if objs[ix] > best_val
            best_val = objs[ix]
            best_Gamma = Gamma
            xhat[:] = x_t
            thresh = best_val^2 / sq_norm_u_vinv
        end
    end     
    Gamma_grid[indmax(objs)] > Gamma_max && println("Gamma_Grid too small")
    return xhat, Gamma_grid[1:len_grid], objs[1:len_grid]
end


### EB Optimization Approaches

##Dirac Expansion
#Evaluation of Stein term by approximating lambda approx
#Serves as an upperbound on choice of bandwidth
#cs are unscaled
function dirac_bias(vs, tau, muhat, lam, cs_unsc, thetas)
    out = 0.
    const n = length(vs)
    vmin = minimum(vs)
    for jx = 1:n
        rj = (vmin + tau)/vmin * vs[jx]/(vs[jx] + tau) * thetas[jx]
        arg = (lam * cs_unsc[jx] - rj) 
        arg *= (vs[jx] + tau) * vmin / (vmin + tau) / sqrt(vs[jx])
        out += exp(-.5 * arg^2) / sqrt(2*pi * vs[jx])
    end
    return out/n
end

##Stein approach using an exact dirac calculation
#only approximations are lam(t,z) -> lam(t) and ULLN
#Serves as an upperbound on performance of a kernel bandwith
function x_stein_exact(cs_unsc, muhat, vs, thetas; tau_step = .01, tau_max = 5.)
    const n = length(vs)
    tau_grid = collect(0.:tau_step:tau_max)
    objs = zeros(length(tau_grid))
    obj_best = -Inf
    tau_best = -1.
    for ix = 1:length(tau_grid)
        xs, lam = x_dual(cs_unsc, shrink(muhat, vs, tau_grid[ix]))
        #compute the value corresponding to evaluating the dirac
        objs[ix] = dot(muhat, xs)/n - dirac_bias(vs, tau_grid[ix], muhat, lam, cs_unsc, thetas)

        if objs[ix] > obj_best
            tau_best = tau_grid[ix]
            obj_best = objs[ix]
        end
    end
    #solving once more is fast compared to rest of algorithm
    return x(cs_unsc, shrink(muhat, vs, tau_best)), tau_grid, objs
end

#Approximate the dirac with an box kernel
#rs is the shrunken value 
function approx_bias(vs, tau, rs, lam, cs, h)
    out = 0.
    const n = length(vs)
    vmin = minimum(vs)
    hj = h
    for jx = 1:n
        hj = h * (vmin + tau)/vmin * sqrt(vs[jx])/(vs[jx] + tau)
        if abs(rs[jx] - lam * cs[jx]) <= hj
            out += 1/sqrt(vs[jx]) / h / 2
        end
    end
    out/n
end

#Stein formula using an kernel approximation for the dirac
function x_stein_box(cs_unsc, muhat, vs, h; 
                     tau_step = .01, tau_max = 5.)
    const n = length(vs)
    tau_grid = collect(0.:tau_step:tau_max)
    objs = zeros(tau_grid)
    obj_best = -Inf
    tau_best = -1.
    for ix = 1:length(tau_grid)
        xs, lam = x_dual(cs_unsc, shrink(muhat, vs, tau_grid[ix]))
        #compute the value corresponding to evaluating the dirac
        #approximates lambda by the finite dual value...
        objs[ix] = dot(muhat, xs)/n - approx_bias(vs, tau_grid[ix], muhat, lam, cs_unsc, h)

        if objs[ix] > obj_best
            tau_best = tau_grid[ix]
            obj_best = objs[ix]
        end
    end
    
    return x(cs_unsc, shrink(muhat, vs, tau_best)), tau_grid, objs
end

##Approximating Stein by First Order (primal) difference
#Lazy direct way.
function primal_approx_bias(vs, tau, muhat, lam, cs, h)
    out = 0.
    const n = length(vs)
    ej = zeros(n)
    for jx = 1:n
        if jx > 1
            ej[jx-1] = 0.
        end
        ej[jx] = .5h

        xs_p = x(cs, shrink(muhat + ej, vs, tau))[jx]
        xs_m = x(cs, shrink(muhat - ej, vs, tau))[jx]

        out += (xs_p - xs_m)/h/vs[jx]
    end
    out/n
end

# a smarter way to do this would leverage the parametric programming
function x_stein_primal(cs_unsc, muhat, vs, h; tau_step = .01, tau_max = 5.)
    tau_grid = collect(0.:tau_step:tau_max)
    objs = zeros(length(tau_grid))
    obj_best = -Inf
    tau_best = -1.
    const n =length(vs)
    for ix = 1:length(tau_grid)
        xs, lam = x_dual(cs_unsc, shrink(muhat, vs, tau_grid[ix]))
        #compute the value corresponding to evaluating the dirac
        #approximates lambda by the finite dual value...
        objs[ix] = dot(muhat, xs)/n - primal_approx_bias(vs, tau_grid[ix], muhat, lam, cs_unsc, h)

        if objs[ix] > obj_best
            tau_best = tau_grid[ix]
            obj_best = objs[ix]
        end
    end
    
    return x(cs_unsc, shrink(muhat, vs, tau_best)), tau_grid, objs
end

##Regularization bias via stein's lemma
function reg_bias(cs, vs, muhats, lam, Gamma)
    out = 0.
    const n = length(vs)
    const sqrt_vmin = sqrt(minimum(vs))
    for jx = 1:n
        if 0 <= muhats[jx] - lam * cs[jx] <= Gamma * sqrt_vmin/vs[jx]
            out += 1
        end
    end
    out / Gamma / sqrt_vmin / n
end

#A smarter way to do this would leverage smoothness in gamma to update solutions
#rather than resolve
function x_stein_reg(cs_unsc, muhat, vs; Gamma_step = .01, Gamma_min = .1, Gamma_max = 10)
    Gamma_grid = collect(Gamma_min:Gamma_step:Gamma_max)

    #preallocated for speed
    const n =length(vs)
    objs = zeros(Gamma_grid)
    xs = zeros(n)

    #updated as we find better values
    obj_best, Gamma_best, lam = -Inf, -1., -1.

    for ix = 1:length(Gamma_grid)
        #use warm-start information
        if ix == 1
            xs[:], lam = x_l2reg_warm(cs_unsc, muhat, vs, Gamma_grid[ix])
        else
            xs[:], lam = x_l2reg_warm(cs_unsc, muhat, vs, Gamma_grid[ix], lambda_0=lam, Gamma_0=Gamma_grid[ix-1])
        end
        objs[ix] = dot(muhat, xs)/n - reg_bias(cs_unsc, vs, muhat, lam, Gamma_grid[ix])

        if objs[ix] > obj_best
            Gamma_best = Gamma_grid[ix]
            obj_best = objs[ix]
        end
    end    
    Gamma_best >= Gamma_max-1e-10 && println("Gamma Grid too small $Gamma_min $Gamma_max")
    return x_l2reg(cs_unsc, muhat, vs, Gamma_best)[1], Gamma_grid, objs
end

function x_LOO_reg(cs_unsc, muhat1, muhat2, vs; 
                    Gamma_step = .01, Gamma_min = .1, Gamma_max = 10)
    #VG What is the nice way to do this?
    xs1, Gamma_grid, objs1 = x_l2reg_CV(cs_unsc, muhat1, vs/2, muhat2, 
                        Gamma_step=Gamma_step, Gamma_min=Gamma_min, 
                        Gamma_max=Gamma_max)
    xs2, Gamma_grid, objs2 = x_l2reg_CV(cs_unsc, muhat2, vs/2, muhat1, 
                        Gamma_step=Gamma_step, Gamma_min=Gamma_min, 
                        Gamma_max=Gamma_max)

    #Depends on fact that Gamma_grid is identical.
    objs = .5 .* (objs1 + objs2)
    Gamma_best = Gamma_grid[indmax(objs)]
    Gamma_best >= Gamma_max-1e-10 && println("Gamma Grid too small $Gamma_min $Gamma_max")

    return x_l2reg(cs_unsc, .5 .* (muhat1 + muhat2), vs, Gamma_best)[1], Gamma_grid, objs
end

#Solves the ellipsoidal robust problem with radius r
# Algorithm seeks the equivalent "Gamma" for the regularized problem that matches the KKT conditions
# Corrects gamma_min, gamma_max if not a bracket
function x_rob(cs, muhat, vs, r; gamma_min=.01, gamma_max =100.)
    const sqrt_vmin = sqrt(minimum(vs))
    function f(Gamma)
        xs, lam = KP.x_l2reg(cs, muhat, vs, Gamma)
        return Gamma * sqrt_vmin * sqrt(sum(xs.^2 ./ vs)) - r
    end

    #Search for a bracket
    #assume this function is increasing.
    MAX_ITER = 100
    iter = 1
    if f(gamma_min) > 0 
        println("Gamma min was too large")
        while f(gamma_min) > 0
            gamma_max = gamma_min
            gamma_min /= 2
            iter += 1
            @assert iter < MAX_ITER "Maximum iterations reached in bracketing robust"
        end
    elseif f(gamma_max) < 0
        println("Gamma_max was too small")
        while f(gamma_max) < 0
            gamma_min = gamma_max
            gamma_max *= 2
            iter += 1
            @assert iter < MAX_ITER "Maximum iterations reached in bracketing robust"
        end
    end

    @assert f(gamma_min) * f(gamma_max) < 0 "Gamma_min, Gamma_max not a bracket"
    Gammastar = fzero(f, gamma_min, gamma_max)
    return KP.x_l2reg(cs, muhat, vs, Gammastar)
end

#Solves the robust problem using FW
#assumes r is radius of ellipse, i.e. rob counterpart has r/n * norm(xs)_vs
function x_robFW(cs, muhat, vs, r; MAX_ITER = 100, TOL=1e-6)
    iter = 1
    const n = length(muhat)
    sqrt_vs = sqrt.(vs)
    function grad!(xs, g)
        g[:] = muhat/n - r/norm(xs ./ sqrt_vs)/n * xs ./ vs 
    end
    g = zeros(n)
    d = zeros(n)
    xp = x(cs, muhat)
    prev_value = -Inf
    while iter < MAX_ITER
        #solve the subproblem to find a direction.
        grad!(xp, g)
        d[:] = x(cs, g)

        #optimize the step
        #notice the negatives because we are maximizing. 
        f(step) = -dot(muhat, xp + step * (d - xp))/n +
                    r/n * norm((xp + step * (d - xp))./ sqrt_vs) 
        step_star = Optim.minimizer(optimize(f, 0, 1))

        #check to stop 
        #println("iter:\t $(iter) diff:\t $(abs(-f(step_star) - prev_value))")

        if abs(-f(step_star) - prev_value) <= TOL
            break
        else
            prev_value = -f(step_star)
        end

        xp += step_star * (d - xp)
        iter += 1
    end
    #Dumping in case of failure.
    if iter > MAX_ITER
        f = open("BadRobust_Example_$(r).csv", "w")
        writecsv(f, [cs muhat vs ])
        close(f)
    end

    @assert iter < MAX_ITER "Maximum Iterations reached in FW"
    return xp
end



end #ends module
