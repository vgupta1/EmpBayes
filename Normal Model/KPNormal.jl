module KP
using Distributions, Roots, Optim, LinearAlgebra, MLDataUtils, Statistics

##VG Consider killing these exports
export x, best_x_tau, x_MM, x_MLE, 
      x_dual, lam, shrink, x_l2reg, 
      x_l2reg_CV, x_sure_MSE, x_stein_exact, 
      x_OR_MSE, x_stein_box

### The main workhorse
#cs_unsc (unscaled) should have entries on order unity.
#The dual return fails with value -1 if the particular run is degenerate.  
function x_dual(cs_unsc, rs)
    n = length(rs)
    ratios = rs ./ cs_unsc
    cs = cs_unsc/n #deliberately makes a copy
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
            # VG DEBUG out[[indx[ix]]] .= (1 - weight)/cs[indx[ix]]
            out[indx[ix]] = (1 - weight)/cs[indx[ix]]
            dual = ratios[indx[ix]]
            break
        end
    end
    out, dual
end

### Some generic helpers 
function shrink(muhat, vs, tau) 
    vmin = minimum(vs)
#    return @. (vmin + tau)/vmin * (vs / (vs + tau)) * muhat
    return (vmin + tau)/vmin * (vs ./ (vs .+ tau)) .* muhat
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
        out[i] = vs[i] * vs[f_indx] * (cs[i] * muhat[f_indx] - cs[f_indx] * muhat[i])
        out[i] /= cs[f_indx] * vs[i] * muhat[i] - cs[i] * vs[f_indx] * muhat[f_indx]
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
	TOL = 1e-10
	n = length(muhat)
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
            xs[indx[ix]] = (1 - total_weight)/cs[indx[ix]]
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
    			xs[indx_p] -= (1. - xs[frac_indx]) * cs[frac_indx]/cs[indx_p]
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

### Empirical Bayes Type Estimators
##EB Method of Moments
function x_MM(cs, muhat, vs)
	tau_mm = 1/mean(@. muhat^2 - 1/vs )
	tau_mm = max.(0, tau_mm)
	tau_mm, x(cs, shrink(muhat, vs, tau_mm))
end

##EB MLE 
#Note, mathematically MLE might not solve.  
#In this case, returns -1 for tau and the xs corresponding to tau0 = 0
function x_MLE(cs, muhat, vs, max_bnd = 1e2)
	#solve the MLE using a rootfinder solver for the variance and then invert
	vars = 1 ./ vs
	deriv_loglik(v) = mean(muhat.^2 ./ (v .+ vars).^2 - 1 ./(v .+ vars))

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
##Selects tau0 to minimize oracle MSE
#We can use the old shrinkage to compute tau.
function x_OR_MSE(cs, muhat, thetas, vs; max_iter = 10)   
    n = length(cs)
    function deriv_f(tau0)
        #use old scaling
        rs = muhat .* vs ./ (vs .+ tau0)
        mean( (thetas - rs) .* rs ./ (vs .+ tau0))
    end   
    #bracketting tau0
    #this code is a bit lazy.
    max_bnd = 1
    sgn = sign(deriv_f(0))
    iter_ = 0
    for iter = 1:max_iter
        if sgn * deriv_f(max_bnd) < 0
            break
        else
            max_bnd *= 2
        end
        iter_ = iter  #to get around scoping problem
    end
    #if iteration limit reached, just return zero
    tau_star = 0.
    if iter_ != max_iter
        tau_star = fzero(deriv_f, 0, max_bnd)
    end
    return x(cs, shrink(muhat, vs, tau_star)), tau_star
end

#Selects tau0 to minimize SURE of MSE
function x_sure_MSE(cs, muhat, vs)
    f_deriv(tau) = dot(vs .* (muhat.^2 * tau .- 1) .- tau, 
                            1. ./(vs .+ tau).^3)

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


## Cross-val procedure for EB class
#dat is assumed column major size(dat) = (n, S) (S for samples)
#vs is the precision of \xi, i.e, one column of dat
#muhat has precision vs * S
function x_kFoldCV(cs, vs, dat, numFolds; 
                    tau_grid = linspace(0, 10, 100))
    #divy the folds 
    S = size(dat, 2)
    n = size(dat, 1)
    train, test = kfolds(S, numFolds)
    objs_all = zeros(length(tau_grid))
    vs_sc = zeros(length(vs))

    xs_t = zeros(n)
    for i = 1:numFolds
        num_train = length(train[i])
        @. vs_sc = vs * num_train

        #form muhat  and leftover
        muhat = vec(mean(dat[:, train[i]], dims=2))
        zeta  = vec(mean(dat[:,  test[i]], dims=2))

        #eval on tau_grid pseudo-manually
        for j = 1:length(objs_all)
            xs_t[:] = x(cs, shrink(muhat, vs_sc, tau_grid[j]))
            objs_all[j] += dot(zeta, xs_t) / n /numFolds
        end
    end
    tau_hat = tau_grid[argmax(objs_all)]
    x(cs, shrink(vec(mean(dat, dims=2)), vs * S, tau_hat)), tau_grid, objs_all
end


#######
### Regularization based Methods
#########
function g_j(Gamma, lam, cj, muhat_j, vj, sqrt_vmin)
    t = muhat_j - lam * cj
    if t < 0
        return 0.
    elseif t < Gamma * sqrt_vmin / vj
        return t * vj / Gamma / sqrt_vmin 
    else
        return 1.
    end
end

#cs are unscaled (each element order unity)
#this is a bad signature.  x_out is overwritten and it returns dual variable
#Try to rationalize this and think about what is most memory efficient in code
function x_l2reg2!(cs, muhat, vs, Gamma, x_out; 
                    lam_guess::Float64 = -1., ROOT_TOL::Float64=1e-6, TOL::Float64=1e-6)
    n = length(cs)
    sqrt_vmin = sqrt(minimum(vs))
    @assert length(x_out) == n "Return vector wrong size $(n) vs. $(length(x_out))"

    #Find lambda by making knapsack tight
    function f(lam)
        mean(cs .* g_j.(Gamma, lam, cs, muhat, vs, sqrt_vmin)) - 1
    end

    #if you have lam_guess, assume gammas are close together
    #use this to guess a reasonable STARTING bracket
    if lam_guess  > 0
        lhs_prev = f(lam_guess)
        if lhs_prev > TOL
            lb = lam_guess
            ub = 2lb
        elseif lhs_prev < -TOL
            ub = lam_guess
            lb = .5*ub
        else #by luck, lam_guess is perfect
            #exit early
            x_out[:] = g_j.(Gamma, lam_guess, cs[jx], muhat[jx], vs[jx], sqrt_vmin)
            return lam_guess
        end
    else #solving from scratch
        lb, ub = 0.1, 1.        
    end

    #Now refine/correct the bracket so it is valid
    #A valid bracket satisfies f(lb) > 0, f(ub) < 0
    iter = 1
    MAX_ITER = 20
    while f(ub)*f(lb) > 0. && iter < MAX_ITER
        f(lb) <= 0 && (lb /= 2.)
        f(ub) >= 0 && (ub *= 2.)
        iter += 1       
    end

    if iter == MAX_ITER
        #maybe everything fits in the knapsack
        if f(0.) <= 0
            for jx = 1:n
                x_out[jx] = g_j(Gamma, 0., cs[jx], muhat[jx], vs[jx], sqrt_vmin)
            end
            return 0.
        else #bracketing failed
            @assert iter < MAX_ITER "Maximum bracketing iterations reached in regularized dual"
        end
    end

    #refine solution. Default to bisection if goes badly.
    lam_out = 0.
    try
        lam_out = find_zero(lam -> f(lam), [lb, ub], FalsePosition(), ftol=ROOT_TOL)
    catch
        println("False Position Failed")
        try
            lam_out = find_zero(lam -> f(lam), [lb, ub], Bisection(), ftol=ROOT_TOL)        
        catch
            println("Bisection also fails??")
            error()
        end
    end 

    x_out[:] = g_j.(Gamma, lam_out, cs, muhat, vs, sqrt_vmin)
    lam_out
end

#Selects Gamma via hold-out validation
#Use muhat/thetas to compute oracle value
#Use muhat1/muhat2 for holdout validation
#output xhat, Gamma_grid, objs
function x_l2reg_CV(cs, muhat, vs, thetas; 
                    Gamma_step = .01, Gamma_min = .1, Gamma_max = 10)
    n = length(muhat)
    Gamma_grid = collect(Gamma_min:Gamma_step:Gamma_max)

    #Exhaustive search
    best_val = -1.
    xhat = zeros(n)
    x_t  = zeros(n)
    lam_t = 0.
    objs = zeros(length(Gamma_grid))

    for (ix, Gamma) in enumerate(Gamma_grid)
        if ix > 1 #warm start
            lam_t = x_l2reg2!(cs, muhat, vs, Gamma, x_t, lam_guess=lam_t)
        else
            lam_t = x_l2reg2!(cs, muhat, vs, Gamma, x_t)
        end

        objs[ix] = dot(thetas, x_t)/n
        if objs[ix] > best_val
            best_val = objs[ix]
            xhat[:] = x_t
        end
    end     
 
    return xhat, Gamma_grid, objs
end


## Cross-val procedure for Regularization Class
#dat is column major size(dat) = (n, S)  (S for samples)
#vs is the \xi, i.e, one column of dat
#muhat has precision vs * S
function x_l2reg_kFoldCV(cs, vs, dat, numFolds; 
                    Gamma_step = .01, Gamma_min = .1, Gamma_max = 10)
    #divy the folds 
    S = size(dat, 2)
    train, test = kfolds(S, numFolds)
    Gamma_grid = collect(Gamma_min:Gamma_step:Gamma_max)
    objs_all = zeros(length(Gamma_grid))
    vs_sc = zeros(length(vs))

    for i = 1:numFolds
        #scale by size of fold, might not be equal for all flds
        num_train = length(train[i])
        vs_sc[:] = vs .* num_train

        #form muhat  and leftover
        muhat = mean(dat[:, train[i]], dims=2)
        zeta  = mean(dat[:,  test[i]], dims=2)

        #eval on GammaGrid using above workhorse
        objs = x_l2reg_CV(cs, muhat, vs_sc, zeta; 
                    Gamma_step=Gamma_step, Gamma_min=Gamma_min, Gamma_max=Gamma_max)[3]
        objs_all[:] += objs ./ numFolds
    end
    #avg and choose Gammahat
    #Crucially gamma_grid same between runs and as above
    Gammahat = Gamma_grid[argmax(objs_all)]
    xs = zeros(size(dat, 1))
    x_l2reg2!(cs, mean(dat, dims=2), vs * S, Gammahat, xs)

    xs, Gamma_grid, objs_all
end


### EB Optimization Approaches

##Dirac Expansion
#Evaluation of Stein term by approximating lambda approx
#Serves as an upperbound on choice of bandwidth
#cs are unscaled
function dirac_bias(vs, tau, muhat, lam, cs_unsc, thetas)
    out = 0.
    n = length(vs)
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
    n = length(vs)
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
    n = length(vs)
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
    n = length(vs)
    tau_grid = collect(0.:tau_step:tau_max)
    objs = fill(0., size(tau_grid))
    obj_best = -Inf
    tau_best = -1.
    for ix = 1:length(tau_grid)
        xs, lam = x_dual(cs_unsc, shrink(muhat, vs, tau_grid[ix]))
        #compute the value corresponding to evaluating the dirac
        #approximates lambda by the finite dual value...
        objs[ix] = dot(muhat, xs)/n - 
                approx_bias(vs, tau_grid[ix], shrink(muhat, vs, tau_grid[ix]), 
                                lam, cs_unsc, h)

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
    n = length(vs)
    ej = zeros(n)
    for jx = 1:n
        if jx > 1
            ej[jx-1] = 0.
        end
        ej[jx] = h / sqrt(vs[jx])

        xs_p = x(cs, shrink(muhat .+ ej, vs, tau))[jx]
        xs_m = x(cs, shrink(muhat .- ej, vs, tau))[jx]

        out += (xs_p - xs_m) / 2 / h / sqrt(vs[jx])
    end
    out/n
end

# a smarter way to do this would leverage the parametric programming
function x_stein_primal(cs_unsc, muhat, vs, h; tau_step = .01, tau_max = 5.)
    tau_grid = collect(0.:tau_step:tau_max)
    objs = zeros(length(tau_grid))
    obj_best = -Inf
    tau_best = -1.
    n =length(vs)
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
    n = length(vs)
    sqrt_vmin = sqrt(minimum(vs))
    for jx = 1:n
        if 0 <= muhats[jx] - lam * cs[jx] <= Gamma * sqrt_vmin/vs[jx]
            out += 1
        end
    end
    out / Gamma / sqrt_vmin / n
end


function x_stein_reg(cs_unsc, muhat, vs; Gamma_step = .01, Gamma_min = .1, Gamma_max = 10)
    Gamma_grid = collect(Gamma_min:Gamma_step:Gamma_max)

    #preallocated for speed
    n =length(vs)
    objs = fill(0., size(Gamma_grid))
    xs = zeros(n)

    #updated as we find better values
    obj_best, Gamma_best, lam = -Inf, -1., -1.

    for ix = 1:length(Gamma_grid)
        #use warm-start information
        if ix > 1
            lam = x_l2reg2!(cs_unsc, muhat, vs, Gamma_grid[ix], xs, lam_guess=lam)
        else
            lam = x_l2reg2!(cs_unsc, muhat, vs, Gamma_grid[ix], xs)
        end
        objs[ix] = dot(muhat, xs)/n - reg_bias(cs_unsc, vs, muhat, lam, Gamma_grid[ix])

        if objs[ix] > obj_best
            Gamma_best = Gamma_grid[ix]
            obj_best = objs[ix]
        end
    end    
    Gamma_best >= Gamma_max-1e-10 && println("Gamma Grid too small ", Gamma_min, " ", Gamma_max)
    lam = x_l2reg2!(cs_unsc, muhat, vs, Gamma_best, xs)
    return xs, Gamma_grid, objs
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

    #Exploits fact that Gamma_grid is identical.
    objs = .5 .* (objs1 .+ objs2)
    Gamma_best = Gamma_grid[argmax(objs)]
    Gamma_best >= Gamma_max-1e-10 && println("Gamma Grid too small ", Gamma_min, " ", Gamma_max)

    lam = x_l2reg2!(cs_unsc, .5 .* (muhat1 .+ muhat2), vs, Gamma_best, xs1)
    return xs1, Gamma_grid, objs
end

#Solves the robust problem using FW
#assumes r is radius of ellipse, i.e. rob counterpart has r/n * norm(xs)_vs
function x_robFW(cs, muhat, vs, r; MAX_ITER = 100, TOL=1e-6)
    iter = 1
    n = length(muhat)
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

        if abs(-f(step_star) - prev_value) <= TOL
            break
        else
            prev_value = -f(step_star)
        end

        xp += step_star * (d - xp)
        iter += 1
    end

    @assert iter < MAX_ITER "Maximum Iterations reached in FW"
    return xp
end




end #ends module
