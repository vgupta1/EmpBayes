##testing the uniform limit


module KP
using Distributions, kNN

function q_beta(cs, rs, beta)
    ratios = (beta + rs) ./ cs
    indx = sortperm(ratios, rev=true)
    out = zeros(Float64, length(indx))
    weight = 0.
    for ix = 1:length(indx)
        #check if next item fits
        if weight + cs[indx[ix]] <= 1
            out[[indx[ix]]] = 1.
            weight += cs[indx[ix]]
        else
            #take fractional part
            out[[indx[ix]]] = (1-weight)/cs[indx[ix]]
            break
        end
    end
    out
end

q(cs, rs) = q_beta(cs, rs, 0.)

# compute the distribution of rewards for various betas
# This tells us something about the oracle value
function perf_betas(thetas, cs, numSims, beta_grid = [0.])
	const n = length(thetas)

	out = zeros(Float64, numSims, length(beta_grid))
	for iSim = 1:numSims
		rs = [rand(Exponential(1/thetas[ix])) for ix = 1:n]
		for (ix, b) in enumerate(beta_grid)
			q = q_beta(cs, rs, b)
			out[iSim, ix] = dot(q, 1./thetas) / n
		end
	end
	out
end

ideal_val(thetas, cs) = dot(q_beta(cs, 1./thetas, 0), 1./thetas) /length(thetas)

function q_linreg(cs, rs, ts)
	a, b = linreg(rs, ts) 
	rhat = a + b*rs
	q(cs, rhat)
end

function q_knn(cs, rs, ts, kernel = :gaussian)
	fit = kernelregression(rs, ts, kernel) #let julia pick bandwidth
	rhat = predict(fit, rs)
	q(cs, rhat)
end

function perf_reg(thetas, cs, numSims)
	const n = length(thetas)
	out = zeros(Float64, n, 3)
	for iSim = 1:numSims
		rs = [rand(Exponential(1/thetas[ix])) for ix = 1:n]
		ts = [rand(Exponential(1/thetas[ix])) for ix = 1:n]

		q = q_linreg(cs, rs, ts)
		out[iSim, 1] = dot(q, 1./thetas) / n

		# q = q_knn(cs, rs, ts)
		# out[iSim, 2] = dot(q, 1./thetas) / n
	end
	out
end

end #ends module