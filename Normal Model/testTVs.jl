include("testCLTHarness.jl")
using KernelDensity, QuadGK

### Helper function for TVs
function TV(dist_a, dist_b, a, b)
    interp_a = InterpKDE(dist_a)
    interp_b = InterpKDE(dist_b)
    quadgk( x-> abs(pdf(interp_a, x) - pdf(interp_b, x)) * .5, a, b)
end

######################
####
# Across an S-grid, for a fixed n, 
#evaluate the average TV for the robustness to non-normality test
# f file to write to 
# numSamples -- number of samples to construct KDE estimate
# o  -- CLTExp() pre-configured 
function test_computeTVs(numSamples, o)
	#gen all the data upfront.  Not memory efficient but whatever.
	n = size(o.data, 1)
	dat = zeros(n, numSamples)
	muhat = zeros(n)
	for i = 1:numSamples
		sim!(o, muhat, true)  #constant precision version	
		dat[:, i] .= muhat
	end

	#one standard normal to be used for everything
	ref_dat = randn(numSamples)
	kde_normal = kde_lscv(ref_dat)

	#now iteratate through across n
	dat_j = zeros(numSamples)
	out = zeros(n)
	for j = 1:n
		dat_j[:] .= vec( (dat[j, :] .- o.thetas[j]) .* sqrt(o.vs[j]) )

		#do integrations across mu += 5 stdevs
		sigma = sqrt(mean(dat_j.^2))

		#build density obj and interpolator
		kde_j = kde_lscv(dat_j)
		out[j] = TV(kde_normal, kde_j, -Inf, Inf)[1]	
#		out[j] = TV(kde_normal, kde_j, min(-7, -7*sigma), max(7, 7 * sigma))[1]	

	end
	return mean(out)
end

##Here's the actual test function
function test_computeTVs(numSamples, s_path_out, n, S_grid, distribution; seed = 8675309)
	param_path = "./Results/param_portExp_mtn2.csv"
	Random.seed!(seed)
	dat, header = readdlm(param_path, ',', header=true)

	@assert n <= size(dat, 1) "Param file too short for n"
	cs = vec(dat[1:n, 3])
	thetas = vec(dat[1:n, 1])
	vs = vec(dat[1:n, 2])

	dist = get_dist(distribution)
	out = zeros(length(S_grid))
	for (i, S) in enumerate(S_grid)
		o = CLTExp(cs, thetas, vs, zeros(n, S), dist)
		out[i] = test_computeTVs(numSamples, o)
	end

	file_name = "$(s_path_out)_$(distribution)_$(seed).csv"
	file_out = open(file_name, "w")
	writedlm(file_out, ["S" "TV"], ',')
	writedlm(file_out, [vec(S_grid) out])
	close(fil_out)

	return file_name
end