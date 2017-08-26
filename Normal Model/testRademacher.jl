include("KPNormal.jl")
using Distributions, KP

##a little experiment to understand rademacher conv
function radTests(file_out, numRuns, n_max, seed=8675309)
	srand(seed)	
	#open and set-up an output file
	f = open("$(file_out)_$(seed).csv", "w")
	writecsv(f, ["Run" "n" "Diff" "tauOR"])
	for iRun =1:numRuns
		#simulate some data
		thetas = rand(n_max)
		thetas = sign(thetas - .5)
		vs = rand(n_max) * 3
		cs = rand(n_max)  * 2 * 10;
		zs = thetas + randn(n_max) ./ sqrt(vs);

		#iterate over the n_grid

		for n = round(Int64, 2.^(5:log2(n_max)))
		    #solve for tauORZ, and compute the requisite difference.
		    qOR, tau_grid_full, objs  = KP.best_x_tau(cs[1:n], zs[1:n], vs[1:n], thetas[1:n])
		    tauOR = tau_grid_full[indmax(objs)]
		    diff = dot(thetas[1:n] - shrink(zs[1:n], vs[1:n], tauOR), qOR)/n

		    writecsv(f, [iRun n diff tauOR])
		end
		flush(f)
	end
	close(f)
end

radTests("temp", 3, 64)
